#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <mudock/chem/autodock_parameters.hpp>
#include <mudock/chem/autodock_protein.hpp>
#include <mudock/grid.hpp>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace mudock {

  //===------------------------------------------------------------------------------------------------------
  // Global parameters for deriving the pre-computation grid
  //===------------------------------------------------------------------------------------------------------

  static constexpr auto boundaries_cutoff   = fp_type{8};
  static constexpr auto resolution          = fp_type{0.5};
  static constexpr auto num_radius_tick     = std::size_t{2048};
  static constexpr auto num_radius_angstrom = fp_type{20.48};

  // this is the list of atom types that we want to consider while building the atom maps
  static constexpr auto target_atom_types = std::array<autodock_ff, 12>{{autodock_ff::A,
                                                                         autodock_ff::C,
                                                                         autodock_ff::H,
                                                                         autodock_ff::HD,
                                                                         autodock_ff::N,
                                                                         autodock_ff::NA,
                                                                         autodock_ff::OA,
                                                                         autodock_ff::SA,
                                                                         autodock_ff::Cl,
                                                                         autodock_ff::F,
                                                                         autodock_ff::S,
                                                                         autodock_ff::Br}};

  //===------------------------------------------------------------------------------------------------------
  // Utility functions to derive the energy table for the vdw component
  //===------------------------------------------------------------------------------------------------------

  // this struct define the vdw energy shape of the interaction between two different atom types
  struct vdw_shape {
    autodock_ff ligand_type;
    autodock_ff protein_type;
    fp_type nbp_r   = 0;     // radius of energy-well minimum
    fp_type nbp_eps = 0;     // depth of energy-well minimum
    fp_type xA      = 12;    // exponent for lennard jones, generally 12
    fp_type xB      = 6;     // exponent for lennard jones, 6 for non-hbonders 10 for h-bonders
    bool hbonder    = false; // if the interaction between the two atoms can form an hbond
  };

  // help function that compute the vdw energy shapes for the target types
  // It returns a multi vector where each cell is the energy shape of the interaction between the two types
  // Dimension 0 -> protein types, Dimension 1 -> ligand types
  static auto compute_vdw_interaction_shapes() {
    md_vector<vdw_shape, 2> vdw_shapes(target_atom_types.size(), target_atom_types.size());
    for (std::size_t lig_type_index = 0; lig_type_index < vdw_shapes.size<1>(); ++lig_type_index) {
      const auto lig_type  = target_atom_types[lig_type_index];
      const auto& lig_desc = get_description(lig_type);
      for (std::size_t prot_type_index = 0; prot_type_index < vdw_shapes.size<0>(); ++prot_type_index) {
        auto& shape           = vdw_shapes.get(prot_type_index, lig_type_index);
        shape.ligand_type     = target_atom_types[lig_type_index];
        shape.protein_type    = target_atom_types[prot_type_index];
        const auto& prot_desc = get_description(shape.protein_type);
        shape.nbp_r           = (lig_desc.Rii + prot_desc.Rii) / fp_type{2};
        shape.nbp_eps         = std::sqrt(lig_desc.epsii * prot_desc.epsii);
        shape.hbonder         = lig_desc.hbond > fp_type{0} ? true : false;
        if (lig_desc.hbond > 2 && (prot_desc.hbond == 1 || prot_desc.hbond == 2)) {
          shape.xB      = 10;
          shape.nbp_r   = lig_desc.Rij_hb;
          shape.nbp_eps = lig_desc.epsij_hb;
        } else if ((lig_desc.hbond == 1 || lig_desc.hbond == 2) && prot_desc.hbond > 2) {
          shape.xB      = 10;
          shape.nbp_r   = prot_desc.Rij_hb;
          shape.nbp_eps = prot_desc.epsij_hb;
        }
      }
    }
    return vdw_shapes;
  }

  // this function will compute the actual values for the vdw energy shapes, as function of the distance
  // It returns a multi vector of float with the energy values
  // Dimension 0 -> Energy values as changing the distance
  // Dimension 1 -> protein types, Dimension 2 -> ligand types
  // NOTE: we apply a smoothing function on the theoretical value by looking at previous value
  static auto compute_vdw_interaction_energies() {
    static const md_vector<vdw_shape, 2> vdw_shapes = compute_vdw_interaction_shapes();
    md_vector<fp_type, 3> energy_table(vdw_shapes.size<0>(), vdw_shapes.size<1>(), num_radius_tick);
    for (std::size_t lig_type_index = 0; lig_type_index < vdw_shapes.size<2>(); ++lig_type_index) {
      for (std::size_t prot_type_index = 0; prot_type_index < vdw_shapes.size<1>(); ++prot_type_index) {
        // compute the  energy shape
        const auto& shape     = vdw_shapes.get(lig_type_index, prot_type_index);
        const auto temp_const = shape.nbp_eps / (shape.xA - shape.xB);
        const auto cA         = temp_const * std::pow(shape.nbp_r, shape.xA) * shape.xB;
        const auto cB         = temp_const * std::pow(shape.nbp_r, shape.xB) * shape.xA;
        auto& voxel           = energy_table.get(0, prot_type_index, lig_type_index);
        voxel                 = fp_type{0};
        for (std::size_t radius_index = 1; radius_index < num_radius_tick; ++radius_index) {
          const auto radius = static_cast<fp_type>(radius_index) *
                              (num_radius_angstrom / static_cast<fp_type>(num_radius_tick));
          const auto rA              = std::pow(radius, shape.xA);
          const auto rB              = std::pow(radius, shape.xB);
          const auto proposed_energy = cA / rA - cB / rB;
          const auto clapped_energy  = std::min(fp_type{100000}, proposed_energy);
          energy_table.get(radius_index, prot_type_index, lig_type_index) = clapped_energy;
        }

        // apply a "smoothing" process that actually discretize the shape
        // NOTE: the diameter of the smoothing is 0.5A
        static constexpr auto smooth_diameter = static_cast<std::size_t>(
            (static_cast<fp_type>(num_radius_tick) / num_radius_angstrom) * fp_type{0.5});
        static constexpr auto smooth_radius = smooth_diameter / std::size_t{2};
        if constexpr (auto* energy_row = &voxel; smooth_diameter < num_radius_tick) {
          for (std::size_t radius_index = 0; radius_index < smooth_radius; radius_index++) {
            auto min_value       = energy_row[0];
            const auto end_index = radius_index + smooth_radius;
            for (std::size_t i = 1; i < end_index; ++i) { min_value = std::min(min_value, energy_row[i]); }
            energy_row[radius_index] = min_value;
          }
          static constexpr auto middle_end = num_radius_tick - smooth_radius;
          for (std::size_t radius_index = smooth_radius; radius_index < middle_end; ++radius_index) {
            const auto begin_index   = radius_index - smooth_radius;
            const auto end_condition = radius_index + smooth_radius;
            auto min_value           = energy_row[begin_index];
            for (std::size_t i = begin_index + std::size_t{1}; i < end_condition; ++i) {
              min_value = std::min(min_value, energy_row[i]);
            }
            energy_row[radius_index] = min_value;
          }
          for (std::size_t radius_index = middle_end; radius_index < num_radius_tick; ++radius_index) {
            const auto begin_index = radius_index - smooth_radius;
            auto min_value         = energy_row[begin_index];
            for (std::size_t i = begin_index + std::size_t{1}; i < num_radius_tick; ++i) {
              min_value = std::min(min_value, energy_row[i]);
            }
            energy_row[radius_index] = min_value;
          }
        } else { // we have really low vlaues
          const auto min_value = *std::min_element(energy_row, energy_row + num_radius_tick);
          for (std::size_t i = 0; i < num_radius_tick; ++i) { energy_row[i] = min_value; }
        }
      }
    }
    return energy_table;
  }

  //===------------------------------------------------------------------------------------------------------
  // Utility function to derive the desolvation energy by varying the distance
  //===------------------------------------------------------------------------------------------------------

  // this function will compute the desolvation energy among the two atoms types
  // NOTE: this function is independent from the atom types
  static auto compute_desolvation_energy() {
    std::array<fp_type, num_radius_tick> energy_table;
    for (std::size_t radius_index = 0; radius_index < num_radius_tick; ++radius_index) {
      const auto radius =
          static_cast<fp_type>(radius_index) * (num_radius_angstrom / static_cast<fp_type>(num_radius_tick));
      energy_table[radius_index] = autodock_parameters::coeff_desolv *
                                   std::exp(-(radius * radius) / (fp_type{2} * fp_type{3.6} * fp_type{3.6}));
    }
    return energy_table;
  }

  //===------------------------------------------------------------------------------------------------------
  // Utility functions to compute the protein HB gemetries
  //===------------------------------------------------------------------------------------------------------

  // this struct represents the geometries to compute the Hbonds
  struct hbond_geometries {
    std::vector<point3D> vector1;
    std::vector<point3D> vector2;
    std::vector<fp_type> exp;
    std::vector<int> disorder;

    inline hbond_geometries(const std::size_t n)
        : vector1(n, point3D{0, 0, 0}), vector2(n, point3D{0, 0, 0}), exp(n, fp_type{0}), disorder(n, 0) {}
  };

  // this function compute the hbon geometries
  hbond_geometries compute_hbon_geometries(const std::span<const fp_type> x,
                                           const std::span<const fp_type> y,
                                           const std::span<const fp_type> z,
                                           const std::span<const std::size_t> hbond,
                                           const std::span<const element> elements,
                                           const molecule_graph_type& graph) {
    const auto num_atoms = x.size();
    assert(y.size() == num_atoms);
    assert(z.size() == num_atoms);
    assert(hbond.size() == num_atoms);
    assert(elements.size() == num_atoms);
    auto result                       = hbond_geometries{num_atoms};
    const auto [atom_begin, atom_end] = boost::vertices(graph);
    for (auto it = atom_begin; it != atom_end; ++it) {
      const auto atom_index  = graph[*it].atom_index;
      const auto hbond_value = hbond[atom_index];
      const auto atom_point  = point3D{x[atom_index], y[atom_index], z[atom_index]};
      if (hbond_value == std::size_t{2}) { // ----------------------------------  D1 hydrogen bond donor
        if (elements[atom_index] != element::H) [[unlikely]]
          throw std::runtime_error("Unexpected atom type");
        const auto num_neighs = boost::out_degree(*it, graph);
        if (num_neighs > std::size_t{1}) [[unlikely]]
          throw std::runtime_error("Unsupported number of neighbors");
        const auto [neigh_begin, neigh_end] = boost::out_edges(*it, graph);
        for (auto neigh = neigh_begin; neigh != neigh_end; ++neigh) {
          const auto neigh_index = graph[neigh->m_target].atom_index;
          const auto neigh_point = point3D{x[neigh_index], y[neigh_index], z[neigh_index]};
          const auto diff        = difference(atom_point, neigh_point);
          const auto d2          = sum_components(square(diff));
          if (d2 < fp_type{1.9}) {
            const auto neigh_element = elements[neigh_index];
            if (neigh_element == element::O || neigh_element == element::S) {
              result.exp[atom_index]      = fp_type{4};
              result.disorder[atom_index] = 1;
            } else {
              result.exp[atom_index]      = fp_type{2};
              result.disorder[atom_index] = 1;
            }
            result.vector1[atom_index] = normalize(diff);
          }
        }
      } else if (hbond_value == std::size_t{5}) { // ----------------------------------  A2 oxygen
        if (elements[atom_index] != element::O) [[unlikely]]
          throw std::runtime_error("Unexpected atom type");
        const auto [neigh_begin, neigh_end] = boost::out_edges(*it, graph);
        auto bond_counter                   = std::size_t{0};
        auto neigh1_vertex                  = std::size_t{0};
        auto neigh1_index                   = std::size_t{0};
        auto neigh2_index                   = std::size_t{0};
        auto neigh1_point                   = point3D{};
        auto neigh2_point                   = point3D{};
        for (auto neigh = neigh_begin; neigh != neigh_end; ++neigh) {
          const auto neigh_vertex  = neigh->m_target;
          const auto neigh_index   = graph[neigh_vertex].atom_index;
          const auto neigh_point   = point3D{x[neigh_index], y[neigh_index], z[neigh_index]};
          const auto diff          = difference(atom_point, neigh_point);
          const auto d2            = sum_components(square(diff));
          const auto neigh_element = elements[neigh_index];
          if ((d2 < fp_type{3.61} && neigh_element != element::H) ||
              (d2 < fp_type{1.69} && neigh_element == element::H)) {
            switch (bond_counter) {
              case std::size_t{0}:
                bond_counter  = 1;
                neigh1_vertex = neigh_vertex;
                neigh1_index  = neigh_index;
                neigh1_point  = neigh_point;
                break;

              case std::size_t{1}:
                bond_counter = 2;
                neigh2_index = neigh_index;
                neigh2_point = neigh_point;
                break;

              default: throw std::runtime_error("Unsupported number of neighbors");
            }
          }
        }

        if (bond_counter == std::size_t{0}) [[unlikely]]
          throw std::runtime_error("Oxygen with no bonded atoms");
        else if (bond_counter == std::size_t{1}) { // in this case we have lone pairs
          // so we need to explore the neighbor's neighbor for Carbonyl Oxygen O=C-X
          if (elements[neigh1_index] != element::C) [[unlikely]]
            throw std::runtime_error("The original autogrid was not expecting a non C atom with a O");
          result.vector1[atom_index] = normalize(difference(atom_point, std::as_const(neigh1_point)));
          auto found                 = false;
          const auto [c_neigh_begin, c_neigh_end] = boost::out_edges(neigh1_vertex, graph);
          for (auto c_neigh = c_neigh_begin; c_neigh != c_neigh_end; ++c_neigh) {
            const auto c_neigh_index   = graph[c_neigh->m_target].atom_index;
            const auto c_neigh_point   = point3D{x[c_neigh_index], y[c_neigh_index], z[c_neigh_index]};
            const auto c_diff          = difference(std::as_const(neigh1_point), c_neigh_point);
            const auto c_d2            = sum_components(square(c_diff));
            const auto c_neigh_element = elements[c_neigh_index];
            if ((c_d2 < fp_type{2.89} && c_neigh_element != element::H) ||
                (c_d2 < fp_type{2.69} && c_neigh_element == element::H)) {
              found = true;
              // C=O cross C-X gives the lone pair plane normal
              result.vector2[atom_index] =
                  normalize(cross_product(std::as_const(result.vector1[atom_index]), c_neigh_point));
            }
          }
          if (!found) [[unlikely]]
            throw std::runtime_error("Unknown carbonyl oxygen");
        } else if (bond_counter == std::size_t{2}) { // in this case we assume that we have either C and H or
          // just two H, maybe it is due chemical rules
          const auto neigh1_element = elements[neigh1_index];
          const auto neigh2_element = elements[neigh2_index];
          if (neigh1_element == element::H && neigh2_element == element::C) {
            result.vector1[atom_index] = normalize(difference(atom_point, neigh1_point));
          } else if (neigh1_element == element::C && neigh2_element == element::H) {
            result.vector1[atom_index] = normalize(difference(atom_point, neigh2_point));
          } else if (neigh1_element == element::H && neigh2_element == element::H) {
            result.vector2[atom_index] = normalize(difference(neigh2_point, neigh1_point));
            result.vector1[atom_index] =
                normalize(add(scale(std::as_const(result.vector2[atom_index]),
                                    sum_components(difference(atom_point, std::as_const(neigh1_point)))),
                              atom_point));
          } else [[unlikely]]
            throw std::runtime_error("The original autogrid was not expecting a non C atom");
        }
      } else if (hbond_value == std::size_t{4}) { // ----------------------------------  A1 nitrogen
        if (elements[atom_index] != element::N) [[unlikely]]
          throw std::runtime_error("Unexpected atom type");
        const auto [neigh_begin, neigh_end] = boost::out_edges(*it, graph);
        auto bond_counter                   = std::size_t{0};
        auto neigh1_index                   = std::size_t{0};
        auto neigh1_point                   = point3D{};
        auto neigh2_point                   = point3D{};
        auto neigh3_point                   = point3D{};
        for (auto neigh = neigh_begin; neigh != neigh_end; ++neigh) {
          const auto neigh_vertex  = neigh->m_target;
          const auto neigh_index   = graph[neigh_vertex].atom_index;
          const auto neigh_point   = point3D{x[neigh_index], y[neigh_index], z[neigh_index]};
          const auto diff          = difference(atom_point, neigh_point);
          const auto d2            = sum_components(square(diff));
          const auto neigh_element = elements[neigh_index];
          if ((d2 < fp_type{2.89} && neigh_element != element::H) ||
              (d2 < fp_type{1.69} && neigh_element == element::H)) {
            switch (bond_counter) {
              case std::size_t{0}:
                bond_counter = 1;
                neigh1_index = neigh_index;
                neigh1_point = neigh_point;
                break;

              case std::size_t{1}:
                bond_counter = 2;
                neigh2_point = neigh_point;
                break;

              case std::size_t{2}:
                bond_counter = 3;
                neigh3_point = neigh_point;
                break;

              default: throw std::runtime_error("Unsupported number of neighbors");
            }
          }
        }

        if (bond_counter == std::size_t{0}) [[unlikely]]
          throw std::runtime_error("Nitrogen with no bonded atoms");
        else if (bond_counter == std::size_t{1}) { // Azide Nitrogen N=C bond vector
          if (elements[neigh1_index] != element::C) [[unlikely]]
            throw std::runtime_error("The original autogrid was not expecting a non C atom with an N");
          result.vector1[atom_index] = normalize(difference(atom_point, neigh1_point));
        } else if (bond_counter == std::size_t{2}) { // two bonds: X1-N=X2
          result.vector1[atom_index] = normalize(
              difference(atom_point, scale(add(neigh1_point, neigh2_point), fp_type{1} / fp_type{2})));
        } else if (bond_counter == std::size_t{3}) { // three bonds
          result.vector1[atom_index] = normalize(difference(
              atom_point,
              scale(
                  add(std::as_const(neigh1_point), std::as_const(neigh2_point), std::as_const(neigh3_point)),
                  fp_type{1} / fp_type{3})));
        }
      }
    }
    return result;
  }

  //===------------------------------------------------------------------------------------------------------
  // Utility function to compute the dielectric Ewald Summation (Mehler and Solmajer, Prot Eng 4, 903-910)
  //===------------------------------------------------------------------------------------------------------

  static auto compute_dielectric_ewds() {
    static constexpr auto lambda   = fp_type{0.003627};
    static constexpr auto epsilon0 = fp_type{78.4};
    static constexpr auto A        = fp_type{-8.5525};
    static constexpr auto B        = epsilon0 - A;
    static constexpr auto rk       = fp_type{7.7839};
    static constexpr auto lambda_B = -lambda * B;

    std::array<fp_type, num_radius_tick> result;
    result[0] = fp_type{1};
    for (std::size_t radius_index = 0; radius_index < num_radius_tick; ++radius_index) {
      const auto radius =
          static_cast<fp_type>(radius_index) * (num_radius_angstrom / static_cast<fp_type>(num_radius_tick));
      result[radius_index] = fp_type{332} / (A + B / (fp_type{1} + rk * std::exp(lambda_B * radius)));
    }
    return result;
  }

  //===------------------------------------------------------------------------------------------------------
  // Implementation of the actual function that computes the grid
  //===------------------------------------------------------------------------------------------------------

  autodock_protein make_autodock_protein(const dynamic_molecule& protein, const molecule_graph_type& graph) {
    // compute the energy tables that are actually independent from the protein
    static const auto vdw_energies         = compute_vdw_interaction_energies();
    static const auto desolvation_energies = compute_desolvation_energy();
    static const auto dielectric_ewds      = compute_dielectric_ewds();

    // get the protein's atoms coordinate
    const auto x = protein.get_x();
    const auto y = protein.get_y();
    const auto z = protein.get_z();
    assert(x.size() == y.size());
    assert(x.size() == z.size());
    const auto num_atoms = x.size();

    // bail out if there are no protein atoms
    if (num_atoms == std::size_t{0}) {
      throw std::runtime_error("Protein without atoms");
    }

    // find the protein boundaries and its center
    auto min    = point3D{x[0], y[0], z[0]};
    auto max    = point3D{x[0], y[0], z[0]};
    auto center = point3D{x[0], y[0], z[0]};
    for (std::size_t i = 1; i < num_atoms; ++i) {
      min.x = std::min(min.x, x[i]);
      min.y = std::min(min.y, y[i]);
      min.z = std::min(min.z, z[i]);
      max.x = std::max(max.x, x[i]);
      max.y = std::max(max.y, y[i]);
      max.z = std::max(max.z, z[i]);
      center.x += x[i];
      center.y += y[i];
      center.z += z[i];
    }
    center.x /= static_cast<fp_type>(num_atoms);
    center.y /= static_cast<fp_type>(num_atoms);
    center.z /= static_cast<fp_type>(num_atoms);

    // cap the boundaries in a box of 8A of radius
    min.x = std::max(min.x, center.x - fp_type{boundaries_cutoff});
    min.y = std::max(min.y, center.y - fp_type{boundaries_cutoff});
    min.z = std::max(min.z, center.z - fp_type{boundaries_cutoff});
    max.x = std::max(max.x, center.x + fp_type{boundaries_cutoff});
    max.y = std::max(max.y, center.y + fp_type{boundaries_cutoff});
    max.z = std::max(max.z, center.z + fp_type{boundaries_cutoff});

    // find out the geometries of HBonds from the protein
    const auto [vector1, vector2, exp, disorder] =
        compute_hbon_geometries(x, y, z, protein.get_num_hbond(), protein.get_elements(), graph);

    // declare the maps that will describe the protein
    autodock_protein result;

    // allocate memory for the different grids
    result.electrostatic = space_grid(min, max, resolution);
    result.desolvation   = space_grid(min, max, resolution);
    for (const auto element: target_atom_types) {
      result.atom_map.emplace(element, space_grid(min, max, resolution));
    }

    // get the remaining protein information
    const auto charge = protein.get_charge();
    const auto volume = protein.get_vol();
    const auto size_x = result.electrostatic.size<0>();
    const auto size_y = result.electrostatic.size<1>();
    const auto size_z = result.electrostatic.size<2>();
    for (std::size_t index_z = 0; index_z < size_z; ++index_z) {
      for (std::size_t index_y = 0; index_y < size_y; ++index_y) {
        for (std::size_t index_x = 0; index_x < size_x; ++index_x) {
          const auto voxel_point    = result.electrostatic.to_coord(index_x, index_y, index_z);
          auto electrostatic_energy = fp_type{0};
          auto desolvation_energy   = fp_type{0};
          for (std::size_t i = 0; i < num_atoms; ++i) {
            const auto atom_point = point3D{x[i], y[i], z[i]};
            const auto d          = distance(atom_point, voxel_point);
            if (d < boundaries_cutoff) {
              const auto radius_index =
                  std::min(num_radius_tick - std::size_t{1},
                           static_cast<std::size_t>(
                               d * (static_cast<fp_type>(num_radius_tick) / num_radius_angstrom)));
              electrostatic_energy +=
                  charge[i] * (fp_type{1} / std::max(d, fp_type{0.5})) * autodock_parameters::coeff_estat;
              desolvation_energy += fp_type{0.01097} * volume[i] * desolvation_energies[radius_index];
            }
          }

          // commit the values in the actual grid maps
          result.electrostatic.get(index_x, index_y, index_z) = electrostatic_energy;
        }
      }
    }

    return result;
  }

} // namespace mudock
