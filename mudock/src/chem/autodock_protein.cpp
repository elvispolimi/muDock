#include <algorithm>
#include <array>
#include <cassert>
#include <csignal>
#include <iomanip>
#include <limits>
#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/chem/autodock_protein.hpp>
#include <mudock/chem/autodock_types.hpp>
#include <mudock/chem/grid_const.hpp>
#include <mudock/chem/mehler_solmajer.hpp>
#include <mudock/grid.hpp>
#include <mudock/grid/space_grid.hpp>
#include <mudock/type_alias.hpp>
#include <mudock/utils.hpp>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mudock {

  autodock_protein::autodock_protein(const point3D min, const point3D max, const fp_type resolution)
      : index(static_cast<int>(std::ceil(std::abs(min.x - max.x) / resolution) + 1),
              static_cast<int>(std::ceil(std::abs(min.y - max.y) / resolution) + 1),
              static_cast<int>(std::ceil(std::abs(min.z - max.z) / resolution) + 1)),
        data(index.size_x(), index.size_y(), index.size_z(), num_autodock_grids()),
        _inv_resolution(1 / resolution),
        _min(min),
        _max(max),
        _center{(max.x - min.x) / fp_type{2} + min.x,
                (max.y - min.y) / fp_type{2} + min.y,
                (max.z - min.z) / fp_type{2} + min.z} {};

  //===------------------------------------------------------------------------------------------------------
  // Global parameters for deriving the pre-computation grid
  //===------------------------------------------------------------------------------------------------------

  static constexpr auto energy_cutoff   = fp_type{8};
  static constexpr auto resolution      = fp_type{0.5};
  static constexpr auto half_resolution = resolution / fp_type{2};

  //===------------------------------------------------------------------------------------------------------
  // Utility functions to compute the protein HB gemetries
  //===------------------------------------------------------------------------------------------------------

  // this struct represents the geometries to compute the Hbonds
  struct hbond_geometries {
    std::vector<point3D> vector1;
    std::vector<point3D> vector2;
    std::vector<int> exp;
    std::vector<int> disorder;

    inline hbond_geometries(const std::size_t n)
        : vector1(n, point3D{0, 0, 0}), vector2(n, point3D{0, 0, 0}), exp(n, fp_type{0}), disorder(n, 0) {}
  };

  // this function compute the hbon geometries
  hbond_geometries compute_hbon_geometries(const std::span<const fp_type> x,
                                           const std::span<const fp_type> y,
                                           const std::span<const fp_type> z,
                                           const std::span<const int> hbond,
                                           const std::span<const element> elements,
                                           const molecule_graph_type& graph) {
    const auto num_atoms = x.size();
    assert(y.size() == num_atoms);
    assert(z.size() == num_atoms);
    assert(hbond.size() == num_atoms);
    assert(elements.size() == num_atoms);
    auto adt_protein                  = hbond_geometries{num_atoms};
    const auto [atom_begin, atom_end] = boost::vertices(graph);
    for (auto it = atom_begin; it != atom_end; ++it) {
      const auto atom_index  = graph[*it].atom_index;
      const auto hbond_value = hbond[atom_index];
      const auto atom_point  = point3D{x[atom_index], y[atom_index], z[atom_index]};
      const auto temp        = atom_index - 20;
      // const auto from_neighbor = temp <= index ? temp : int{0};
      const auto from_neighbor = std::max(temp, 0);
      const auto to_neighbor   = std::min(atom_index + 20, static_cast<int>(num_atoms));
      if (hbond_value == std::size_t{2}) { // ----------------------------------  D1 hydrogen bond donor
        if (elements[atom_index] != element::H) [[unlikely]]
          throw std::runtime_error("Unexpected atom type");
        // const auto num_neighs = boost::out_degree(*it, graph);
        // if (num_neighs > std::size_t{1}) [[unlikely]]
        //   throw std::runtime_error("Unsupported number of neighbors");
        // const auto [neigh_begin, neigh_end] = boost::out_edges(*it, graph);
        for (int neigh_index = from_neighbor; neigh_index < to_neighbor; ++neigh_index)
          if (neigh_index != atom_index) {
            // for (auto neigh = neigh_begin; neigh != neigh_end; ++neigh) {
            // const auto neigh_index = graph[neigh->m_target].atom_index;
            const auto neigh_point = point3D{x[neigh_index], y[neigh_index], z[neigh_index]};
            const auto diff        = difference(atom_point, neigh_point);
            const auto d2          = sum_components(square(diff));
            if (d2 < fp_type{1.9}) {
              const auto neigh_element = elements[neigh_index];
              if (neigh_element == element::O || neigh_element == element::S) {
                adt_protein.exp[atom_index]      = 4;
                adt_protein.disorder[atom_index] = 1;
              } else {
                adt_protein.exp[atom_index]      = 2;
                adt_protein.disorder[atom_index] = 1;
              }
              adt_protein.vector1[atom_index] = normalize(diff);
              // }
            }
          }
      } else if (hbond_value == std::size_t{5}) { // ----------------------------------  A2 oxygen
        // if (elements[atom_index] != element::O) [[unlikely]]
        //   throw std::runtime_error("Unexpected atom type");
        // const auto [neigh_begin, neigh_end] = boost::out_edges(*it, graph);
        auto bond_counter = 0;
        // auto neigh1_vertex = std::size_t{0};
        auto neigh1_index = 0;
        auto neigh1_point = point3D{};
        auto neigh2_point = point3D{};
        // for (auto neigh = neigh_begin; neigh != neigh_end; ++neigh) {
        for (int neigh_index = from_neighbor; neigh_index < to_neighbor; ++neigh_index)
          if (neigh_index != atom_index) {
            // const auto neigh_vertex  = neigh->m_target;
            // const auto neigh_index   = graph[neigh_vertex].atom_index;
            const auto neigh_point   = point3D{x[neigh_index], y[neigh_index], z[neigh_index]};
            const auto diff          = difference(atom_point, neigh_point);
            const auto d2            = sum_components(square(diff));
            const auto neigh_element = elements[neigh_index];
            if ((d2 < fp_type{3.61} && neigh_element != element::H) ||
                (d2 < fp_type{1.69} && neigh_element == element::H)) {
              switch (bond_counter) {
                case std::size_t{0}:
                  bond_counter = 1;
                  // neigh1_vertex = neigh_vertex;
                  neigh1_index = neigh_index;
                  neigh1_point = neigh_point;
                  break;

                case std::size_t{1}:
                  bond_counter = 2;
                  neigh2_point = neigh_point;
                  break;

                default: throw std::runtime_error("Unsupported number of neighbors");
              }
            }
          }

        // if (bond_counter == std::size_t{0}) [[unlikely]]
        //   throw std::runtime_error("Oxygen with no bonded atoms");
        if (bond_counter == std::size_t{1}) { // in this case we have lone pairs
          // so we need to explore the neighbor's neighbor for Carbonyl Oxygen O=C-X
          if (elements[neigh1_index] != element::C) [[unlikely]]
            throw std::runtime_error("The original autogrid was not expecting a non C atom with a O");
          adt_protein.vector1[atom_index] = normalize(difference(atom_point, std::as_const(neigh1_point)));
          // auto found                      = false;
          // const auto [c_neigh_begin, c_neigh_end] = boost::out_edges(neigh1_vertex, graph);
          // for (auto c_neigh = c_neigh_begin; c_neigh != c_neigh_end; ++c_neigh) {
          for (int other_index = from_neighbor; other_index < to_neighbor; ++other_index) {
            if ((other_index != neigh1_index) && (other_index != atom_index)) {
              // const auto c_neigh_index   = graph[c_neigh->m_target].atom_index;
              const auto c_neigh_point   = point3D{x[other_index], y[other_index], z[other_index]};
              const auto c_diff          = difference(std::as_const(neigh1_point), c_neigh_point);
              const auto c_d2            = sum_components(square(c_diff));
              const auto c_norm          = normalize(c_diff);
              const auto c_neigh_element = elements[other_index];
              if ((c_d2 < fp_type{2.89} && c_neigh_element != element::H) ||
                  (c_d2 < fp_type{1.69} && c_neigh_element == element::H)) {
                // found = true;
                // C=O cross C-X gives the lone pair plane normal
                adt_protein.vector2[atom_index] =
                    normalize(cross_product(std::as_const(adt_protein.vector1[atom_index]), c_norm));
              }
            }
            // if (!found) [[unlikely]]
            //   throw std::runtime_error("Unknown carbonyl oxygen");
          }
        } else if (bond_counter == std::size_t{2}) { // in this case we assume that we have either C and H or
          // just two H, maybe it is due chemical rules
          // const auto neigh1_element = elements[neigh1_index];
          // const auto neigh2_element = elements[neigh2_index];
          // if (neigh1_element == element::H && neigh2_element == element::C) {
          //   adt_protein.vector1[atom_index] = normalize(difference(atom_point, neigh1_point));
          // } else if (neigh1_element == element::C && neigh2_element == element::H) {
          //   adt_protein.vector1[atom_index] = normalize(difference(atom_point, neigh2_point));
          // } else if (neigh1_element == element::H && neigh2_element == element::H) {
          //   adt_protein.vector2[atom_index] = normalize(difference(neigh2_point, neigh1_point));
          //   const auto p                    = scale(std::as_const(adt_protein.vector2[atom_index]),
          //                        sum_components(difference(atom_point, std::as_const(neigh1_point))));
          //   adt_protein.vector1[atom_index] = normalize(add(p, atom_point));
          // } else [[unlikely]]
          //   throw std::runtime_error("The original autogrid was not expecting a non C atom");
          // rvector2[index] = normalize(point3D{receptor.x(index_2), receptor.y(index_2), receptor.z(index_2)},
          adt_protein.vector2[atom_index] = normalize(neigh2_point, neigh1_point);

          const point3D diff = difference(atom_point, neigh1_point);
          const auto rdot    = (diff * adt_protein.vector2[atom_index]).sum_components();

          adt_protein.vector1[atom_index] = normalize(
              difference(atom_point, add(scale(adt_protein.vector2[atom_index], rdot), neigh1_point)));
        }
      } else if (hbond_value == std::size_t{4}) { // ----------------------------------  A1 nitrogen
        if (elements[atom_index] != element::N) [[unlikely]]
          throw std::runtime_error("Unexpected atom type");
        // const auto [neigh_begin, neigh_end] = boost::out_edges(*it, graph);
        auto bond_counter = std::size_t{0};
        auto neigh1_index = std::size_t{0};
        auto neigh1_point = point3D{};
        auto neigh2_point = point3D{};
        auto neigh3_point = point3D{};
        // for (auto neigh = neigh_begin; neigh != neigh_end; ++neigh) {
        for (int neigh_index = from_neighbor; neigh_index < to_neighbor; ++neigh_index)
          if (neigh_index != atom_index) {
            // const auto neigh_vertex  = neigh->m_target;
            // const auto neigh_index   = graph[neigh_vertex].atom_index;
            const auto neigh_point   = point3D{x[neigh_index], y[neigh_index], z[neigh_index]};
            const auto diff          = difference(atom_point, neigh_point);
            const auto d2            = sum_components(square(diff));
            const auto neigh_element = elements[neigh_index];
            if ((d2 < fp_type{3.61} && neigh_element != element::H) ||
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
          adt_protein.vector1[atom_index] = normalize(difference(atom_point, neigh1_point));
        } else if (bond_counter == std::size_t{2}) { // two bonds: X1-N=X2
          adt_protein.vector1[atom_index] = normalize(
              difference(atom_point, scale(add(neigh1_point, neigh2_point), fp_type{1} / fp_type{2})));
        } else if (bond_counter == std::size_t{3}) { // three bonds
          const auto p1 = add(std::as_const(neigh1_point), std::as_const(neigh2_point));
          const auto p2 = scale(add(p1, std::as_const(neigh3_point)), fp_type{1} / fp_type{3});
          adt_protein.vector1[atom_index] = normalize(difference(atom_point, p2));
        }
      }
    }
    return adt_protein;
  }

  struct voxel_scratchpad {
    fp_type hbondmin{999999};
    fp_type hbondmax{-999999};
    fp_type hbondflag{false};
    fp_type energy{0};
  };

  //===------------------------------------------------------------------------------------------------------
  // Implementation of the actual function that computes the grid
  //===------------------------------------------------------------------------------------------------------

  autodock_protein make_autodock_protein(const dynamic_molecule& protein) {
    static const auto vdw_shapes = compute_vdw_interaction_shapes();
    auto graph                   = make_graph(protein.get_bonds(), protein.num_atoms());

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
    auto min = point3D{std::ranges::min(x), std::ranges::min(y), std::ranges::min(z)};
    std::transform(min.begin(), min.end(), min.begin(), [](const auto p) {
      return std::floor((p - (grid_spacing + cutoff_distance)) * inv_spacing) / inv_spacing;
    });
    auto max = point3D{std::ranges::max(x), std::ranges::max(y), std::ranges::max(z)};
    std::transform(max.begin(), max.end(), max.begin(), [](const auto p) {
      return std::ceil((p + (grid_spacing + cutoff_distance)) * inv_spacing) / inv_spacing;
    });

    // find out the geometries of HBonds from the protein
    const auto [vector1, vector2, exp, disorder] =
        compute_hbon_geometries(x, y, z, protein.get_num_hbond(), protein.get_elements(), graph);

    // declare the maps that will describe the protein
    autodock_protein adt_protein{min, max, resolution};

    // get the remaining protein information
    const auto charge         = protein.get_charge();
    const auto volume         = protein.get_vol();
    const auto num_hbonds     = protein.get_num_hbond();
    const auto autodock_types = protein.get_autodock_type();
    const auto size_x         = adt_protein.get_eletrostatic().size<0>();
    const auto size_y         = adt_protein.get_eletrostatic().size<1>();
    const auto size_z         = adt_protein.get_eletrostatic().size<2>();
    for (std::size_t index_z = 0; index_z < size_z; ++index_z) {
      for (std::size_t index_y = 0; index_y < size_y; ++index_y) {
        for (std::size_t index_x = 0; index_x < size_x; ++index_x) {
          std::array<voxel_scratchpad, num_autodock_ff_grids()> voxel_scratchs;
          // get the voxel point
          const auto voxel_point = adt_protein.get_eletrostatic().to_coord(static_cast<fp_type>(index_x),
                                                                           static_cast<fp_type>(index_y),
                                                                           static_cast<fp_type>(index_z));

          // add electrostatic and desolvation energy + find the nearest Hbond
          auto nearest_H_index      = std::size_t{0};
          auto nearest_H_distance   = point3D{x[0], y[0], z[0]}.distance(voxel_point);
          auto nearest_H_valid      = num_hbonds[0] == 1 || num_hbonds[0] == 2;
          auto electrostatic_energy = fp_type{0};
          auto desolvation_energy   = fp_type{0};
          for (std::size_t i = 0; i < num_atoms; ++i) {
            // compute properties of the given atom
            const auto atom_point = point<fp_type, 3>{x[i], y[i], z[i]};
            const auto d          = atom_point.distance(voxel_point);
            // const auto inv_d      = fp_type{1} / d;
            const auto inv_dmax = fp_type{1} / std::max(fp_type{0.5}, d);
            const auto indx_r   = std::min<int>(std::floor(d * lookup_resolution), num_radius_tick_elect);

            // add the electrostatic constribution
            electrostatic_energy += charge[i] * inv_dmax * adt_protein.electrostatic_energies[indx_r] *
                                    autodock_parameters::coeff_estat;

            // find the nearest Hbond
            const auto is_valid = num_hbonds[i] == 1 || num_hbonds[i] == 2;
            if ((d < nearest_H_distance && is_valid) || (!nearest_H_valid && is_valid)) {
              nearest_H_distance = d;
              nearest_H_index    = i;
              nearest_H_valid    = is_valid;
            }

            // add the desolvation contribute
            if (d <= energy_cutoff) {
              const auto radius_index = std::min(num_radius_tick_desolv - std::size_t{1},
                                                 static_cast<std::size_t>(d * lookup_resolution));
              desolvation_energy +=
                  fp_type{0.01097} * volume[i] * adt_protein.desolvation_energies[radius_index];
            }
          }

          // commit the values in the actual grid maps
          adt_protein.get_eletrostatic().get(index_x, index_y, index_z) = electrostatic_energy;
          adt_protein.get_desolvation().get(index_x, index_y, index_z)  = desolvation_energy;

          // find out the Hbond parameters
          if (nearest_H_valid) {
            for (std::size_t i = 0; i < num_atoms; ++i) {
              const auto atom_point = point3D{x[i], y[i], z[i]};
              const auto diff       = (atom_point - voxel_point).normalize();
              const auto d          = atom_point.distance(voxel_point);
              if (d <= cutoff_distance) {
                auto racc = fp_type{1}, rdon = fp_type{1}, Hramp = fp_type{1}, cos_theta = fp_type{0};
                switch (num_hbonds[i]) {
                  case std::size_t{2}:
                    cos_theta = -diff.product(vector1[i]).sum_components();
                    if (cos_theta <= fp_type{0}) {
                      racc = fp_type{0};
                    } else {
                      switch (exp[i]) {
                        case 0: racc = cos_theta; break;
                        case 1: racc = cos_theta; break;
                        case 2: racc = cos_theta * cos_theta; break;
                        case 4:
                          racc = cos_theta * cos_theta;
                          racc = racc * racc;
                          break;
                        default: throw std::runtime_error("Unexpected exponent " + std::to_string(exp[i]));
                      }
                      if (i == nearest_H_index) {
                        Hramp = fp_type{1};
                      } else {
                        cos_theta = sum_components(product(vector1[nearest_H_index], vector1[i]));
                        cos_theta = std::min(fp_type{1}, std::max(cos_theta, fp_type{-1}));
                        Hramp     = fp_type{0.5} -
                                fp_type{0.5} * std::cos(std::acos(cos_theta) * fp_type{120} / fp_type{90});
                      }
                    }
                    break;

                  case std::size_t{4}:
                    cos_theta = -diff.product(vector1[i]).sum_components();
                    if (cos_theta <= fp_type{0}) {
                      rdon = fp_type{0};
                    } else {
                      rdon = cos_theta * cos_theta;
                    }
                    break;

                  case std::size_t{5}: {
                    cos_theta = -diff.product(vector1[i]).sum_components();
                    const auto t0 =
                        math::pi_halved -
                        std::acos(
                            std::clamp(diff.product(vector2[i]).sum_components(), fp_type{-1}, fp_type{1}));
                    const auto cross = diff.cross(vector2[i]).normalize();
                    fp_type ti       = cross.product(vector1[i]).sum_components();

                    /* rdon expressions from Goodford */
                    rdon = 0.;
                    if (cos_theta >= fp_type{0}) {
                      ti = std::clamp(ti, fp_type{-1}, fp_type{1});
                      ti = std::acos(ti) - math::pi_halved;
                      if (ti < 0) {
                        ti = -ti;
                      }
                      /* the 2.0*ti can be replaced by (ti + ti) in: rdon = (0.9 + 0.1*sin(2.0*ti))*cos(t0);*/
                      rdon = (fp_type{0.9} + fp_type{0.1} * std::sin(ti + ti)) * std::cos(t0);
                    } else if (cos_theta >= fp_type{-0.34202}) {
                      /* 0.34202 = cos (100 deg) */
                      // TODO @Davide ok here the fp_type?
                      rdon = fp_type{562.25} *
                             std::pow(fp_type{0.116978} - cos_theta * cos_theta, fp_type{3}) * std::cos(t0);
                    }
                    break;
                  }
                  default: break;
                }

                const auto& protein_desc = get_description(autodock_types[i]);

                const auto indx_n =
                    std::min<int>(std::floor(d * lookup_resolution), num_radius_tick_desolv - 1);
                const int indx_r =
                    std::min<int>(std::floor(d * lookup_resolution), num_radius_tick_elect - 1);

                for (int grid_index = 0; grid_index < num_autodock_ff_grids(); ++grid_index) {
                  auto& voxel_scratch                = voxel_scratchs[grid_index];
                  const autodock_grid_type grid_type = static_cast<autodock_grid_type>(grid_index);
                  const auto& ligand_desc            = get_description(autodock_ff_from_grid(grid_type));

                  const auto vdw_shape = vdw_shapes.get(static_cast<int>(autodock_types[i]), grid_index);
                  const auto vdw_hb_value =
                      adt_protein.vdw_energies.get(indx_n, static_cast<int>(autodock_types[i]), grid_index);

                  if (vdw_shape.hbonder) {
                    fp_type rsph = vdw_hb_value / fp_type{100};
                    rsph         = std::clamp(rsph, fp_type{0}, fp_type{1});
                    if ((ligand_desc.hbond == 3 || ligand_desc.hbond == 5)         /*AS or A2*/
                        && (protein_desc.hbond == 1 || protein_desc.hbond == 2)) { /*DS or D1*/
                      voxel_scratch.energy += vdw_hb_value * Hramp * (racc + (fp_type{1} - racc) * rsph);
                    } else if ((ligand_desc.hbond == 4)                                   /*A1*/
                               && (protein_desc.hbond == 1 || protein_desc.hbond == 2)) { /*DS,D1*/
                      voxel_scratch.hbondmin  = std::min(voxel_scratch.hbondmin,
                                                        vdw_hb_value * (racc + (fp_type{1} - racc) * rsph));
                      voxel_scratch.hbondmax  = std::max(voxel_scratch.hbondmax,
                                                        vdw_hb_value * (racc + (fp_type{1} - racc) * rsph));
                      voxel_scratch.hbondflag = true;
                    } else if ((ligand_desc.hbond == 1 || ligand_desc.hbond == 2) &&
                               (protein_desc.hbond > 2)) { /*DS,D1 vs AS,A1,A2*/
                      const fp_type temp_hbond_enrg = vdw_hb_value * (rdon + (fp_type{1} - rdon) * rsph);
                      voxel_scratch.hbondmin        = std::min(voxel_scratch.hbondmin, temp_hbond_enrg);
                      voxel_scratch.hbondmax        = std::max(voxel_scratch.hbondmax, temp_hbond_enrg);
                      voxel_scratch.hbondflag       = true;
                    } else { /*end of is_hbonder*/
                      voxel_scratch.energy += vdw_hb_value;
                    }
                  } else {
                    voxel_scratch.energy += vdw_hb_value;
                  } /* end hbonder tests */

                  voxel_scratch.energy +=
                      ligand_desc.solpar * protein_desc.vol * adt_protein.desolvation_energies[indx_r] +
                      (protein_desc.solpar + solpar_q * std::fabs(charge[i])) * ligand_desc.vol *
                          adt_protein.desolvation_energies[indx_r];
                }
              }
            }

            for (int grid_index = 0; grid_index < num_autodock_ff_grids(); ++grid_index) {
              auto& voxel_scratch = voxel_scratchs[grid_index];
              auto energy         = voxel_scratch.energy + voxel_scratch.hbondmin + voxel_scratch.hbondmax;

              if (std::fabs(energy) < precision)
                energy = 0;
              adt_protein.get_atom_map(static_cast<autodock_grid_type>(grid_index))
                  .get(index_x, index_y, index_z) = energy;
            }
          }
        }
      }
    }
    return adt_protein;
  }
} // namespace mudock
