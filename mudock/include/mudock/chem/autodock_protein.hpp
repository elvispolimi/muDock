#pragma once

#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/chem/autodock_layer.hpp>
#include <mudock/chem/autodock_parameters.hpp>
#include <mudock/chem/pdbqt_forcefield_override.hpp>
#include <mudock/chem/autodock_types.hpp>
#include <mudock/chem/grid_const.hpp>
#include <mudock/grid.hpp>
#include <mudock/grid/mdspan.hpp>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {
  static constexpr auto num_radius_tick_desolv = std::size_t{2048};
  static constexpr auto num_radius_tick_elect  = std::size_t{16384};
  static constexpr auto lookup_resolution      = fp_type{100}; /* Used in distance look-up table. */

  //===------------------------------------------------------------------------------------------------------
  // Utility functions to derive the energy table for the vdw component
  //===------------------------------------------------------------------------------------------------------

  // this struct define the vdw energy shape of the interaction between two different atom types
  struct vdw_shape {
    // autodock_ff ligand_type;
    // autodock_ff protein_type;
    fp_type nbp_r   = 0;     // radius of energy-well minimum
    fp_type nbp_eps = 0;     // depth of energy-well minimum
    fp_type xA      = 12;    // exponent for lennard jones, generally 12
    fp_type xB      = 6;     // exponent for lennard jones, 6 for non-hbonders 10 for h-bonders
    bool hbonder    = false; // if the interaction between the two atoms can form an hbond
  };

  // help function that compute the vdw energy shapes for the target types
  // It returns a multi vector where each cell is the energy shape of the interaction between the two types
  // Dimension 0 -> protein types, Dimension 1 -> ligand types
  static inline auto compute_vdw_interaction_shapes() {
    md_vector<vdw_shape, 2> vdw_shapes(num_autodock_ff(), num_autodock_ff_grids());
    for (std::size_t lig_type_index = 0; lig_type_index < vdw_shapes.size<1>(); ++lig_type_index) {
      // const auto lig_type  = target_atom_types[lig_type_index];
      const auto& lig_desc =
          get_description(autodock_ff_from_grid(static_cast<autodock_grid_type>(lig_type_index)));
      for (std::size_t prot_type_index = 0; prot_type_index < vdw_shapes.size<0>(); ++prot_type_index) {
        auto& shape = vdw_shapes.get(prot_type_index, lig_type_index);
        // shape.ligand_type     = target_atom_types[lig_type_index];
        // shape.protein_type    = target_atom_types[prot_type_index];
        const auto& prot_desc = get_description(static_cast<autodock_ff>(prot_type_index));
        shape.nbp_r           = (lig_desc.Rii + prot_desc.Rii) / fp_type{2};
        shape.nbp_eps         = std::sqrt(lig_desc.epsii * autodock_parameters::coeff_vdW * prot_desc.epsii *
                                          autodock_parameters::coeff_vdW);
        shape.hbonder         = lig_desc.hbond > 0 ? true : false;
        if (lig_desc.hbond > 2 && (prot_desc.hbond == 1 || prot_desc.hbond == 2)) {
          shape.xB      = 10;
          shape.nbp_r   = lig_desc.Rij_hb;
          shape.nbp_eps = lig_desc.epsij_hb * autodock_parameters::coeff_hbond;
        } else if ((lig_desc.hbond == 1 || lig_desc.hbond == 2) && prot_desc.hbond > 2) {
          shape.xB      = 10;
          shape.nbp_r   = prot_desc.Rij_hb;
          shape.nbp_eps = prot_desc.epsij_hb * autodock_parameters::coeff_hbond;
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
  static inline auto compute_vdw_interaction_energies() {
    static const md_vector<vdw_shape, 2> vdw_shapes = compute_vdw_interaction_shapes();
    md_vector<fp_type, 3> energy_table(num_radius_tick_desolv, vdw_shapes.size<0>(), vdw_shapes.size<1>());
    for (std::size_t lig_type_index = 0; lig_type_index < vdw_shapes.size<1>(); ++lig_type_index) {
      for (std::size_t prot_type_index = 0; prot_type_index < vdw_shapes.size<0>(); ++prot_type_index) {
        // compute the  energy shape
        const auto& shape     = vdw_shapes.get(prot_type_index, lig_type_index);
        const auto temp_const = shape.nbp_eps / (shape.xA - shape.xB);
        const auto cA         = temp_const * std::pow(shape.nbp_r, shape.xA) * shape.xB;
        const auto cB         = temp_const * std::pow(shape.nbp_r, shape.xB) * shape.xA;
        // auto& voxel           = energy_table.get(0, prot_type_index, lig_type_index);
        // voxel                 = fp_type{100000};
        // for (std::size_t radius_index = 1; radius_index < num_radius_tick_desolv; ++radius_index) {
        //   const auto radius          = static_cast<fp_type>(radius_index) * (1 / lookup_resolution);
        //   const auto rA              = std::pow(radius, );
        //   const auto rB              = std::pow(radius, shape.xB);
        //   const auto proposed_energy = cA / rA - cB / rB;
        //   const auto clapped_energy  = std::min(fp_type{100000}, proposed_energy);
        //   energy_table.get(radius_index, prot_type_index, lig_type_index) = clapped_energy;
        // }

        energy_table.get(0, prot_type_index, lig_type_index) = EINTCLAMP;
        for (int indx_r = 1; indx_r < NEINT - 1; ++indx_r) {
          const fp_type r  = static_cast<fp_type>(indx_r) / A_DIV;
          const fp_type rA = std::pow(r, shape.xA);
          const fp_type rB = std::pow(r, shape.xB);

          energy_table.get(indx_r, prot_type_index, lig_type_index) =
              std::min(EINTCLAMP, (cA / rA - cB / rB));
        }
        energy_table.get(NEINT - 1, prot_type_index, lig_type_index) = 0;

        // apply a "smoothing" process that actually discretize the shape
        // NOTE: the diameter of the smoothing is 0.5A
        // static constexpr auto smooth_diameter = static_cast<std::size_t>(lookup_resolution * fp_type{0.5});
        // static constexpr auto smooth_radius   = smooth_diameter / std::size_t{2};
        // if constexpr (auto* energy_row = &voxel; smooth_diameter < num_radius_tick_desolv) {
        //   for (std::size_t radius_index = 0; radius_index < smooth_radius; radius_index++) {
        //     auto min_value       = energy_row[0];
        //     const auto end_index = radius_index + smooth_radius;
        //     for (std::size_t i = 1; i < end_index; ++i) { min_value = std::min(min_value, energy_row[i]); }
        //     energy_row[radius_index] = min_value;
        //   }
        //   static constexpr auto middle_end = num_radius_tick_desolv - smooth_radius;
        //   for (std::size_t radius_index = smooth_radius; radius_index < middle_end; ++radius_index) {
        //     const auto begin_index   = radius_index - smooth_radius;
        //     const auto end_condition = radius_index + smooth_radius;
        //     auto min_value           = energy_row[begin_index];
        //     for (std::size_t i = begin_index + std::size_t{1}; i < end_condition; ++i) {
        //       min_value = std::min(min_value, energy_row[i]);
        //     }
        //     energy_row[radius_index] = min_value;
        //   }
        //   for (std::size_t radius_index = middle_end; radius_index < num_radius_tick_desolv; ++radius_index) {
        //     const auto begin_index = radius_index - smooth_radius;
        //     auto min_value         = energy_row[begin_index];
        //     for (std::size_t i = begin_index + std::size_t{1}; i < num_radius_tick_desolv; ++i) {
        //       min_value = std::min(min_value, energy_row[i]);
        //     }
        //     energy_row[radius_index] = min_value;
        //   }
        // } else { // we have really low vlaues
        //   const auto min_value = *std::min_element(energy_row, energy_row + num_radius_tick_desolv);
        //   for (std::size_t i = 0; i < num_radius_tick_elect; ++i) { energy_row[i] = min_value; }
        // }
        const int i_smooth = static_cast<int>(std::floor(r_smooth * A_DIV / inv_spacing));
        std::vector<fp_type> energy_smooth;
        energy_smooth.resize(NEINT, EINTCLAMP);
        if (i_smooth > 0) {
          for (int indx_r = 0; indx_r < NEINT; ++indx_r) {
            const auto temp = indx_r - i_smooth;
            int j           = temp <= int{0} ? int{0} : temp;
            for (; j < std::min(int{NEINT}, indx_r + i_smooth + 1); j++)
              energy_smooth[indx_r] =
                  std::min(energy_smooth[indx_r], energy_table.get(j, prot_type_index, lig_type_index));
          }
          for (int indx_r = 0; indx_r < NEINT; indx_r++) {
            energy_table.get(indx_r, prot_type_index, lig_type_index) = energy_smooth[indx_r];
          }
        } /* endif smoothing */
      }
    }
    return energy_table;
  }

  //===------------------------------------------------------------------------------------------------------
  // Utility function to derive the desolvation energy by varying the distance
  //===------------------------------------------------------------------------------------------------------

  // this function will compute the desolvation energy among the two atoms types
  // NOTE: this function is independent from the atom types
  static inline auto compute_desolvation_energy() {
    std::array<fp_type, num_radius_tick_desolv> energy_table;
    for (std::size_t radius_index = 0; radius_index < num_radius_tick_desolv; ++radius_index) {
      const auto radius = static_cast<fp_type>(radius_index) / lookup_resolution;
      energy_table[radius_index] =
          autodock_parameters::coeff_desolv *
          std::exp(-(radius * radius) /
                   (static_cast<fp_type>(2) * static_cast<fp_type>(3.6) * static_cast<fp_type>(3.6)));
    }
    return energy_table;
  }

  //===------------------------------------------------------------------------------------------------------
  // Utility function to compute the dielectric Ewald Summation (Mehler and Solmajer, Prot Eng 4, 903-910)
  //===------------------------------------------------------------------------------------------------------

  static inline auto compute_electostatic_energy() {
    constexpr auto lambda   = static_cast<fp_type>(0.003627);
    constexpr auto epsilon0 = static_cast<fp_type>(78.4);
    constexpr auto A        = static_cast<fp_type>(-8.5525);
    constexpr auto B        = epsilon0 - A;
    constexpr auto rk       = static_cast<fp_type>(7.7839);
    constexpr auto lambda_B = -lambda * B;

    std::array<fp_type, num_radius_tick_elect> adt_protein;
    adt_protein[0] = fp_type{332};
    for (std::size_t radius_index = 1; radius_index < num_radius_tick_elect; ++radius_index) {
      const auto radius         = static_cast<fp_type>(radius_index) / lookup_resolution;
      adt_protein[radius_index] = fp_type{332} / (A + B / (fp_type{1} + rk * std::exp(lambda_B * radius)));
    }
    return adt_protein;
  }

  struct autodock_grid {
    autodock_grid(const point3D min, const point3D max, const fp_type resolution)
        : index(min.difference(max)
                    .apply(std::abs<fp_type>)
                    .divide({resolution})
                    .apply([](fp_type x) { return std::ceil(x); })
                    .add({fp_type{1}})),
          data(index.size_x(), index.size_y(), index.size_z(), num_autodock_grids()),
          _inv_resolution(1 / resolution),
          _min(min),
          _max(max),
          _center{max.difference(min).divide({fp_type{2}}).add(min)} {};
    autodock_grid() {};

    inline auto get_atom_map(const autodock_grid_type type) {
      return space_grid_view<fp_type>{
          _min,
          _max,
          _center,
          _inv_resolution,
          data.get_slice(md_index<4>{index.size_x(), index.size_y(), index.size_z(), static_cast<int>(type)},
                         index)};
    };

    inline auto get_atom_map(const autodock_grid_type type) const {
      return space_grid_view<const fp_type>{
          _min,
          _max,
          _center,
          _inv_resolution,
          data.get_slice(md_index<4>{index.size_x(), index.size_y(), index.size_z(), static_cast<int>(type)},
                         index)};
    };

    [[nodiscard]] inline auto get_map_flat_size() const { return index.flat_size(); };
    [[nodiscard]] inline const auto* get_maps_pointer() const { return data.data(); };
    [[nodiscard]] inline auto* get_min_p() const { return _min.data(); };
    [[nodiscard]] inline auto* get_max_p() const { return _max.data(); };
    [[nodiscard]] inline auto* get_center_p() const { return _center.data(); };
    [[nodiscard]] inline auto get_center() const { return _center; };
    [[nodiscard]] inline auto get_min() const { return _min; };
    [[nodiscard]] inline auto get_max() const { return _max; };

    [[nodiscard]] inline auto get_size_x() const { return index.size_x(); };
    [[nodiscard]] inline auto get_size_xy() const { return index.size_xy(); };
    [[nodiscard]] inline auto get_size_xyz() const { return index.flat_size(); };

    md_index<3> index;
    md_container<std::vector<fp_type>, 4> data;
    fp_type _inv_resolution = 2;
    point<fp_type, 3> _min, _max, _center;
  };

  /**
   * This data structure represents how autodock encode the target protein. The idea is to represent different
   * the 3D space using a set discretized grids. Each grid model one property. We always have two grid to
   * represent electrostatic and desolvation components. Then, we can have a set of grid that represent a
   * contribution for different ligand's atom types.
   */
  struct autodock_protein: public autodock_dynamic_layer {
    inline auto get_atom_map(const autodock_grid_type type) { return adt_grid.get_atom_map(type); };

    inline auto get_atom_map(const autodock_grid_type type) const { return adt_grid.get_atom_map(type); };
    [[nodiscard]] inline auto get_eletrostatic() { return get_atom_map(autodock_grid_type::ELEC); };
    [[nodiscard]] inline auto get_desolvation() { return get_atom_map(autodock_grid_type::DESOLV); };
    [[nodiscard]] inline auto get_map_flat_size() const { return adt_grid.get_map_flat_size(); };
    [[nodiscard]] inline const auto* get_maps_pointer() const { return adt_grid.get_maps_pointer(); };
    [[nodiscard]] inline auto* get_min_p() const { return adt_grid.get_min_p(); };
    [[nodiscard]] inline auto* get_max_p() const { return adt_grid.get_max_p(); };
    [[nodiscard]] inline auto* get_center_p() const { return adt_grid.get_center_p(); };
    [[nodiscard]] inline auto get_center() const { return adt_grid.get_center(); };
    [[nodiscard]] inline auto get_min() const { return adt_grid.get_min(); };
    [[nodiscard]] inline auto get_max() const { return adt_grid.get_max(); };

    [[nodiscard]] inline auto get_size_x() const { return adt_grid.get_size_x(); };
    [[nodiscard]] inline auto get_size_xy() const { return adt_grid.get_size_xy(); };
    [[nodiscard]] inline auto get_size_xyz() const { return adt_grid.get_map_flat_size(); };

    autodock_protein(const point3D min,
                     const point3D max,
                     const fp_type resolution,
                     dynamic_molecule& _molecule,
                     std::function<void(dynamic_molecule&)> f = {})
        : autodock_dynamic_layer(_molecule, f), adt_grid(min, max, resolution) {
      check_source_path_pdbqt_forcefield_override(*this, _molecule);
      autodock_protein::prepare();
    };
    autodock_protein(dynamic_molecule& _molecule, std::function<void(dynamic_molecule&)> f = {})
        : autodock_dynamic_layer(_molecule, f) {
      check_source_path_pdbqt_forcefield_override(*this, _molecule);
      autodock_protein::prepare();
    };
    autodock_protein(dynamic_molecule& _molecule, std::function<void(autodock_dynamic_layer&)> f)
        : autodock_dynamic_layer(_molecule, f) {
      autodock_protein::prepare();
    };
    const std::array<fp_type, num_radius_tick_elect> electrostatic_energies = {compute_electostatic_energy()};
    const std::array<fp_type, num_radius_tick_desolv> desolvation_energies  = {compute_desolvation_energy()};
    const md_vector<fp_type, 3> vdw_energies = {compute_vdw_interaction_energies()};

  private:
    autodock_grid adt_grid;
    void prepare();
  };
} // namespace mudock
