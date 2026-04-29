#include <mudock/chem/mehler_solmajer.hpp>
#include <mudock/compute/adt_score_kernel.hpp>
#include <mudock/compute/devices_memory.hpp>
#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/devices.hpp>
#include <mudock/sycl_implementation/adt_score_sycl.hpp>
#include <mudock/sycl_implementation/invoke_kernel_sycl.hpp>
#include <mudock/sycl_implementation/queue_sycl.hpp>
#include <mudock/sycl_implementation/sycl_texture.hpp>
#include <mudock/sycl_implementation/sycl_utils.hpp>
#include <mudock/utils.hpp>
#include <stdexcept>

#define FLATTENED_3D(x, y, z, index_x, index_xy) (index_xy * (z) + (y) * index_x + (x))

#define BUCKET_MULTIPLIER 3

#ifndef MUDOCK_SYCL_WG_SIZE
  #define MUDOCK_SYCL_WG_SIZE 32
#endif

namespace mudock {
  inline fp_type trilinear_interpolation_sycl(const fp_type* __restrict__ map,
                                              const fp_type* __restrict__ coeffs,
                                              const int& map_index_x,
                                              const int& map_index_xy) {
    fp_type value{0};

    value = coeffs[0] * map[0] + value;
    value = coeffs[1] * map[map_index_xy] + value;
    value = coeffs[2] * map[map_index_x] + value;
    value = coeffs[3] * map[map_index_x + map_index_xy] + value;
    value = coeffs[4] * map[1] + value;
    value = coeffs[5] * map[1 + map_index_xy] + value;
    value = coeffs[6] * map[1 + map_index_x] + value;
    value = coeffs[7] * map[1 + map_index_x + map_index_xy] + value;

    return value;
  }

  constexpr int k_max_devices = 16;
  device_memory_array<k_max_devices, sycl_texture_devices>* get_sycl_texture_memory() {
    // Intentionally leaked to avoid static destruction after SYCL runtime teardown.
    static auto* storage = new device_memory_array<k_max_devices, sycl_texture_devices>();
    return storage;
  }

  void init_device(const int dev,
                   const device_type dev_type,
                   const int map_index_xyz,
                   const fp_type* map_grids) {
    // Thread-safe, exactly-once init per device:
    auto* texture_memory = get_sycl_texture_memory();
    texture_memory->init(dev, dev_type, map_index_xyz, num_autodock_grids(), map_grids);
  }

  template<int MAX_ATOMS>
  struct calc_energy {
    void operator()(sycl::nd_item<3> it,
                    const int atom_stride,
                    const int scores_per_ligand,
                    const fp_type* scratch_x,
                    const fp_type* scratch_y,
                    const fp_type* scratch_z,
                    const fp_type* vol_b,
                    const fp_type* solpar_b,
                    const fp_type* charge_b,
                    const int* num_atoms_b,
                    const int* num_rotamers_b,
                    const int* num_nonbonds_b,
                    const int* __restrict__ nonbond_a1,
                    const int* __restrict__ nonbond_a2,
                    const fp_type* __restrict__ nonbond_cA,
                    const fp_type* __restrict__ nonbond_cB,
                    const int* __restrict__ nonbond_xB,
                    const float* __restrict__ minimum,
                    const float* __restrict__ maximum,
                    const float* __restrict__ center,
                    const int map_index_x,
                    const int map_index_xy,
                    const int map_index_xyz,
                    const fp_type* __restrict__ grid_maps,
                    const int* __restrict__ map_tex_indexes,
                    fp_type* __restrict__ scores) const {
      const int workgroup_id         = static_cast<int>(it.get_group(0));
      const int ligand_id            = static_cast<int>(workgroup_id);
      const int workitem_id_in_group = static_cast<int>(it.get_local_id(0));
      const auto sub_group           = it.get_sub_group();
      assert(it.get_local_range(0) == MUDOCK_SYCL_WG_SIZE &&
             "SYCL WG size and the number of thread per block does not coincide");

      const fp_type* electro_map = grid_maps + map_index_xyz * static_cast<int>(autodock_grid_type::ELEC);
      const fp_type* desolv_map  = grid_maps + map_index_xyz * static_cast<int>(autodock_grid_type::DESOLV);

      const int num_atoms    = num_atoms_b[ligand_id];
      const int num_nonbonds = num_nonbonds_b[ligand_id + 1] - num_nonbonds_b[ligand_id];
      const int num_rotamers = num_rotamers_b[ligand_id];
      const int stride       = ligand_id * atom_stride;

      const fp_type* l_scratch_x = scratch_x + stride * scores_per_ligand;
      const fp_type* l_scratch_y = scratch_y + stride * scores_per_ligand;
      const fp_type* l_scratch_z = scratch_z + stride * scores_per_ligand;
      const fp_type* l_vol       = vol_b + stride;
      const fp_type* l_solpar    = solpar_b + stride;
      const fp_type* l_charge    = charge_b + stride;
      // Point to the next population buffer
      const auto* l_atom_tex_indexes = map_tex_indexes + stride;
      const int* l_nonbond_a1        = nonbond_a1 + num_nonbonds_b[ligand_id];
      const int* l_nonbond_a2        = nonbond_a2 + num_nonbonds_b[ligand_id];
      const fp_type* l_nonbond_cA    = nonbond_cA + num_nonbonds_b[ligand_id];
      const fp_type* l_nonbond_cB    = nonbond_cB + num_nonbonds_b[ligand_id];
      const int* l_nonbond_xB        = nonbond_xB + num_nonbonds_b[ligand_id];

      fp_type* scores_l = scores + ligand_id * scores_per_ligand;

      for (int scores_index = 0; scores_index < scores_per_ligand; ++scores_index) {
        // Copy original coordinates
        const fp_type* ligand_x = l_scratch_x + scores_index * atom_stride;
        const fp_type* ligand_y = l_scratch_y + scores_index * atom_stride;
        const fp_type* ligand_z = l_scratch_z + scores_index * atom_stride;

        // Calculate energy
        fp_type elect_total_trilinear = 0;
        fp_type emap_total_trilinear  = 0;
        fp_type dmap_total_trilinear  = 0;
        MUDOCK_PRAGMA_UNROLL(MUDOCK_UNROLL_FACTOR)
        for (int atom_index = workitem_id_in_group; atom_index < MAX_ATOMS;
             atom_index += MUDOCK_SYCL_WG_SIZE) {
          if (atom_index < num_atoms) {
            fp_type coord_tex[3]{ligand_x[atom_index], ligand_y[atom_index], ligand_z[atom_index]};
            if (coord_tex[0] < minimum[0] || coord_tex[0] > maximum[0] || coord_tex[1] < minimum[1] ||
                coord_tex[1] > maximum[1] || coord_tex[2] < minimum[2] || coord_tex[2] > maximum[2]) {
              // Is outside
              const auto diff_x          = coord_tex[0] - center[0];
              const auto diff_y          = coord_tex[1] - center[1];
              const auto diff_z          = coord_tex[2] - center[2];
              const fp_type distance_two = diff_x * diff_x + diff_y * diff_y + diff_z * diff_z;

              const fp_type epenalty = distance_two * ENERGYPENALTY;
              elect_total_trilinear += epenalty;
              emap_total_trilinear += epenalty;
            } else {
              // Is inside
              // Center atom coordinates on the grid center
              coord_tex[0]            = (coord_tex[0] - minimum[0]) * inv_spacing,
              coord_tex[1]            = (coord_tex[1] - minimum[1]) * inv_spacing;
              coord_tex[2]            = (coord_tex[2] - minimum[2]) * inv_spacing;
              const auto& charge      = l_charge[atom_index];
              const fp_type* atom_map = grid_maps + l_atom_tex_indexes[atom_index];

              const int u0      = static_cast<int>(coord_tex[0]);
              const fp_type p0u = coord_tex[0] - static_cast<fp_type>(u0);
              const fp_type p1u = fp_type{1} - p0u;

              const int v0      = static_cast<int>(coord_tex[1]);
              const fp_type p0v = coord_tex[1] - static_cast<fp_type>(v0);
              const fp_type p1v = fp_type{1} - p0v;

              const int w0      = static_cast<int>(coord_tex[2]);
              const fp_type p0w = coord_tex[2] - static_cast<fp_type>(w0);
              const fp_type p1w = fp_type{1} - p0w;

              const fp_type pu[2] = {p1u, p0u};
              const fp_type pv[2] = {p1v, p0v};
              const fp_type pw[2] = {p1w, p0w};

              // Compute coefficients
              const fp_type coeffs[8] = {pu[0] * pv[0] * pw[0],
                                         pu[0] * pv[0] * pw[1],
                                         pu[0] * pv[1] * pw[0],
                                         pu[0] * pv[1] * pw[1],
                                         pu[1] * pv[0] * pw[0],
                                         pu[1] * pv[0] * pw[1],
                                         pu[1] * pv[1] * pw[0],
                                         pu[1] * pv[1] * pw[1]};

              // Precompute flattened indices
              const int base_index = FLATTENED_3D(u0, v0, w0, map_index_x, map_index_xy);
              // Trilinear Interpolationp
              elect_total_trilinear +=
                  trilinear_interpolation_sycl(electro_map + base_index, coeffs, map_index_x, map_index_xy) *
                  charge;
              emap_total_trilinear +=
                  trilinear_interpolation_sycl(atom_map + base_index, coeffs, map_index_x, map_index_xy);
              dmap_total_trilinear +=
                  trilinear_interpolation_sycl(desolv_map + base_index, coeffs, map_index_x, map_index_xy) *
                  std::fabs(charge);
            }
          }
        }

        fp_type elect_total_eintcal{0}, emap_total_eintcal{0}, dmap_total_eintcal{0};
        if (num_rotamers > 0)
          for (int nonbond_list = workitem_id_in_group; nonbond_list < num_nonbonds;
               nonbond_list += MUDOCK_SYCL_WG_SIZE) {
            const int& a1 = l_nonbond_a1[nonbond_list];
            const int& a2 = l_nonbond_a2[nonbond_list];

            const auto diff_x                = ligand_x[a1] - ligand_x[a2];
            const auto diff_y                = ligand_y[a1] - ligand_y[a2];
            const auto diff_z                = ligand_z[a1] - ligand_z[a2];
            const fp_type distance_two       = diff_x * diff_x + diff_y * diff_y + diff_z * diff_z;
            const fp_type distance_two_clamp = sycl::max(distance_two, RMIN_ELEC_SQUARE);
            const fp_type distance           = sycl::sqrt(distance_two_clamp);

            //  Calculate  Electrostatic  Energy
            const fp_type epsilon =
                mehler_solmajer::A +
                mehler_solmajer::B /
                    (fp_type{1} + mehler_solmajer::rk * sycl::exp(mehler_solmajer::lambda_B * distance));
            const fp_type r_dielectric = fp_type{1} / (distance * epsilon);
            const fp_type e_elec =
                l_charge[a1] * l_charge[a2] * ELECSCALE * autodock_parameters::coeff_estat * r_dielectric;
            elect_total_eintcal += e_elec;

            // Calcuate desolv
            const fp_type nb_desolv = (l_vol[a2] * (l_solpar[a1] + qsolpar * sycl::fabs(l_charge[a1])) +
                                       l_vol[a1] * (l_solpar[a2] + qsolpar * sycl::fabs(l_charge[a2])));

            const fp_type e_desolv = autodock_parameters::coeff_desolv *
                                     sycl::exp(fp_type{-0.5} / (sigma * sigma) * distance_two_clamp) *
                                     nb_desolv;
            dmap_total_eintcal += e_desolv;

            fp_type e_vdW_Hb{0};
            if (distance_two_clamp < nbc2) {
              //  Find internal energy parameters, i.e.  epsilon and r-equilibrium values...
              //  Lennard-Jones and Hydrogen Bond Potentials
              // This can be precomputed as in intnbtable.cc
              const int xA = xA_default;
              const int xB = l_nonbond_xB[nonbond_list];

              if (xA != xB) {
                const fp_type cA = l_nonbond_cA[nonbond_list];
                const fp_type cB = l_nonbond_cB[nonbond_list];

                const auto log_distance = sycl::log(distance);
                const fp_type rA        = sycl::exp(static_cast<fp_type>(xA) * log_distance);
                const fp_type rB        = sycl::exp(static_cast<fp_type>(xB) * log_distance);

                e_vdW_Hb = sycl::min(EINTCLAMP, (cA / rA - cB / rB));
              }
            }
            emap_total_eintcal += e_vdW_Hb;
          }
        fp_type total_energy = emap_total_eintcal + elect_total_eintcal + dmap_total_eintcal +
                               emap_total_trilinear + elect_total_trilinear + dmap_total_trilinear;
        total_energy         = sycl::reduce_over_group(sub_group, total_energy, std::plus<fp_type>());

        if (workitem_id_in_group == 0) {
          const fp_type tors_free_energy =
              static_cast<fp_type>(num_rotamers) * autodock_parameters::coeff_tors;
          scores_l[scores_index] = total_energy + tors_free_energy;
        }
      }
    };
  };

  template<>
  void adt_score_kernel<queue_sycl>::operator()() {
    const int dev_id    = q->get_id();
    const auto dev_type = q->get_dev_type();
    init_device(dev_id, dev_type, map_index_xyz, grid_maps);

    constexpr_switch_bucket<0, reorder_buffer<static_molecule>::get_num_atom_clusters(), 1>(
        [&](const auto atom_index) {
          const auto max_atoms = reorder_buffer<static_molecule>::atoms_clusters[atom_index];
          q->invoke_kernel<calc_energy<max_atoms>>(batch_ligands,
                                                   MUDOCK_SYCL_WG_SIZE,
                                                   batch_atoms,
                                                   scores_per_ligand,
                                                   x_scratch_b,
                                                   y_scratch_b,
                                                   z_scratch_b,
                                                   vols_b,
                                                   solpars_b,
                                                   charges_b,
                                                   num_atoms_b,
                                                   num_rotamers_b,
                                                   num_nonbonds_b,
                                                   nonbond_a1_b,
                                                   nonbond_a2_b,
                                                   nonbond_cA_b,
                                                   nonbond_cB_b,
                                                   nonbond_xB_b,
                                                   minimum,
                                                   maximum,
                                                   center,
                                                   map_index_x,
                                                   map_index_xy,
                                                   map_index_xyz,
                                                   get_sycl_texture_memory()->v[dev_id].data->tex_dev,
                                                   map_offsets_b,
                                                   scores_b);
        },
        batch_atoms,
        reorder_buffer<static_molecule>::atoms_clusters.data());

  }; // namespace mudock

  template<>
  batch_multiple get_adt_score_batch_multiple<queue_sycl>(const int atoms, std::shared_ptr<queue_sycl> q_b) {
    batch_multiple bucket_multiple{};
    constexpr_for<0, reorder_buffer<static_molecule>::get_num_atom_clusters(), 1>([&](const auto atom_index) {
      const auto n_atoms = reorder_buffer<static_molecule>::atoms_clusters[atom_index];
      if (atoms == n_atoms)
        bucket_multiple = get_kernel_batch_multiple_sycl<calc_energy<n_atoms>>(q_b, "adt_score::calc_energy");
    });
    if (bucket_multiple.total_multiple() <= 0)
      throw std::runtime_error(
          "Compilation error: there is a bucket of atoms number which it is not handled.");
    return normalize_batch_multiple(bucket_multiple);
  }
} // namespace mudock
