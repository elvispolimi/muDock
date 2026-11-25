
#include <mudock/cpp_implementation/ds_score_cpp.hpp>

namespace mudock {

  inline void calc_drugscore(const int batch_atoms,
                             const int batch_ligands,
                             const int scores_per_ligand,
                             const fp_type* __restrict__ batch_ligand_x,
                             const fp_type* __restrict__ batch_ligand_y,
                             const fp_type* __restrict__ batch_ligand_z,
                             const ds_type* __restrict__ batch_ligand_ds_types,
                             const int* batch_ligand_num_atoms,
                             const fp_type* __restrict__ protein_x,
                             const fp_type* __restrict__ protein_y,
                             const fp_type* __restrict__ protein_z,
                             const ds_type* __restrict__ protein_ds_types,
                             const int protein_num_atoms,
                             const fp_type* __restrict__ ds_scores,
                             fp_type* __restrict__ batch_scores) {
    // Some useful constants
    const unsigned int nd_idx     = ds2index(ds_type::ND);
    const unsigned int out_stride = outer_stride();
    const unsigned int int_stride = inner_stride();
    const unsigned int min_bkt    = min_bucket();
    const unsigned int max_bkt    = max_bucket();

    MUDOCK_CPP_MARKER_REGISTER("CPP ds_score_kernel_core");

    for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
      const int ligand_stride    = ligand_index * batch_atoms * scores_per_ligand;
      const int ligand_num_atoms = batch_ligand_num_atoms[ligand_index];

      const fp_type* __restrict__ pool_ligand_x   = batch_ligand_x + ligand_stride;
      const fp_type* __restrict__ pool_ligand_y   = batch_ligand_y + ligand_stride;
      const fp_type* __restrict__ pool_ligand_z   = batch_ligand_z + ligand_stride;
      const ds_type* __restrict__ ligand_ds_types = batch_ligand_ds_types + ligand_index * batch_atoms;

      fp_type* __restrict__ pool_scores = batch_scores + ligand_index * scores_per_ligand;

      for (int score_index{0}; score_index < scores_per_ligand; ++score_index) {
        const int pool_stride = score_index * batch_atoms;

        const fp_type* __restrict__ ligand_x = pool_ligand_x + pool_stride;
        const fp_type* __restrict__ ligand_y = pool_ligand_y + pool_stride;
        const fp_type* __restrict__ ligand_z = pool_ligand_z + pool_stride;

        MUDOCK_CPP_MARKER_START("CPP ds_score_kernel_core");

        fp_type drugscore{0};

#pragma omp simd
        for (int ligand_it = 0; ligand_it < ligand_num_atoms; ++ligand_it) {
          // Get outer molecule atom type
          const ds_type ligand_ds = ligand_ds_types[ligand_it];
          // Skip non-defined types
          if (ligand_ds == ds_type::ND) {
            continue;
          }

          // Get other outer molecule atom info
          const unsigned int ligand_idx = ds2index(ligand_ds);
          const fp_type ligand_atom_x   = ligand_x[ligand_it];
          const fp_type ligand_atom_y   = ligand_y[ligand_it];
          const fp_type ligand_atom_z   = ligand_z[ligand_it];

          for (int protein_it = 0; protein_it < protein_num_atoms; ++protein_it) {
            // Get inner molecule atom info
            const ds_type protein_ds       = protein_ds_types[protein_it];
            const unsigned int protein_idx = ds2index(protein_ds);
            const fp_type protein_atom_x   = protein_x[protein_it];
            const fp_type protein_atom_y   = protein_y[protein_it];
            const fp_type protein_atom_z   = protein_z[protein_it];

            // Calculate distance between ligand and protein atom
            const fp_type dist_x   = protein_atom_x - ligand_atom_x;
            const fp_type dist_y   = protein_atom_y - ligand_atom_y;
            const fp_type dist_z   = protein_atom_z - ligand_atom_z;
            const fp_type dist_two = dist_x * dist_x + dist_y * dist_y + dist_z * dist_z;
            const fp_type dist     = std::sqrt(dist_two);

            // Bucketize distance and check if in range of scoring dataset
            const unsigned int dist_bucket = bucketize(dist);
            const bool range_valid         = (dist_bucket >= min_bkt) && (dist_bucket <= max_bkt);
            const bool type_valid          = protein_idx != nd_idx;
            const bool valid               = range_valid && type_valid;

            // Get flat index for score lookup
            const unsigned int bucket        = dist_bucket - min_bkt;
            const unsigned int temp_flat_idx = ligand_idx * out_stride + protein_idx * int_stride + bucket;
            const unsigned int flat_idx      = temp_flat_idx * valid; // zero if not valid (safe index)

            // Retrieve score
            const fp_type temp_score = ds_scores[flat_idx];
            const fp_type score      = temp_score * valid; // zero if not valid, (fix for safe index score)

            // if (protein_ds != ds_type::ND)
            //   mudock::info(std::format("MUDOCK - {} - {} - {}",
            //                            get_description(ligand_ds).symbol,
            //                            get_description(protein_ds).symbol,
            //                            dist_bucket));

            drugscore += score;
          }
        }

        MUDOCK_CPP_MARKER_STOP("CPP ds_score_kernel_core");

        pool_scores[score_index] = drugscore;
      }
    }
  }

  inline void calc_drugscore_precomp(const int batch_atoms,
                                     const int batch_ligands,
                                     const int scores_per_ligand,
                                     const fp_type* __restrict__ batch_ligand_x,
                                     const fp_type* __restrict__ batch_ligand_y,
                                     const fp_type* __restrict__ batch_ligand_z,
                                     const ds_type* __restrict__ batch_ligand_ds_types,
                                     const int* batch_ligand_num_atoms,
                                     const fp_type* __restrict__ minimum,
                                     const fp_type* __restrict__ maximum,
                                     const int map_index_x,
                                     const int map_index_xy,
                                     const int map_index_xyz,
                                     const fp_type* __restrict__ grid_maps,
                                     fp_type* __restrict__ batch_scores) {
    MUDOCK_CPP_MARKER_REGISTER("CPP ds_score_kernel_core");

    for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
      const int ligand_stride    = ligand_index * batch_atoms * scores_per_ligand;
      const int ligand_num_atoms = batch_ligand_num_atoms[ligand_index];

      const fp_type* __restrict__ pool_ligand_x   = batch_ligand_x + ligand_stride;
      const fp_type* __restrict__ pool_ligand_y   = batch_ligand_y + ligand_stride;
      const fp_type* __restrict__ pool_ligand_z   = batch_ligand_z + ligand_stride;
      const ds_type* __restrict__ ligand_ds_types = batch_ligand_ds_types + ligand_index * batch_atoms;

      fp_type* __restrict__ pool_scores = batch_scores + ligand_index * scores_per_ligand;

      for (int score_index{0}; score_index < scores_per_ligand; ++score_index) {
        const int pool_stride = score_index * batch_atoms;

        const fp_type* __restrict__ ligand_x = pool_ligand_x + pool_stride;
        const fp_type* __restrict__ ligand_y = pool_ligand_y + pool_stride;
        const fp_type* __restrict__ ligand_z = pool_ligand_z + pool_stride;

        MUDOCK_CPP_MARKER_START("CPP ds_score_kernel_core");

        fp_type drugscore{0};

#pragma omp simd reduction(+ : drugscore)
        for (int ligand_it = 0; ligand_it < ligand_num_atoms; ++ligand_it) {
          // Get ligand atom type
          const ds_type ligand_ds   = ligand_ds_types[ligand_it];
          const bool valid_type     = (ligand_ds != ds_type::ND);
          const unsigned int ds_idx = ds2index(ligand_ds);

          // Get ligand coordinates (load once, aligned access)
          const fp_type ligand_atom_x = ligand_x[ligand_it];
          const fp_type ligand_atom_y = ligand_y[ligand_it];
          const fp_type ligand_atom_z = ligand_z[ligand_it];

          // Compute bounds check as mask (no branches)
          const bool in_bounds_x  = (ligand_atom_x >= minimum[0]) & (ligand_atom_x <= maximum[0]);
          const bool in_bounds_y  = (ligand_atom_y >= minimum[1]) & (ligand_atom_y <= maximum[1]);
          const bool in_bounds_z  = (ligand_atom_z >= minimum[2]) & (ligand_atom_z <= maximum[2]);
          const bool valid_bounds = in_bounds_x & in_bounds_y & in_bounds_z;

          const fp_type grid_coord_x = (ligand_atom_x - minimum[0]) * inv_spacing;
          const fp_type grid_coord_y = (ligand_atom_y - minimum[1]) * inv_spacing;
          const fp_type grid_coord_z = (ligand_atom_z - minimum[2]) * inv_spacing;

          const int atom_index_x  = static_cast<int>(grid_coord_x + fp_type(0.5));
          const fp_type weight_x0 = grid_coord_x - static_cast<fp_type>(atom_index_x);
          const fp_type weight_x1 = fp_type(1) - weight_x0;

          const int atom_index_y  = static_cast<int>(grid_coord_y + fp_type(0.5));
          const fp_type weight_y0 = grid_coord_y - static_cast<fp_type>(atom_index_y);
          const fp_type weight_y1 = fp_type(1) - weight_y0;

          const int atom_index_z  = static_cast<int>(grid_coord_z + fp_type(0.5));
          const fp_type weight_z0 = grid_coord_z - static_cast<fp_type>(atom_index_z);
          const fp_type weight_z1 = fp_type(1) - weight_z0;

          // Compute coefficients
          const fp_type c0 = weight_x1 * weight_y1 * weight_z1;
          const fp_type c1 = weight_x1 * weight_y1 * weight_z0;
          const fp_type c2 = weight_x1 * weight_y0 * weight_z1;
          const fp_type c3 = weight_x1 * weight_y0 * weight_z0;
          const fp_type c4 = weight_x0 * weight_y1 * weight_z1;
          const fp_type c5 = weight_x0 * weight_y1 * weight_z0;
          const fp_type c6 = weight_x0 * weight_y0 * weight_z1;
          const fp_type c7 = weight_x0 * weight_y0 * weight_z0;

          const int base_index = map_index_xy * atom_index_z + atom_index_y * map_index_x + atom_index_x;
          const int map_offset = (base_index + map_index_xyz * ds_idx) * valid_type; // safe-zero if not valid

          // Inline trilinear interpolation for better vectorization
          fp_type value = fp_type(0);
          value += c0 * grid_maps[map_offset];
          value += c1 * grid_maps[map_offset + map_index_xy];
          value += c2 * grid_maps[map_offset + map_index_x];
          value += c3 * grid_maps[map_offset + map_index_x + map_index_xy];
          value += c4 * grid_maps[map_offset + 1];
          value += c5 * grid_maps[map_offset + 1 + map_index_xy];
          value += c6 * grid_maps[map_offset + 1 + map_index_x];
          value += c7 * grid_maps[map_offset + 1 + map_index_x + map_index_xy];

          drugscore += value * (valid_type & valid_bounds);
        }

        MUDOCK_CPP_MARKER_STOP("CPP ds_score_kernel_core");

        pool_scores[score_index] = drugscore;
      }
    }
  }

  template<>
  void ds_score_kernel<queue_cpp>::operator()() {
    q->invoke_kernel<this->ds_region_name>(calc_drugscore,
                                           batch_atoms,
                                           batch_ligands,
                                           scores_per_ligand,
                                           batch_ligand_x,
                                           batch_ligand_y,
                                           batch_ligand_z,
                                           batch_ligand_ds_types,
                                           batch_ligand_num_atoms,
                                           protein_x,
                                           protein_y,
                                           protein_z,
                                           protein_ds_types,
                                           protein_num_atoms,
                                           ds_scores,
                                           scores_b);
  }

  template<>
  void ds_score_precomp_kernel<queue_cpp>::operator()() {
    q->invoke_kernel<this->ds_region_name>(calc_drugscore_precomp,
                                           batch_atoms,
                                           batch_ligands,
                                           scores_per_ligand,
                                           batch_ligand_x,
                                           batch_ligand_y,
                                           batch_ligand_z,
                                           batch_ligand_ds_types,
                                           batch_ligand_num_atoms,
                                           minimum,
                                           maximum,
                                           map_index_x,
                                           map_index_xy,
                                           map_index_xyz,
                                           grid_maps,
                                           scores_b);
  }
} // namespace mudock
