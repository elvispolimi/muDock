#include <mudock/alpaka_implementation/geom_transform_alpaka.hpp>

#include <alpaka/alpaka.hpp>

#include <cmath>
#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/utils.hpp>

#ifndef MUDOCK_ALPAKA_BLOCK_SIZE
  #define MUDOCK_ALPAKA_BLOCK_SIZE 32
#endif

namespace mudock {
  namespace {
    template<int MAX_ATOMS>
    ALPAKA_FN_ACC void translate_molecule_alpaka(fp_type* x,
                                                 fp_type* y,
                                                 fp_type* z,
                                                 const fp_type offset_x,
                                                 const fp_type offset_y,
                                                 const fp_type offset_z,
                                                 const int num_atoms) {
      for (int i = 0; i < MAX_ATOMS; ++i) {
        if (i < num_atoms) {
          x[i] += offset_x;
          y[i] += offset_y;
          z[i] += offset_z;
        }
      }
    }

    template<int MAX_ATOMS>
    ALPAKA_FN_ACC void rotate_molecule_alpaka(fp_type* x,
                                              fp_type* y,
                                              fp_type* z,
                                              const fp_type angle_x,
                                              const fp_type angle_y,
                                              const fp_type angle_z,
                                              const int num_atoms) {
      fp_type c_x{0}, c_y{0}, c_z{0};
      for (int i = 0; i < MAX_ATOMS; ++i) {
        if (i < num_atoms) {
          c_x += x[i];
          c_y += y[i];
          c_z += z[i];
        }
      }
      c_x /= num_atoms;
      c_y /= num_atoms;
      c_z /= num_atoms;

      const auto rad_x = deg_to_rad(angle_x), rad_y = deg_to_rad(angle_y), rad_z = deg_to_rad(angle_z);
      const auto cx = ::cos(rad_x), sx = ::sin(rad_x);
      const auto cy = ::cos(rad_y), sy = ::sin(rad_y);
      const auto cz = ::cos(rad_z), sz = ::sin(rad_z);

      const auto m00 = cy * cz;
      const auto m01 = sx * sy * cz - cx * sz;
      const auto m02 = cx * sy * cz + sx * sz;
      const auto m10 = cy * sz;
      const auto m11 = sx * sy * sz + cx * cz;
      const auto m12 = cx * sy * sz - sx * cz;
      const auto m20 = -sy;
      const auto m21 = sx * cy;
      const auto m22 = cx * cy;

      for (int i = 0; i < MAX_ATOMS; ++i) {
        if (i < num_atoms) {
          const auto translated_x = x[i] - c_x, translated_y = y[i] - c_y, translated_z = z[i] - c_z;
          x[i] = translated_x * m00 + translated_y * m01 + translated_z * m02 + c_x;
          y[i] = translated_x * m10 + translated_y * m11 + translated_z * m12 + c_y;
          z[i] = translated_x * m20 + translated_y * m21 + translated_z * m22 + c_z;
        }
      }
    }

    template<int MAX_ATOMS>
    ALPAKA_FN_ACC void rotate_fragment_alpaka(fp_type* x,
                                              fp_type* y,
                                              fp_type* z,
                                              const int* bitmask,
                                              const int start_index,
                                              const int stop_index,
                                              const fp_type angle,
                                              const int num_atoms) {
      const auto origx = x[start_index], origy = y[start_index], origz = z[start_index];
      const auto destx = x[stop_index], desty = y[stop_index], destz = z[stop_index];
      const auto u = destx - origx;
      const auto v = desty - origy;
      const auto w = destz - origz;

      const auto u2 = u * u, v2 = v * v, w2 = w * w;
      const auto l2 = u2 + v2 + w2;
      const auto l = ::sqrt(l2);

      const auto rad = deg_to_rad(angle);
      const auto s = ::sin(rad), c = ::cos(rad);
      const auto one_minus_c = fp_type{1} - c;
      const auto ls          = l * s;

      const auto m00 = (u2 + (v2 + w2) * c) / l2;
      const auto m01 = (u * v * one_minus_c - w * l * s) / l2;
      const auto m02 = (u * w * one_minus_c + v * l * s) / l2;
      const auto m03 =
          ((origx * (v2 + w2) - u * (origy * v + origz * w)) * one_minus_c + (origy * w - origz * v) * ls) / l2;

      const auto m10 = (u * v * one_minus_c + w * ls) / l2;
      const auto m11 = (v2 + (u2 + w2) * c) / l2;
      const auto m12 = (v * w * one_minus_c - u * ls) / l2;
      const auto m13 =
          ((origy * (u2 + w2) - v * (origx * u + origz * w)) * one_minus_c + (origz * u - origx * w) * ls) / l2;

      const auto m20 = (u * w * one_minus_c - v * ls) / l2;
      const auto m21 = (v * w * one_minus_c + u * ls) / l2;
      const auto m22 = (w2 + (u2 + v2) * c) / l2;
      const auto m23 =
          ((origz * (u2 + v2) - w * (origx * u + origy * v)) * one_minus_c + (origx * v - origy * u) * ls) / l2;

      for (int i = 0; i < MAX_ATOMS; ++i) {
        if (i < num_atoms && bitmask[i] != 0) {
          const auto prev_x = x[i], prev_y = y[i], prev_z = z[i];
          x[i] = prev_x * m00 + prev_y * m01 + prev_z * m02 + m03;
          y[i] = prev_x * m10 + prev_y * m11 + prev_z * m12 + m13;
          z[i] = prev_x * m20 + prev_y * m21 + prev_z * m22 + m23;
        }
      }
    }

    template<int MAX_ATOMS>
    struct apply_alpaka {
      template<typename TAcc>
      ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                    const int chromosome_number,
                                    const int atom_stride,
                                    const fp_type* original_x,
                                    const fp_type* original_y,
                                    const fp_type* original_z,
                                    fp_type* scratch_x,
                                    fp_type* scratch_y,
                                    fp_type* scratch_z,
                                    const chromosome* chromosomes,
                                    const int* fragments,
                                    const int* ligand_fragments_start,
                                    const int* fragments_start_index,
                                    const int* fragments_stop_index,
                                    const int* frag_indices_start,
                                    const int* num_rotamers_b,
                                    const int* num_atoms_b) const {
        const int ligand_id = static_cast<int>(alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u]);
        const int thread_id = static_cast<int>(alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u]);
        if (thread_id != 0) {
          return;
        }

        const int num_atoms    = num_atoms_b[ligand_id];
        const int num_rotamers = num_rotamers_b[ligand_id];

        const fp_type* l_original_x = original_x + ligand_id * atom_stride;
        const fp_type* l_original_y = original_y + ligand_id * atom_stride;
        const fp_type* l_original_z = original_z + ligand_id * atom_stride;
        fp_type* l_scratch_x        = scratch_x + ligand_id * atom_stride * chromosome_number;
        fp_type* l_scratch_y        = scratch_y + ligand_id * atom_stride * chromosome_number;
        fp_type* l_scratch_z        = scratch_z + ligand_id * atom_stride * chromosome_number;
        const chromosome* chromosomes_b = chromosomes + ligand_id * chromosome_number;
        const auto* l_fragments         = fragments + ligand_fragments_start[ligand_id];
        const auto* l_frag_start_atom_index = fragments_start_index + frag_indices_start[ligand_id];
        const auto* l_frag_stop_atom_index  = fragments_stop_index + frag_indices_start[ligand_id];

        for (int chromosome_index = 0; chromosome_index < chromosome_number; ++chromosome_index) {
          const chromosome& l_chromosomes = chromosomes_b[chromosome_index];
          fp_type* x_scratch_chromosome   = l_scratch_x + chromosome_index * atom_stride;
          fp_type* y_scratch_chromosome   = l_scratch_y + chromosome_index * atom_stride;
          fp_type* z_scratch_chromosome   = l_scratch_z + chromosome_index * atom_stride;

          for (int atom_index = 0; atom_index < MAX_ATOMS; ++atom_index) {
            if (atom_index < num_atoms) {
              x_scratch_chromosome[atom_index] = l_original_x[atom_index];
              y_scratch_chromosome[atom_index] = l_original_y[atom_index];
              z_scratch_chromosome[atom_index] = l_original_z[atom_index];
            }
          }

          translate_molecule_alpaka<MAX_ATOMS>(x_scratch_chromosome,
                                               y_scratch_chromosome,
                                               z_scratch_chromosome,
                                               l_chromosomes[0],
                                               l_chromosomes[1],
                                               l_chromosomes[2],
                                               num_atoms);
          rotate_molecule_alpaka<MAX_ATOMS>(x_scratch_chromosome,
                                            y_scratch_chromosome,
                                            z_scratch_chromosome,
                                            l_chromosomes[3],
                                            l_chromosomes[4],
                                            l_chromosomes[5],
                                            num_atoms);

          for (int i = 0; i < num_rotamers; ++i) {
            const int* bitmask = l_fragments + i * num_atoms;
            rotate_fragment_alpaka<MAX_ATOMS>(x_scratch_chromosome,
                                              y_scratch_chromosome,
                                              z_scratch_chromosome,
                                              bitmask,
                                              l_frag_start_atom_index[i],
                                              l_frag_stop_atom_index[i],
                                              l_chromosomes[6 + i],
                                              num_atoms);
          }
        }
      }
    };
  } // namespace

  template<>
  void geom_kernel<queue_alpaka>::operator()() {
    constexpr_switch_bucket<0, reorder_buffer<static_molecule>::get_num_atom_clusters(), 1>(
        [&](const auto atom_index) {
          const auto max_atoms = reorder_buffer<static_molecule>::atoms_clusters[atom_index];
          q->invoke_kernel<apply_alpaka<max_atoms>>(batch_ligands,
                                                    MUDOCK_ALPAKA_BLOCK_SIZE,
                                                    chromsomes_per_ligand,
                                                    batch_atoms,
                                                    x_coords_b,
                                                    y_coords_b,
                                                    z_coords_b,
                                                    x_scratch_b,
                                                    y_scratch_b,
                                                    z_scratch_b,
                                                    chromosomes_b,
                                                    ligand_fragments_b,
                                                    ligand_fragments_start_b,
                                                    frag_start_indices_b,
                                                    frag_stop_indices_b,
                                                    frag_indices_start_b,
                                                    num_rotamers_b,
                                                    num_atoms_b);
        },
        batch_atoms,
        reorder_buffer<static_molecule>::atoms_clusters.data());
  }
} // namespace mudock
