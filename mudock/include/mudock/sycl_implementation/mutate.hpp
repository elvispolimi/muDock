#pragma once

#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>
#include <mudock/utils.hpp>
#include <sycl/sycl.hpp>

#ifndef MUDOCK_SYCL_WG_SIZE
  #define MUDOCK_SYCL_WG_SIZE 32
#endif

namespace mudock {

  template<int MAX_ATOMS>
  inline void translate_molecule_sycl(fp_type* __restrict__ x,
                                      fp_type* __restrict__ y,
                                      fp_type* __restrict__ z,
                                      const fp_type* offset_x,
                                      const fp_type* offset_y,
                                      const fp_type* offset_z,
                                      const int num_atoms,
                                      sycl::nd_item<3> it) {
    MUDOCK_PRAGMA_UNROLL(MUDOCK_UNROLL_FACTOR)
    for (int i = static_cast<int>(it.get_local_id(0)); i < MAX_ATOMS; i += MUDOCK_SYCL_WG_SIZE) {
      if (i < num_atoms) {
        x[i] += *offset_x;
        y[i] += *offset_y;
        z[i] += *offset_z;
      }
    }
  };

  template<int MAX_ATOMS>
  inline void rotate_molecule_sycl(fp_type* __restrict__ x,
                                   fp_type* __restrict__ y,
                                   fp_type* __restrict__ z,
                                   const fp_type* angle_x,
                                   const fp_type* angle_y,
                                   const fp_type* angle_z,
                                   const int num_atoms,
                                   sycl::nd_item<3> it) { // compute the angles sine and cosine
    // compute the molecule center of mass
    const auto& sub_group = it.get_sub_group();

    fp_type c_x{0}, c_y{0}, c_z{0};
    MUDOCK_PRAGMA_UNROLL(MUDOCK_UNROLL_FACTOR)
    for (int i = static_cast<int>(it.get_local_id(0)); i < MAX_ATOMS; i += MUDOCK_SYCL_WG_SIZE) {
      if (i < num_atoms) {
        c_x += x[i];
        c_y += y[i];
        c_z += z[i];
      }
    }
    c_x = sycl::reduce_over_group(sub_group, c_x, sycl::plus<fp_type>()) / static_cast<fp_type>(num_atoms);
    c_y = sycl::reduce_over_group(sub_group, c_y, sycl::plus<fp_type>()) / static_cast<fp_type>(num_atoms);
    c_z = sycl::reduce_over_group(sub_group, c_z, sycl::plus<fp_type>()) / static_cast<fp_type>(num_atoms);

    const auto rad_x = deg_to_rad(*angle_x), rad_y = deg_to_rad(*angle_y), rad_z = deg_to_rad(*angle_z);
    const auto cx = sycl::cos(rad_x), sx = sycl::sin(rad_x);
    const auto cy = sycl::cos(rad_y), sy = sycl::sin(rad_y);
    const auto cz = sycl::cos(rad_z), sz = sycl::sin(rad_z);

    // compute the rotation matrix defined as Rz*Ry*Rx
    const auto m00 = cy * cz;
    const auto m01 = sx * sy * cz - cx * sz;
    const auto m02 = cx * sy * cz + sx * sz;
    const auto m10 = cy * sz;
    const auto m11 = sx * sy * sz + cx * cz;
    const auto m12 = cx * sy * sz - sx * cz;
    const auto m20 = -sy;
    const auto m21 = sx * cy;
    const auto m22 = cx * cy;

    // apply the rotation matrix
    MUDOCK_PRAGMA_UNROLL(MUDOCK_UNROLL_FACTOR)
    for (int i = static_cast<int>(it.get_local_id(0)); i < MAX_ATOMS; i += MUDOCK_SYCL_WG_SIZE) {
      if (i < num_atoms) {
        const auto translated_x = x[i] - c_x, translated_y = y[i] - c_y, translated_z = z[i] - c_z;
        x[i] = translated_x * m00 + translated_y * m01 + translated_z * m02 + c_x;
        y[i] = translated_x * m10 + translated_y * m11 + translated_z * m12 + c_y;
        z[i] = translated_x * m20 + translated_y * m21 + translated_z * m22 + c_z;
      }
    }
  };

  template<int MAX_ATOMS>
  inline void rotate_fragment_sycl(fp_type* __restrict__ x,
                                   fp_type* __restrict__ y,
                                   fp_type* __restrict__ z,
                                   const int* bitmask,
                                   const int start_index,
                                   const int stop_index,
                                   const fp_type* angle,
                                   const int num_atoms,
                                   sycl::nd_item<3> it) { // compute the axis vector (and some properties)
    const auto origx = x[start_index], origy = y[start_index], origz = z[start_index];
    const auto destx = x[stop_index], desty = y[stop_index], destz = z[stop_index];
    const auto u = destx - origx;
    const auto v = desty - origy;
    const auto w = destz - origz;

    const auto u2 = u * u, v2 = v * v, w2 = w * w;
    const auto l2 = u * u + v * v + w * w;
    // Check if origin and dest coincide
    // No need to continue the intramolecular energy will be very high
    // if (std::isinf(l2) || l2 == fp_type{0} || std::isnan(l2))
    //   // TODO print error?
    //   return;
    const auto l = sycl::sqrt(l2);

    // compute the angle sine and cosine
    const auto rad = deg_to_rad(*angle);
    const auto s = sycl::sin(rad), c = sycl::cos(rad);
    const auto one_minus_c = fp_type{1} - c;
    const auto ls          = l * s;

    // compute the rotation matrix (rodrigues' rotation formula)
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

    // apply the rotation matrix
    MUDOCK_PRAGMA_UNROLL(MUDOCK_UNROLL_FACTOR)
    for (int i = static_cast<int>(it.get_local_id(0)); i < MAX_ATOMS; i += MUDOCK_SYCL_WG_SIZE) {
      if (i < num_atoms && bitmask[i] != 0) {
        const auto prev_x = x[i], prev_y = y[i], prev_z = z[i];
        x[i] = prev_x * m00 + prev_y * m01 + prev_z * m02 + m03;
        y[i] = prev_x * m10 + prev_y * m11 + prev_z * m12 + m13;
        z[i] = prev_x * m20 + prev_y * m21 + prev_z * m22 + m23;
      }
    }
  };

  template<int MAX_ATOMS>
  struct apply_sycl {
    void operator()(sycl::nd_item<3> it,
                    const int chromosome_number,
                    const int atom_stride,
                    const fp_type* __restrict__ original_x,
                    const fp_type* __restrict__ original_y,
                    const fp_type* __restrict__ original_z,
                    fp_type* __restrict__ scratch_x,
                    fp_type* __restrict__ scratch_y,
                    fp_type* __restrict__ scratch_z,
                    const chromosome* __restrict__ chromosomes,
                    const int* __restrict__ fragments,
                    const int* __restrict__ ligand_fragments_start,
                    const int* __restrict__ fragments_start_index,
                    const int* __restrict__ fragments_stop_index,
                    const int* __restrict__ frag_indices_start,
                    const int* __restrict__ num_rotamers_b,
                    const int* __restrict__ num_atoms_b) const {
      const int ligand_id       = static_cast<int>(it.get_group(0));
      const int local_thread_id = static_cast<int>(it.get_local_id(0));
      assert(it.get_local_range(0) == MUDOCK_SYCL_WG_SIZE &&
             "SYCL WG size and the number of thread per block does not coincide");

      const int num_atoms    = num_atoms_b[ligand_id];
      const int num_rotamers = num_rotamers_b[ligand_id];

      const fp_type* __restrict__ l_original_x = original_x + ligand_id * atom_stride;
      const fp_type* __restrict__ l_original_y = original_y + ligand_id * atom_stride;
      const fp_type* __restrict__ l_original_z = original_z + ligand_id * atom_stride;
      fp_type* __restrict__ l_scratch_x        = scratch_x + ligand_id * atom_stride * chromosome_number;
      fp_type* __restrict__ l_scratch_y        = scratch_y + ligand_id * atom_stride * chromosome_number;
      fp_type* __restrict__ l_scratch_z        = scratch_z + ligand_id * atom_stride * chromosome_number;
      const chromosome* chromosomes_b          = chromosomes + ligand_id * chromosome_number;
      const auto* __restrict__ l_fragments     = fragments + ligand_fragments_start[ligand_id];
      const auto* __restrict__ l_frag_start_atom_index =
          fragments_start_index + frag_indices_start[ligand_id];
      const auto* __restrict__ l_frag_stop_atom_index = fragments_stop_index + frag_indices_start[ligand_id];

      for (int chromosome_index = 0; chromosome_index < chromosome_number; ++chromosome_index) {
        const chromosome& l_chromosomes            = chromosomes_b[chromosome_index];
        fp_type* __restrict__ x_scratch_chromosome = l_scratch_x + chromosome_index * atom_stride;
        fp_type* __restrict__ y_scratch_chromosome = l_scratch_y + chromosome_index * atom_stride;
        fp_type* __restrict__ z_scratch_chromosome = l_scratch_z + chromosome_index * atom_stride;

        // Copy original coordinates
        MUDOCK_PRAGMA_UNROLL(MUDOCK_UNROLL_FACTOR)
        for (int atom_index = local_thread_id; atom_index < MAX_ATOMS; atom_index += MUDOCK_SYCL_WG_SIZE) {
          if (atom_index < num_atoms) {
            x_scratch_chromosome[atom_index] = l_original_x[atom_index];
            y_scratch_chromosome[atom_index] = l_original_y[atom_index];
            z_scratch_chromosome[atom_index] = l_original_z[atom_index];
          }
        }
        // apply rigid transformations
        translate_molecule_sycl<MAX_ATOMS>(x_scratch_chromosome,
                                           y_scratch_chromosome,
                                           z_scratch_chromosome,
                                           &l_chromosomes[0],
                                           &l_chromosomes[1],
                                           &l_chromosomes[2],
                                           num_atoms,
                                           it);
        rotate_molecule_sycl<MAX_ATOMS>(x_scratch_chromosome,
                                        y_scratch_chromosome,
                                        z_scratch_chromosome,
                                        &l_chromosomes[3],
                                        &l_chromosomes[4],
                                        &l_chromosomes[5],
                                        num_atoms,
                                        it);

        // change the molecule shape
        MUDOCK_PRAGMA_UNROLL(MUDOCK_UNROLL_FACTOR)
        for (int i = 0; i < num_rotamers; ++i) {
          const int* bitmask = l_fragments + i * num_atoms;
          rotate_fragment_sycl<MAX_ATOMS>(x_scratch_chromosome,
                                          y_scratch_chromosome,
                                          z_scratch_chromosome,
                                          bitmask,
                                          l_frag_start_atom_index[i],
                                          l_frag_stop_atom_index[i],
                                          &l_chromosomes[6 + i],
                                          num_atoms,
                                          it);
        }
      }
    }
  };
} // namespace mudock
