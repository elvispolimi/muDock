#pragma once

#include <alpaka/alpaka.hpp>
#include <cmath>
#include <mudock/type_alias.hpp>
#include <mudock/utils.hpp>

namespace mudock {

  template<int MAX_ATOMS, int BLOCK_SIZE, typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE void translate_molecule_alpaka(TAcc const& acc,
                                               fp_type* __restrict__ x,
                                               fp_type* __restrict__ y,
                                               fp_type* __restrict__ z,
                                               const fp_type offset_x,
                                               const fp_type offset_y,
                                               const fp_type offset_z,
                                               const int num_atoms) {
    const int thread_id = static_cast<int>(alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u]);

    ALPAKA_UNROLL(MUDOCK_ATOM_LOOP_UNROLL_FACTOR(MAX_ATOMS, BLOCK_SIZE))
    for (int i = 0; i < MAX_ATOMS; i += BLOCK_SIZE) {
      const int atom_index = i + thread_id;
      if (atom_index < num_atoms) {
        x[atom_index] += offset_x;
        y[atom_index] += offset_y;
        z[atom_index] += offset_z;
      }
    }
  }

  template<int MAX_ATOMS, int BLOCK_SIZE, typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE void rotate_molecule_alpaka(TAcc const& acc,
                                            fp_type* __restrict__ x,
                                            fp_type* __restrict__ y,
                                            fp_type* __restrict__ z,
                                            const fp_type angle_x,
                                            const fp_type angle_y,
                                            const fp_type angle_z,
                                            const int num_atoms) {
    const int thread_id = static_cast<int>(alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u]);

    fp_type c_x{0}, c_y{0}, c_z{0};
    ALPAKA_UNROLL(MUDOCK_ATOM_LOOP_UNROLL_FACTOR(MAX_ATOMS, BLOCK_SIZE))
    for (int i = 0; i < MAX_ATOMS; i += BLOCK_SIZE) {
      const int atom_index = i + thread_id;
      if (atom_index < num_atoms) {
        c_x += x[atom_index];
        c_y += y[atom_index];
        c_z += z[atom_index];
      }
    }

    ALPAKA_UNROLL(MUDOCK_UNROLL_FACTOR)
    for (int offset = BLOCK_SIZE / 2; offset > 0; offset /= 2) {
      c_x += alpaka::warp::shfl_down(acc, c_x, offset, BLOCK_SIZE);
      c_y += alpaka::warp::shfl_down(acc, c_y, offset, BLOCK_SIZE);
      c_z += alpaka::warp::shfl_down(acc, c_z, offset, BLOCK_SIZE);
    }
    c_x = alpaka::warp::shfl(acc, c_x, 0, BLOCK_SIZE);
    c_y = alpaka::warp::shfl(acc, c_y, 0, BLOCK_SIZE);
    c_z = alpaka::warp::shfl(acc, c_z, 0, BLOCK_SIZE);

    c_x /= static_cast<fp_type>(num_atoms);
    c_y /= static_cast<fp_type>(num_atoms);
    c_z /= static_cast<fp_type>(num_atoms);

    const auto rad_x = deg_to_rad(angle_x), rad_y = deg_to_rad(angle_y), rad_z = deg_to_rad(angle_z);

    const fp_type cx = alpaka::math::cos(acc, rad_x);
    const fp_type sx = alpaka::math::sin(acc, rad_x);
    const fp_type cy = alpaka::math::cos(acc, rad_y);
    const fp_type sy = alpaka::math::sin(acc, rad_y);
    const fp_type cz = alpaka::math::cos(acc, rad_z);
    const fp_type sz = alpaka::math::sin(acc, rad_z);

    const fp_type m00 = cy * cz;
    const fp_type m01 = sx * sy * cz - cx * sz;
    const fp_type m02 = cx * sy * cz + sx * sz;
    const fp_type m10 = cy * sz;
    const fp_type m11 = sx * sy * sz + cx * cz;
    const fp_type m12 = cx * sy * sz - sx * cz;
    const fp_type m20 = -sy;
    const fp_type m21 = sx * cy;
    const fp_type m22 = cx * cy;

    ALPAKA_UNROLL(MUDOCK_ATOM_LOOP_UNROLL_FACTOR(MAX_ATOMS, BLOCK_SIZE))
    for (int i = 0; i < MAX_ATOMS; i += BLOCK_SIZE) {
      const int atom_index = i + thread_id;
      if (atom_index < num_atoms) {
        const auto translated_x = x[atom_index] - c_x;
        const auto translated_y = y[atom_index] - c_y;
        const auto translated_z = z[atom_index] - c_z;

        x[atom_index] = translated_x * m00 + translated_y * m01 + translated_z * m02 + c_x;
        y[atom_index] = translated_x * m10 + translated_y * m11 + translated_z * m12 + c_y;
        z[atom_index] = translated_x * m20 + translated_y * m21 + translated_z * m22 + c_z;
      }
    }
  }

  template<int MAX_ATOMS, int BLOCK_SIZE, typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE void rotate_fragment_alpaka(TAcc const& acc,
                                            fp_type* __restrict__ x,
                                            fp_type* __restrict__ y,
                                            fp_type* __restrict__ z,
                                            const int* __restrict__ bitmask,
                                            const int start_index,
                                            const int stop_index,
                                            const fp_type angle,
                                            const int num_atoms) {
    const int thread_id = static_cast<int>(alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u]);

    const auto origx = x[start_index], origy = y[start_index], origz = z[start_index];
    const auto destx = x[stop_index], desty = y[stop_index], destz = z[stop_index];
    const auto u = destx - origx;
    const auto v = desty - origy;
    const auto w = destz - origz;

    const auto u2 = u * u, v2 = v * v, w2 = w * w;
    const auto l2   = u2 + v2 + w2;
    const fp_type l = alpaka::math::sqrt(acc, l2);

    const auto rad            = deg_to_rad(angle);
    const fp_type s           = alpaka::math::sin(acc, rad);
    const fp_type c           = alpaka::math::cos(acc, rad);
    const fp_type one_minus_c = fp_type{1} - c;
    const fp_type ls          = l * s;

    const fp_type m00 = (u2 + (v2 + w2) * c) / l2;
    const fp_type m01 = (u * v * one_minus_c - w * l * s) / l2;
    const fp_type m02 = (u * w * one_minus_c + v * l * s) / l2;
    const fp_type m03 =
        ((origx * (v2 + w2) - u * (origy * v + origz * w)) * one_minus_c + (origy * w - origz * v) * ls) / l2;
    const fp_type m10 = (u * v * one_minus_c + w * ls) / l2;
    const fp_type m11 = (v2 + (u2 + w2) * c) / l2;
    const fp_type m12 = (v * w * one_minus_c - u * ls) / l2;
    const fp_type m13 =
        ((origy * (u2 + w2) - v * (origx * u + origz * w)) * one_minus_c + (origz * u - origx * w) * ls) / l2;
    const fp_type m20 = (u * w * one_minus_c - v * ls) / l2;
    const fp_type m21 = (v * w * one_minus_c + u * ls) / l2;
    const fp_type m22 = (w2 + (u2 + v2) * c) / l2;
    const fp_type m23 =
        ((origz * (u2 + v2) - w * (origx * u + origy * v)) * one_minus_c + (origx * v - origy * u) * ls) / l2;

    ALPAKA_UNROLL(MUDOCK_ATOM_LOOP_UNROLL_FACTOR(MAX_ATOMS, BLOCK_SIZE))
    for (int i = 0; i < MAX_ATOMS; i += BLOCK_SIZE) {
      const int atom_index = i + thread_id;
      if (atom_index < num_atoms && bitmask[atom_index] != 0) {
        const auto prev_x = x[atom_index], prev_y = y[atom_index], prev_z = z[atom_index];
        x[atom_index] = prev_x * m00 + prev_y * m01 + prev_z * m02 + m03;
        y[atom_index] = prev_x * m10 + prev_y * m11 + prev_z * m12 + m13;
        z[atom_index] = prev_x * m20 + prev_y * m21 + prev_z * m22 + m23;
      }
    }
  }

} // namespace mudock
