#pragma once

#include <mudock/grid.hpp>
#include <mudock/molecule.hpp>
#include <mudock/molecule/containers.hpp>
#include <mudock/type_alias.hpp>
#include <span>

namespace mudock {

  __device__ inline void translate_molecule_cuda(fp_type* __restrict__ x,
                                          fp_type* __restrict__ y,
                                          fp_type* __restrict__ z,
                                          const fp_type* offset_x,
                                          const fp_type* offset_y,
                                          const fp_type* offset_z,
                                          const int num_atoms) {
    for (int i = threadIdx.x; i < num_atoms; i += blockDim.x) {
      x[i] += *offset_x;
      y[i] += *offset_y;
      z[i] += *offset_z;
    }
  }

  __device__ inline void rotate_molecule_cuda(fp_type* __restrict__ x,
                                       fp_type* __restrict__ y,
                                       fp_type* __restrict__ z,
                                       const fp_type* angle_x,
                                       const fp_type* angle_y,
                                       const fp_type* angle_z,
                                       const int num_atoms) {
    // compute the angles sine and cosine
    const auto rad_x = deg_to_rad(*angle_x), rad_y = deg_to_rad(*angle_y), rad_z = deg_to_rad(*angle_z);
    const auto cx = std::cos(rad_x), sx = std::sin(rad_x);
    const auto cy = std::cos(rad_y), sy = std::sin(rad_y);
    const auto cz = std::cos(rad_z), sz = std::sin(rad_z);

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
    for (int i = threadIdx.x; i < num_atoms; i += blockDim.x) {
      const auto prev_x = x[i], prev_y = y[i], prev_z = z[i];
      x[i] = prev_x * m00 + prev_y * m01 + prev_z * m02;
      y[i] = prev_x * m10 + prev_y * m11 + prev_z * m12;
      z[i] = prev_x * m20 + prev_y * m21 + prev_z * m22;
    }
  }

  __device__ inline void rotate_fragment_cuda(fp_type* __restrict__ x,
                                       fp_type* __restrict__ y,
                                       fp_type* __restrict__ z,
                                       const int* bitmask,
                                       const int start_index,
                                       const int stop_index,
                                       const fp_type* angle,
                                       const int num_atoms) {
    // compute the axis vector (and some properties)
    const auto origx = x[start_index], origy = y[start_index], origz = z[start_index];
    const auto destx = x[stop_index], desty = y[stop_index], destz = z[stop_index];
    const auto u = destx - origx;
    const auto v = desty - origy;
    const auto w = destz - origz;

    const auto u2 = u * u, v2 = v * v, w2 = w * w;
    const auto l2 = u * u + v * v + w * w;
    // Check if origin and dest coincide
    // No need to continue the intramolecular energy will be very high
    if (isinf(l2) || l2 == fp_type{0} || isnan(l2))
      // TODO print error?
      return;
    const auto l = std::sqrt(l2);

    // compute the angle sine and cosine
    const auto rad = deg_to_rad(*angle);
    const auto s = std::sin(rad), c = std::cos(rad);
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
    for (int i = threadIdx.x; i < num_atoms; i += blockDim.x) {
      if (bitmask[i] == 1) {
        const auto prev_x = x[i], prev_y = y[i], prev_z = z[i];
        x[i] = prev_x * m00 + prev_y * m01 + prev_z * m02 + m03;
        y[i] = prev_x * m10 + prev_y * m11 + prev_z * m12 + m13;
        z[i] = prev_x * m20 + prev_y * m21 + prev_z * m22 + m23;
      }
    }
  }
} // namespace mudock
