#pragma once

#include <mudock/grid/point3D.hpp>
#include <mudock/molecule.hpp>
#include <mudock/molecule/containers.hpp>
#include <mudock/type_alias.hpp>
#include <span>
#include <sycl/sycl.hpp>

namespace mudock {

  void translate_molecule_sycl(fp_type* __restrict__ x,
                               fp_type* __restrict__ y,
                               fp_type* __restrict__ z,
                               const fp_type* offset_x,
                               const fp_type* offset_y,
                               const fp_type* offset_z,
                               const int num_atoms,
                               sycl::nd_item<1> it) {
    for (int i = it.get_local_id(0); i < num_atoms; i += it.get_local_range(0)) {
      x[i] += *offset_x;
      y[i] += *offset_y;
      z[i] += *offset_z;
    }
  };

  void rotate_molecule_sycl(fp_type* __restrict__ x,
                            fp_type* __restrict__ y,
                            fp_type* __restrict__ z,
                            const fp_type* angle_x,
                            const fp_type* angle_y,
                            const fp_type* angle_z,
                            const int num_atoms,
                            sycl::nd_item<1> it) { // compute the angles sine and cosine
    // compute the molecule center of mass
    const auto& sub_group = it.get_sub_group();

    fp_type c_x{0}, c_y{0}, c_z{0};
    for (int i = it.get_local_id(0); i < num_atoms; i += it.get_local_range(0)) {
      c_x += x[i];
      c_y += y[i];
      c_z += z[i];
    }
    c_x = sycl::reduce_over_group(sub_group, c_x, sycl::plus<fp_type>()) / num_atoms;
    c_y = sycl::reduce_over_group(sub_group, c_x, sycl::plus<fp_type>()) / num_atoms;
    c_z = sycl::reduce_over_group(sub_group, c_x, sycl::plus<fp_type>()) / num_atoms;

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
    for (int i = it.get_local_id(0); i < num_atoms; i += it.get_local_range(0)) {
      const auto translated_x = x[i] - c_x, translated_y = y[i] - c_y, translated_z = z[i] - c_z;
      x[i] = translated_x * m00 + translated_y * m01 + translated_z * m02 + c_x;
      y[i] = translated_x * m10 + translated_y * m11 + translated_z * m12 + c_y;
      z[i] = translated_x * m20 + translated_y * m21 + translated_z * m22 + c_z;
    }
  };

  void rotate_fragment_sycl(fp_type* __restrict__ x,
                            fp_type* __restrict__ y,
                            fp_type* __restrict__ z,
                            const int* bitmask,
                            const int start_index,
                            const int stop_index,
                            const fp_type* angle,
                            const int num_atoms,
                            sycl::nd_item<1> it) { // compute the axis vector (and some properties)
    const auto origx = x[start_index], origy = y[start_index], origz = z[start_index];
    const auto destx = x[stop_index], desty = y[stop_index], destz = z[stop_index];
    const auto u = destx - origx;
    const auto v = desty - origy;
    const auto w = destz - origz;

    const auto u2 = u * u, v2 = v * v, w2 = w * w;
    const auto l2 = u * u + v * v + w * w;
    // Check if origin and dest coincide
    // No need to continue the intramolecular energy will be very high
    if (std::isinf(l2) || l2 == fp_type{0} || std::isnan(l2))
      // TODO print error?
      return;
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
    for (int i = it.get_local_id(0); i < num_atoms; i += it.get_local_range(0)) {
      if (bitmask[i] != 0) {
        const auto prev_x = x[i], prev_y = y[i], prev_z = z[i];
        x[i] = prev_x * m00 + prev_y * m01 + prev_z * m02 + m03;
        y[i] = prev_x * m10 + prev_y * m11 + prev_z * m12 + m13;
        z[i] = prev_x * m20 + prev_y * m21 + prev_z * m22 + m23;
      }
    }
  };

} // namespace mudock
