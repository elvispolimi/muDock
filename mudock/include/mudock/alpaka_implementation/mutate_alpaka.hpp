#pragma once

#include <alpaka/alpaka.hpp>
#include <cmath>
#include <mudock/type_alias.hpp>
#include <mudock/utils.hpp>

namespace mudock {

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
    c_x /= static_cast<fp_type>(num_atoms);
    c_y /= static_cast<fp_type>(num_atoms);
    c_z /= static_cast<fp_type>(num_atoms);

    const auto rad_x = deg_to_rad(angle_x), rad_y = deg_to_rad(angle_y), rad_z = deg_to_rad(angle_z);
    const fp_type cx = static_cast<fp_type>(::cos(rad_x));
    const fp_type sx = static_cast<fp_type>(::sin(rad_x));
    const fp_type cy = static_cast<fp_type>(::cos(rad_y));
    const fp_type sy = static_cast<fp_type>(::sin(rad_y));
    const fp_type cz = static_cast<fp_type>(::cos(rad_z));
    const fp_type sz = static_cast<fp_type>(::sin(rad_z));

    const fp_type m00 = cy * cz;
    const fp_type m01 = sx * sy * cz - cx * sz;
    const fp_type m02 = cx * sy * cz + sx * sz;
    const fp_type m10 = cy * sz;
    const fp_type m11 = sx * sy * sz + cx * cz;
    const fp_type m12 = cx * sy * sz - sx * cz;
    const fp_type m20 = -sy;
    const fp_type m21 = sx * cy;
    const fp_type m22 = cx * cy;

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
    const auto l2   = u2 + v2 + w2;
    const fp_type l = static_cast<fp_type>(::sqrt(l2));

    const auto rad            = deg_to_rad(angle);
    const fp_type s           = static_cast<fp_type>(::sin(rad));
    const fp_type c           = static_cast<fp_type>(::cos(rad));
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

    for (int i = 0; i < MAX_ATOMS; ++i) {
      if (i < num_atoms && bitmask[i] != 0) {
        const auto prev_x = x[i], prev_y = y[i], prev_z = z[i];
        x[i] = prev_x * m00 + prev_y * m01 + prev_z * m02 + m03;
        y[i] = prev_x * m10 + prev_y * m11 + prev_z * m12 + m13;
        z[i] = prev_x * m20 + prev_y * m21 + prev_z * m22 + m23;
      }
    }
  }

} // namespace mudock
