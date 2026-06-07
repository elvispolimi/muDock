#include <cassert>
#include <mudock/cpp_implementation/center_of_mass.hpp>
#include <mudock/cpp_implementation/mutate_cpp.hpp>
#include <mudock/cpp_implementation/vectorization.hpp>
#include <mudock/grid.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {
  template<>
  void translate_molecule<cpu_vectorization::AUTO>(fp_type* __restrict__ x,
                                                   fp_type* __restrict__ y,
                                                   fp_type* __restrict__ z,
                                                   const int num_atoms,
                                                   const fp_type offset_x,
                                                   const fp_type offset_y,
                                                   const fp_type offset_z) {
#pragma omp simd
    for (int i = 0; i < num_atoms; ++i) {
      x[i] += offset_x;
      y[i] += offset_y;
      z[i] += offset_z;
    }
  }

  template<>
  void rotate_molecule<cpu_vectorization::AUTO>(fp_type* __restrict__ x,
                                                fp_type* __restrict__ y,
                                                fp_type* __restrict__ z,
                                                const int num_atoms,
                                                const fp_type angle_x,
                                                const fp_type angle_y,
                                                const fp_type angle_z) {
    // compute the molecule center of mass
    point3D c{fp_type{0}};
#pragma omp simd
    for (int i = 0; i < num_atoms; i++) {
      c.x() += x[i];
      c.y() += y[i];
      c.z() += z[i];
    }
    c.x() /= static_cast<fp_type>(num_atoms);
    c.y() /= static_cast<fp_type>(num_atoms);
    c.z() /= static_cast<fp_type>(num_atoms);

    // compute the angles sine and cosine
    const auto rad_x = deg_to_rad(angle_x), rad_y = deg_to_rad(angle_y), rad_z = deg_to_rad(angle_z);
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
#pragma omp simd
    for (int i = 0; i < num_atoms; ++i) {
      const auto translated_x = x[i] - c.x(), translated_y = y[i] - c.y(), translated_z = z[i] - c.z();
      x[i] = translated_x * m00 + translated_y * m01 + translated_z * m02 + c.x();
      y[i] = translated_x * m10 + translated_y * m11 + translated_z * m12 + c.y();
      z[i] = translated_x * m20 + translated_y * m21 + translated_z * m22 + c.z();
    }
  }

  template<>
  void rotate_fragment<cpu_vectorization::AUTO>(fp_type* __restrict__ x,
                                                fp_type* __restrict__ y,
                                                fp_type* __restrict__ z,
                                                const int num_atoms,
                                                const int* __restrict__ frag_mask,
                                                const int start_index,
                                                const int stop_index,
                                                const fp_type angle) {
    // compute the axis vector (and some properties)
    const auto origx = x[start_index], origy = y[start_index], origz = z[start_index];
    const auto destx = x[stop_index], desty = y[stop_index], destz = z[stop_index];
    const auto u  = destx - origx;
    const auto v  = desty - origy;
    const auto w  = destz - origz;
    const auto u2 = u * u, v2 = v * v, w2 = w * w;
    const auto l2 = u * u + v * v + w * w;
    const auto l  = std::sqrt(l2);

    // compute the angle sine and cosine
    const auto rad = deg_to_rad(angle);
    const auto s = std::sin(rad), c = std::cos(rad);
    const auto one_minus_c = fp_type{1} - c;
    const auto ls          = l * s;

    // Precompute common sub-expressions to reduce redundant calculations
    const auto inv_l2 = fp_type{1} / l2;
    const auto us_vc  = u * v * one_minus_c;
    const auto uw_vc  = u * w * one_minus_c;
    const auto vw_vc  = v * w * one_minus_c;

    // compute the rotation matrix (rodrigues' rotation formula)
    const auto m00 = (u2 + (v2 + w2) * c) * inv_l2;
    const auto m01 = (us_vc - w * l * s) * inv_l2;
    const auto m02 = (uw_vc + v * l * s) * inv_l2;
    const auto m03 =
        ((origx * (v2 + w2) - u * (origy * v + origz * w)) * one_minus_c + (origy * w - origz * v) * ls) *
        inv_l2;

    const auto m10 = (us_vc + w * ls) * inv_l2;
    const auto m11 = (v2 + (u2 + w2) * c) * inv_l2;
    const auto m12 = (vw_vc - u * ls) * inv_l2;
    const auto m13 =
        ((origy * (u2 + w2) - v * (origx * u + origz * w)) * one_minus_c + (origz * u - origx * w) * ls) *
        inv_l2;

    const auto m20 = (uw_vc - v * ls) * inv_l2;
    const auto m21 = (vw_vc + u * ls) * inv_l2;
    const auto m22 = (w2 + (u2 + v2) * c) * inv_l2;
    const auto m23 =
        ((origz * (u2 + v2) - w * (origx * u + origy * v)) * one_minus_c + (origx * v - origy * u) * ls) *
        inv_l2;

// apply the rotation matrix
#pragma omp simd
    for (int i = 0; i < num_atoms; ++i) {
      if (frag_mask[i] != 0) {
        const auto prev_x = x[i], prev_y = y[i], prev_z = z[i];
        x[i] = prev_x * m00 + prev_y * m01 + prev_z * m02 + m03;
        y[i] = prev_x * m10 + prev_y * m11 + prev_z * m12 + m13;
        z[i] = prev_x * m20 + prev_y * m21 + prev_z * m22 + m23;
      }
    }
  }
} // namespace mudock
