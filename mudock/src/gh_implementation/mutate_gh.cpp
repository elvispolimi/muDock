#include <hwy/highway.h>
#include <mudock/gh_implementation/mutate_gh.hpp>
#include <mudock/grid/pi.hpp>

namespace mudock {
  template<>
  void translate_molecule<cpu_vectorization::GH>(fp_type* __restrict__ x,
                                                 fp_type* __restrict__ y,
                                                 fp_type* __restrict__ z,
                                                 const int num_atoms,
                                                 const fp_type offset_x,
                                                 const fp_type offset_y,
                                                 const fp_type offset_z) {
    // Define SIMD type for fp_type (e.g., float or double)
    const HWY_FULL(fp_type) d;
    const auto num_loops = static_cast<size_t>((num_atoms + +Lanes(d) - 1) / Lanes(d));

    // Load offsets as SIMD vectors
    const auto v_offset_x = Set(d, offset_x);
    const auto v_offset_y = Set(d, offset_y);
    const auto v_offset_z = Set(d, offset_z);

    // Process in SIMD lanes
    for (size_t i = 0; i < num_loops * Lanes(d); i += Lanes(d)) {
      const auto remaining = std::min<size_t>(Lanes(d), static_cast<size_t>(num_atoms) - i);
      // Load elements from x, y, and z arrays
      auto vx = LoadN(d, x + i, remaining);
      auto vy = LoadN(d, y + i, remaining);
      auto vz = LoadN(d, z + i, remaining);

      // Add offsets to each component
      vx = Add(vx, v_offset_x);
      vy = Add(vy, v_offset_y);
      vz = Add(vz, v_offset_z);

      // Store results back into the arrays
      StoreN(vx, d, x + i, remaining);
      StoreN(vy, d, y + i, remaining);
      StoreN(vz, d, z + i, remaining);
    }
  }

  template<>
  void rotate_molecule<cpu_vectorization::GH>(fp_type* __restrict__ x,
                                              fp_type* __restrict__ y,
                                              fp_type* __restrict__ z,
                                              const int num_atoms,
                                              const fp_type angle_x,
                                              const fp_type angle_y,
                                              const fp_type angle_z) {
    // Define SIMD type for fp_type (e.g.,x, y, z, num_atoms, angle_x, angle_y, angle_z float or double)
    const HWY_FULL(fp_type) d;
    const auto num_loops = static_cast<size_t>((num_atoms + +Lanes(d) - 1) / Lanes(d));

    // Compute the molecule center of mass
    // Initialize accumulators
    auto sum_x = Zero(d);
    auto sum_y = Zero(d);
    auto sum_z = Zero(d);

    // Process in SIMD lanes
    for (size_t i = 0; i < num_loops * Lanes(d); i += Lanes(d)) {
      const auto remaining = std::min<size_t>(Lanes(d), static_cast<size_t>(num_atoms) - i);
      // Load elements from x, y, and z arrays
      auto vx = LoadN(d, x + i, remaining);
      auto vy = LoadN(d, y + i, remaining);
      auto vz = LoadN(d, z + i, remaining);

      // Accumulate sums
      sum_x = Add(sum_x, vx);
      sum_y = Add(sum_y, vy);
      sum_z = Add(sum_z, vz);
    }

    // Horizontal reduction to compute the final sums
    const auto total_x = ReduceSum(d, sum_x);
    const auto total_y = ReduceSum(d, sum_y);
    const auto total_z = ReduceSum(d, sum_z);
    // Compute center of mass
    const fp_type c_x = total_x / static_cast<fp_type>(num_atoms);
    const fp_type c_y = total_y / static_cast<fp_type>(num_atoms);
    const fp_type c_z = total_z / static_cast<fp_type>(num_atoms);

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

    // Load center of mass as SIMD vectors
    const auto v_c_x = Set(d, c_x);
    const auto v_c_y = Set(d, c_y);
    const auto v_c_z = Set(d, c_z);
    // Factor
    const auto v_m00 = Set(d, m00);
    const auto v_m01 = Set(d, m01);
    const auto v_m02 = Set(d, m02);
    const auto v_m10 = Set(d, m10);
    const auto v_m11 = Set(d, m11);
    const auto v_m12 = Set(d, m12);
    const auto v_m20 = Set(d, m20);
    const auto v_m21 = Set(d, m21);
    const auto v_m22 = Set(d, m22);
    // Process in SIMD lanes
    // TODO check the equal comparison
    for (size_t i = 0; i < num_loops * Lanes(d); i += Lanes(d)) {
      const auto remaining = std::min<size_t>(Lanes(d), static_cast<size_t>(num_atoms) - i);
      // Load elements from x, y, and z arrays
      auto translate_x = LoadN(d, x + i, remaining);
      auto translate_y = LoadN(d, y + i, remaining);
      auto translate_z = LoadN(d, z + i, remaining);

      translate_x = Sub(translate_x, v_c_x);
      translate_y = Sub(translate_y, v_c_y);
      translate_z = Sub(translate_z, v_c_z);

      const auto v_t_x =
          MulAdd(translate_x, v_m00, MulAdd(translate_y, v_m01, MulAdd(translate_z, v_m02, v_c_x)));
      const auto v_t_y =
          MulAdd(translate_x, v_m10, MulAdd(translate_y, v_m11, MulAdd(translate_z, v_m12, v_c_y)));
      const auto v_t_z =
          MulAdd(translate_x, v_m20, MulAdd(translate_y, v_m21, MulAdd(translate_z, v_m22, v_c_z)));

      // Store results back into the arrays
      StoreN(v_t_x, d, x + i, remaining);
      StoreN(v_t_y, d, y + i, remaining);
      StoreN(v_t_z, d, z + i, remaining);
    }
  }

  template<>
  void rotate_fragment<cpu_vectorization::GH>(fp_type* __restrict__ x,
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

    // Define SIMD type for fp_type (e.g., float or double)
    const HWY_FULL(fp_type) d;
    const HWY_FULL(int) d_mask;

    const auto num_loops = static_cast<size_t>((num_atoms + +Lanes(d) - 1) / Lanes(d));
    ;
    // Factor
    const auto v_m00 = Set(d, m00);
    const auto v_m01 = Set(d, m01);
    const auto v_m02 = Set(d, m02);
    const auto v_m03 = Set(d, m03);
    const auto v_m10 = Set(d, m10);
    const auto v_m11 = Set(d, m11);
    const auto v_m12 = Set(d, m12);
    const auto v_m13 = Set(d, m13);
    const auto v_m20 = Set(d, m20);
    const auto v_m21 = Set(d, m21);
    const auto v_m22 = Set(d, m22);
    const auto v_m23 = Set(d, m23);
    // Process in SIMD lanes
    for (size_t i = 0; i < num_loops * Lanes(d); i += Lanes(d)) {
      const auto remaining = std::min<size_t>(Lanes(d), static_cast<size_t>(num_atoms) - i);

      // Load integer mask values
      const auto int_mask_values = LoadN(d_mask, frag_mask + i, remaining);
      // Create mask by comparing to zero (assuming 0 means false, non-zero means true)
      const auto int_mask = Ne(int_mask_values, Zero(d_mask));
      // Convert integer mask to float mask
      const auto m = RebindMask(d, int_mask);

      // Load elements from x, y, and z arrays
      auto v_x = LoadN(d, x + i, remaining);
      auto v_y = LoadN(d, y + i, remaining);
      auto v_z = LoadN(d, z + i, remaining);

      const auto v_t_x = MulAdd(v_x, v_m00, MulAdd(v_y, v_m01, MulAdd(v_z, v_m02, v_m03)));
      const auto v_t_y = MulAdd(v_x, v_m10, MulAdd(v_y, v_m11, MulAdd(v_z, v_m12, v_m13)));
      const auto v_t_z = MulAdd(v_x, v_m20, MulAdd(v_y, v_m21, MulAdd(v_z, v_m22, v_m23)));

      // Store results back into the arrays
      BlendedStore(v_t_x, m, d, x + i);
      BlendedStore(v_t_y, m, d, y + i);
      BlendedStore(v_t_z, m, d, z + i);
    }
  }
} // namespace mudock
