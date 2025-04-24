#include <mudock/cpp_implementation/geometric_transformations_xsimd.hpp>
#include <mudock/grid/pi.hpp>
#include <mudock/type_alias.hpp>
#include <xsimd/xsimd.hpp>

namespace mudock {
  // FIXME unaligned memory requires fixing memory allocation in virtual screening
  // FIXME check tail processing of last elements of the arrays (GH exposes loadN methods)
  template<>
  void translate_molecule<cpu_vectorization::XSIMD>(fp_type* __restrict__ x,
                                                    fp_type* __restrict__ y,
                                                    fp_type* __restrict__ z,
                                                    const int num_atoms,
                                                    const fp_type offset_x,
                                                    const fp_type offset_y,
                                                    const fp_type offset_z) {
    // Define SIMD batch type for fp_type
    using batch_type = xsimd::batch<fp_type>;

    // Get the SIMD batch size
    constexpr auto simd_size = batch_type::size;

    // Load offsets as SIMD vectors
    const batch_type v_offset_x(offset_x);
    const batch_type v_offset_y(offset_y);
    const batch_type v_offset_z(offset_z);

    // Process in SIMD lanes
    size_t i = 0;
    for (; i + simd_size <= static_cast<size_t>(num_atoms); i += simd_size) {
      // Load elements from x, y, and z arrays
      auto vx = batch_type::load_unaligned(x + i);
      auto vy = batch_type::load_unaligned(y + i);
      auto vz = batch_type::load_unaligned(z + i);

      // Add offsets to each component
      vx += v_offset_x;
      vy += v_offset_y;
      vz += v_offset_z;

      // Store results back
      vx.store_unaligned(x + i);
      vy.store_unaligned(y + i);
      vz.store_unaligned(z + i);
    }

    // Process remaining elements
    for (; i < static_cast<size_t>(num_atoms); ++i) {
      x[i] += offset_x;
      y[i] += offset_y;
      z[i] += offset_z;
    }
  }

  template<>
  void rotate_molecule<cpu_vectorization::XSIMD>(fp_type* __restrict__ x,
                                                 fp_type* __restrict__ y,
                                                 fp_type* __restrict__ z,
                                                 const int num_atoms,
                                                 const fp_type angle_x,
                                                 const fp_type angle_y,
                                                 const fp_type angle_z) {
    // Define SIMD batch type
    using batch_type           = xsimd::batch<fp_type>;
    constexpr size_t simd_size = batch_type::size;

    // Compute center of mass
    batch_type sum_x(0), sum_y(0), sum_z(0);
    fp_type total_x = 0, total_y = 0, total_z = 0;

    // Process full SIMD lanes
    size_t i = 0;
    for (; i + simd_size <= static_cast<size_t>(num_atoms); i += simd_size) {
      const auto vx = batch_type::load_unaligned(x + i);
      const auto vy = batch_type::load_unaligned(y + i);
      const auto vz = batch_type::load_unaligned(z + i);

      sum_x += vx;
      sum_y += vy;
      sum_z += vz;
    }

    // Horizontal reduction
    auto partial_sums = xsimd::reduce_add(sum_x);
    total_x += partial_sums;
    partial_sums = xsimd::reduce_add(sum_y);
    total_y += partial_sums;
    partial_sums = xsimd::reduce_add(sum_z);
    total_z += partial_sums;

    // Process remaining elements
    for (; i < static_cast<size_t>(num_atoms); ++i) {
      total_x += x[i];
      total_y += y[i];
      total_z += z[i];
    }

    // Compute center of mass
    const fp_type c_x = total_x / num_atoms;
    const fp_type c_y = total_y / num_atoms;
    const fp_type c_z = total_z / num_atoms;

    // Compute rotation matrix
    const auto rad_x = deg_to_rad(angle_x), rad_y = deg_to_rad(angle_y), rad_z = deg_to_rad(angle_z);
    const auto cx = std::cos(rad_x), sx = std::sin(rad_x);
    const auto cy = std::cos(rad_y), sy = std::sin(rad_y);
    const auto cz = std::cos(rad_z), sz = std::sin(rad_z);

    const auto m00 = cy * cz;
    const auto m01 = sx * sy * cz - cx * sz;
    const auto m02 = cx * sy * cz + sx * sz;
    const auto m10 = cy * sz;
    const auto m11 = sx * sy * sz + cx * cz;
    const auto m12 = cx * sy * sz - sx * cz;
    const auto m20 = -sy;
    const auto m21 = sx * cy;
    const auto m22 = cx * cy;

    // Create SIMD constants
    const batch_type v_c_x(c_x), v_c_y(c_y), v_c_z(c_z);
    const batch_type v_m00(m00), v_m01(m01), v_m02(m02);
    const batch_type v_m10(m10), v_m11(m11), v_m12(m12);
    const batch_type v_m20(m20), v_m21(m21), v_m22(m22);

    // Apply rotation
    i = 0;
    for (; i + simd_size <= static_cast<size_t>(num_atoms); i += simd_size) {
      const auto vx = batch_type::load_unaligned(x + i);
      const auto vy = batch_type::load_unaligned(y + i);
      const auto vz = batch_type::load_unaligned(z + i);

      // Center coordinates
      const auto tx = vx - v_c_x;
      const auto ty = vy - v_c_y;
      const auto tz = vz - v_c_z;

      // Matrix multiplication (fused multiply-add)
      const auto rx = xsimd::fma(tx, v_m00, xsimd::fma(ty, v_m01, xsimd::fma(tz, v_m02, v_c_x)));
      const auto ry = xsimd::fma(tx, v_m10, xsimd::fma(ty, v_m11, xsimd::fma(tz, v_m12, v_c_y)));
      const auto rz = xsimd::fma(tx, v_m20, xsimd::fma(ty, v_m21, xsimd::fma(tz, v_m22, v_c_z)));

      rx.store_unaligned(x + i);
      ry.store_unaligned(y + i);
      rz.store_unaligned(z + i);
    }

    // Process remaining elements
    for (; i < static_cast<size_t>(num_atoms); ++i) {
      // Center coordinates
      const fp_type tx = x[i] - c_x;
      const fp_type ty = y[i] - c_y;
      const fp_type tz = z[i] - c_z;

      // Apply rotation
      x[i] = m00 * tx + m01 * ty + m02 * tz + c_x;
      y[i] = m10 * tx + m11 * ty + m12 * tz + c_y;
      z[i] = m20 * tx + m21 * ty + m22 * tz + c_z;
    }
  }

  template<>
  void rotate_fragment<cpu_vectorization::XSIMD>(fp_type* __restrict__ x,
                                                 fp_type* __restrict__ y,
                                                 fp_type* __restrict__ z,
                                                 const int num_atoms,
                                                 const int* __restrict__ frag_mask,
                                                 const int start_index,
                                                 const int stop_index,
                                                 const fp_type angle) {
    using float_batch          = xsimd::batch<fp_type>;
    using int_batch            = xsimd::batch<int>;
    constexpr size_t simd_size = float_batch::size;

    // Compute rotation axis and properties
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

    // Create SIMD constants
    const float_batch v_m00(m00), v_m01(m01), v_m02(m02), v_m03(m03);
    const float_batch v_m10(m10), v_m11(m11), v_m12(m12), v_m13(m13);
    const float_batch v_m20(m20), v_m21(m21), v_m22(m22), v_m23(m23);

    // Process in SIMD lanes
    size_t i = 0;
    for (; i + simd_size <= static_cast<size_t>(num_atoms); i += simd_size) {
      // Load coordinates and mask
      const auto v_x      = float_batch::load_unaligned(x + i);
      const auto v_y      = float_batch::load_unaligned(y + i);
      const auto v_z      = float_batch::load_unaligned(z + i);
      const auto int_mask = int_batch::load_unaligned(frag_mask + i);
      const auto mask     = xsimd::to_float(int_mask) != 0;

      // Apply rotation only to masked elements
      const auto v_tx =
          xsimd::select(mask,
                        xsimd::fma(v_x, v_m00, xsimd::fma(v_y, v_m01, xsimd::fma(v_z, v_m02, v_m03))),
                        v_x);

      const auto v_ty =
          xsimd::select(mask,
                        xsimd::fma(v_x, v_m10, xsimd::fma(v_y, v_m11, xsimd::fma(v_z, v_m12, v_m13))),
                        v_y);

      const auto v_tz =
          xsimd::select(mask,
                        xsimd::fma(v_x, v_m20, xsimd::fma(v_y, v_m21, xsimd::fma(v_z, v_m22, v_m23))),
                        v_z);

      // Store results
      v_tx.store_unaligned(x + i);
      v_ty.store_unaligned(y + i);
      v_tz.store_unaligned(z + i);
    }

    // Process remaining elements
    for (; i < static_cast<size_t>(num_atoms); ++i) {
      if (frag_mask[i]) {
        const auto tx = x[i], ty = y[i], tz = z[i];
        x[i] = m00 * tx + m01 * ty + m02 * tz + m03;
        y[i] = m10 * tx + m11 * ty + m12 * tz + m13;
        z[i] = m20 * tx + m21 * ty + m22 * tz + m23;
      }
    }
  }
} // namespace mudock
