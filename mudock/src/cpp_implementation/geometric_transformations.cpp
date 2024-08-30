#include <cassert>
#include <mudock/cpp_implementation/center_of_mass.hpp>
#include <mudock/cpp_implementation/geometric_transformations.hpp>
#include <mudock/grid.hpp>

namespace mudock {

  void translate_molecule(std::span<fp_type> x,
                          std::span<fp_type> y,
                          std::span<fp_type> z,
                          const fp_type offset_x,
                          const fp_type offset_y,
                          const fp_type offset_z) {
    const auto num_atoms = x.size();
    assert(y.size() == num_atoms);
    assert(z.size() == num_atoms);
    for (std::size_t i = 0; i < num_atoms; ++i) {
      x[i] += offset_x;
      y[i] += offset_y;
      z[i] += offset_z;
    }
  }

  void rotate_molecule(std::span<fp_type> x,
                       std::span<fp_type> y,
                       std::span<fp_type> z,
                       const fp_type angle_x,
                       const fp_type angle_y,
                       const fp_type angle_z) {
    // get the molecule number of atoms
    const auto num_atoms = x.size();
    assert(y.size() == num_atoms);
    assert(z.size() == num_atoms);

    // compute the molecule center of mass
    const auto c = compute_center_of_mass(x, y, z);

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
    for (std::size_t i = 0; i < num_atoms; ++i) {
      const auto translated_x = x[i] - c.x, translated_y = y[i] - c.y, translated_z = z[i] - c.z;
      x[i] = translated_x * m00 + translated_y * m01 + translated_z * m02 + c.x;
      y[i] = translated_x * m10 + translated_y * m11 + translated_z * m12 + c.y;
      z[i] = translated_x * m20 + translated_y * m21 + translated_z * m22 + c.z;
    }
  }

  void rotate_fragment(std::span<fp_type> x,
                       std::span<fp_type> y,
                       std::span<fp_type> z,
                       std::span<const typename fragments<static_containers>::value_type> bitmask,
                       const std::size_t start_index,
                       const std::size_t stop_index,
                       const fp_type angle) {
    // get the molecule number of atoms
    const auto num_atoms = x.size();
    assert(y.size() == num_atoms);
    assert(z.size() == num_atoms);
    assert(start_index < num_atoms);
    assert(stop_index < num_atoms);

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
    for (std::size_t i = 0; i < num_atoms; ++i) {
      if (bitmask[i] == 1) {
        const auto prev_x = x[i], prev_y = y[i], prev_z = z[i];
        x[i] = prev_x * m00 + prev_y * m01 + prev_z * m02 + m03;
        y[i] = prev_x * m10 + prev_y * m11 + prev_z * m12 + m13;
        z[i] = prev_x * m20 + prev_y * m21 + prev_z * m22 + m23;
      }
    }
  }

} // namespace mudock
