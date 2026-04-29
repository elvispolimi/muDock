#pragma once

#include <array>
#include <cmath>
#include <mudock/molecule.hpp>
#include <stdexcept>
#include <vector>

namespace mudock {

  namespace detail {
    inline std::array<fp_type, 4> dominant_quaternion(const std::array<std::array<fp_type, 4>, 4>& k_matrix) {
      std::array<fp_type, 4> q{fp_type{1}, fp_type{0}, fp_type{0}, fp_type{0}};

      for (int iteration = 0; iteration < 64; ++iteration) {
        std::array<fp_type, 4> next{fp_type{0}, fp_type{0}, fp_type{0}, fp_type{0}};
        for (int row = 0; row < 4; ++row) {
          for (int col = 0; col < 4; ++col) {
            next[row] += k_matrix[row][col] * q[col];
          }
        }

        const fp_type norm =
            std::sqrt(next[0] * next[0] + next[1] * next[1] + next[2] * next[2] + next[3] * next[3]);
        if (norm <= fp_type{1e-12}) {
          return {fp_type{1}, fp_type{0}, fp_type{0}, fp_type{0}};
        }

        for (auto& value: next) {
          value /= norm;
        }
        q = next;
      }

      return q;
    }
  } // namespace detail

  template<class lhs_molecule, class rhs_molecule>
    requires is_molecule<lhs_molecule> && is_molecule<rhs_molecule>
  fp_type aligned_rmsd(const lhs_molecule& lhs, const rhs_molecule& rhs) {
    if (lhs.num_atoms() != rhs.num_atoms()) {
      throw std::invalid_argument("RMSD requires ligands with the same number of atoms");
    }
    if (lhs.num_atoms() <= 0) {
      throw std::invalid_argument("RMSD requires at least one atom");
    }

    const int num_atoms = lhs.num_atoms();
    std::vector<std::array<fp_type, 3>> lhs_centered(static_cast<std::size_t>(num_atoms));
    std::vector<std::array<fp_type, 3>> rhs_centered(static_cast<std::size_t>(num_atoms));

    std::array<fp_type, 3> lhs_centroid{fp_type{0}, fp_type{0}, fp_type{0}};
    std::array<fp_type, 3> rhs_centroid{fp_type{0}, fp_type{0}, fp_type{0}};
    for (int atom_index = 0; atom_index < num_atoms; ++atom_index) {
      lhs_centroid[0] += lhs.x(atom_index);
      lhs_centroid[1] += lhs.y(atom_index);
      lhs_centroid[2] += lhs.z(atom_index);
      rhs_centroid[0] += rhs.x(atom_index);
      rhs_centroid[1] += rhs.y(atom_index);
      rhs_centroid[2] += rhs.z(atom_index);
    }
    for (int axis = 0; axis < 3; ++axis) {
      lhs_centroid[axis] /= static_cast<fp_type>(num_atoms);
      rhs_centroid[axis] /= static_cast<fp_type>(num_atoms);
    }

    fp_type sxx{0}, sxy{0}, sxz{0};
    fp_type syx{0}, syy{0}, syz{0};
    fp_type szx{0}, szy{0}, szz{0};
    for (int atom_index = 0; atom_index < num_atoms; ++atom_index) {
      const std::array<fp_type, 3> lhs_point{lhs.x(atom_index) - lhs_centroid[0],
                                             lhs.y(atom_index) - lhs_centroid[1],
                                             lhs.z(atom_index) - lhs_centroid[2]};
      const std::array<fp_type, 3> rhs_point{rhs.x(atom_index) - rhs_centroid[0],
                                             rhs.y(atom_index) - rhs_centroid[1],
                                             rhs.z(atom_index) - rhs_centroid[2]};

      lhs_centered[static_cast<std::size_t>(atom_index)] = lhs_point;
      rhs_centered[static_cast<std::size_t>(atom_index)] = rhs_point;

      sxx += lhs_point[0] * rhs_point[0];
      sxy += lhs_point[0] * rhs_point[1];
      sxz += lhs_point[0] * rhs_point[2];
      syx += lhs_point[1] * rhs_point[0];
      syy += lhs_point[1] * rhs_point[1];
      syz += lhs_point[1] * rhs_point[2];
      szx += lhs_point[2] * rhs_point[0];
      szy += lhs_point[2] * rhs_point[1];
      szz += lhs_point[2] * rhs_point[2];
    }

    const std::array<std::array<fp_type, 4>, 4> k_matrix{
        std::array<fp_type, 4>{sxx + syy + szz, syz - szy, szx - sxz, sxy - syx},
        std::array<fp_type, 4>{syz - szy, sxx - syy - szz, sxy + syx, szx + sxz},
        std::array<fp_type, 4>{szx - sxz, sxy + syx, -sxx + syy - szz, syz + szy},
        std::array<fp_type, 4>{sxy - syx, szx + sxz, syz + szy, -sxx - syy + szz}};

    const auto quaternion = detail::dominant_quaternion(k_matrix);
    const fp_type w       = quaternion[0];
    const fp_type x       = quaternion[1];
    const fp_type y       = quaternion[2];
    const fp_type z       = quaternion[3];

    const fp_type r00 = w * w + x * x - y * y - z * z;
    const fp_type r01 = fp_type{2} * (x * y - w * z);
    const fp_type r02 = fp_type{2} * (x * z + w * y);
    const fp_type r10 = fp_type{2} * (x * y + w * z);
    const fp_type r11 = w * w - x * x + y * y - z * z;
    const fp_type r12 = fp_type{2} * (y * z - w * x);
    const fp_type r20 = fp_type{2} * (x * z - w * y);
    const fp_type r21 = fp_type{2} * (y * z + w * x);
    const fp_type r22 = w * w - x * x - y * y + z * z;

    fp_type sum_sq{0};
    for (int atom_index = 0; atom_index < num_atoms; ++atom_index) {
      const auto& lhs_point = lhs_centered[static_cast<std::size_t>(atom_index)];
      const auto& rhs_point = rhs_centered[static_cast<std::size_t>(atom_index)];

      const fp_type rot_x = r00 * lhs_point[0] + r01 * lhs_point[1] + r02 * lhs_point[2];
      const fp_type rot_y = r10 * lhs_point[0] + r11 * lhs_point[1] + r12 * lhs_point[2];
      const fp_type rot_z = r20 * lhs_point[0] + r21 * lhs_point[1] + r22 * lhs_point[2];

      const fp_type dx = rot_x - rhs_point[0];
      const fp_type dy = rot_y - rhs_point[1];
      const fp_type dz = rot_z - rhs_point[2];
      sum_sq += dx * dx + dy * dy + dz * dz;
    }

    return std::sqrt(sum_sq / static_cast<fp_type>(num_atoms));
  }

} // namespace mudock
