#pragma once

#include <mudock/cpp_implementation/center_of_mass.hpp>
#include <mudock/grid/point3D.hpp>
#include <mudock/molecule.hpp>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>

namespace mudock {

  enum class ligand_placement_mode { center_bbox, preserve, point, probe };

  inline ligand_placement_mode parse_ligand_placement_mode(const std::string_view value) {
    if (value == "center_bbox") return ligand_placement_mode::center_bbox;
    if (value == "preserve") return ligand_placement_mode::preserve;
    if (value == "point") return ligand_placement_mode::point;
    if (value == "probe") return ligand_placement_mode::probe;
    throw std::runtime_error("Invalid placement mode '" + std::string(value) +
                             "'. Expected center_bbox|preserve|point|probe.");
  }

  inline std::string_view to_string(const ligand_placement_mode mode) {
    switch (mode) {
      case ligand_placement_mode::center_bbox: return "center_bbox";
      case ligand_placement_mode::preserve: return "preserve";
      case ligand_placement_mode::point: return "point";
      case ligand_placement_mode::probe: return "probe";
    }
    throw std::runtime_error("Invalid ligand placement mode");
  }

  struct ligand_placement {
    ligand_placement_mode mode = ligand_placement_mode::center_bbox;
    std::optional<point3D> target;
  };

  inline point3D apply_ligand_placement(static_molecule& ligand,
                                        const ligand_placement& placement,
                                        const point3D protein_center) {
    if (placement.mode == ligand_placement_mode::preserve) return point3D{};

    const auto target = placement.mode == ligand_placement_mode::center_bbox ?
                            std::optional<point3D>{protein_center} : placement.target;
    if (!target.has_value()) throw std::runtime_error("Ligand placement target is missing");

    const auto centroid = compute_centroid(ligand.x(), ligand.y(), ligand.z(), ligand.num_atoms());
    const auto offset   = target.value() - centroid;
    for (int atom = 0; atom < ligand.num_atoms(); ++atom) {
      ligand.x(atom) += offset.x();
      ligand.y(atom) += offset.y();
      ligand.z(atom) += offset.z();
    }
    return offset;
  }

} // namespace mudock
