#include <cmath>
#include <mudock/compute/ligand_placement.hpp>

namespace {
  mudock::static_molecule make_ligand() {
    mudock::static_molecule ligand;
    ligand.resize(3, 0);
    ligand.x(0) = 0;
    ligand.y(0) = 0;
    ligand.z(0) = 0;
    ligand.x(1) = 2;
    ligand.y(1) = 0;
    ligand.z(1) = 0;
    ligand.x(2) = 0;
    ligand.y(2) = 2;
    ligand.z(2) = 0;
    return ligand;
  }

  bool close(const mudock::fp_type lhs, const mudock::fp_type rhs) {
    return std::fabs(lhs - rhs) < static_cast<mudock::fp_type>(0.00001);
  }
} // namespace

int main() {
  const mudock::point3D target{mudock::fp_type{10}, mudock::fp_type{20}, mudock::fp_type{30}};

  auto centered = make_ligand();
  mudock::apply_ligand_placement(
      centered, mudock::ligand_placement{mudock::ligand_placement_mode::center_bbox, std::nullopt}, target);
  const auto centered_centroid =
      mudock::compute_centroid(centered.x(), centered.y(), centered.z(), centered.num_atoms());
  if (!close(centered_centroid.x(), target.x()) || !close(centered_centroid.y(), target.y()) ||
      !close(centered_centroid.z(), target.z()) ||
      !close(centered.x(0) - centered.x(1), mudock::fp_type{-2}))
    return 1;

  auto preserved = make_ligand();
  mudock::apply_ligand_placement(
      preserved, mudock::ligand_placement{mudock::ligand_placement_mode::preserve, std::nullopt}, target);
  if (preserved.x(1) != mudock::fp_type{2} || preserved.y(2) != mudock::fp_type{2}) return 1;

  auto pointed = make_ligand();
  mudock::ligand_placement point_placement{mudock::ligand_placement_mode::point, target};
  mudock::apply_ligand_placement(pointed, point_placement, mudock::point3D{});
  const auto pointed_centroid =
      mudock::compute_centroid(pointed.x(), pointed.y(), pointed.z(), pointed.num_atoms());
  if (!close(pointed_centroid.x(), target.x()) || !close(pointed_centroid.y(), target.y()) ||
      !close(pointed_centroid.z(), target.z()))
    return 1;

  auto probed = make_ligand();
  mudock::ligand_placement probe_placement{mudock::ligand_placement_mode::probe, target};
  mudock::apply_ligand_placement(probed, probe_placement, mudock::point3D{});
  const auto probed_centroid =
      mudock::compute_centroid(probed.x(), probed.y(), probed.z(), probed.num_atoms());
  if (!close(probed_centroid.x(), target.x()) || !close(probed_centroid.y(), target.y()) ||
      !close(probed_centroid.z(), target.z()))
    return 1;
  return 0;
}
