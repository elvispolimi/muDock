#include <mudock/compute/vinardo_affinity.hpp>

namespace mudock {

fp_type vinardo_affinity(const fp_type inter_score, const unsigned num_tors) {
  constexpr fp_type w_rot = fp_type{0.02};
  return inter_score / (fp_type{1} + w_rot * static_cast<fp_type>(num_tors));
}

std::vector<fp_type> vinardo_affinities_from_poses(std::span<const vinardo_pose_score> poses,
                                                   const unsigned num_tors) {
  std::vector<fp_type> affinities;
  affinities.reserve(poses.size());
  if (poses.empty()) {
    return affinities;
  }

  const auto best_intra_score = poses.front().c_intra;
  affinities.push_back(vinardo_affinity(poses.front().c_inter, num_tors));

  for (std::size_t i = 1; i < poses.size(); ++i) {
    const auto conformation_score = poses[i].c_inter + poses[i].c_intra;
    affinities.push_back(vinardo_affinity(conformation_score - best_intra_score, num_tors));
  }
  return affinities;
}

} // namespace mudock
