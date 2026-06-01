#pragma once

#include <mudock/type_alias.hpp>
#include <span>
#include <vector>

namespace mudock {

struct vinardo_pose_score {
  fp_type c_inter;
  fp_type c_intra;
};

fp_type vinardo_affinity(fp_type inter_score, unsigned num_tors);

std::vector<fp_type> vinardo_affinities_from_poses(std::span<const vinardo_pose_score> poses,
                                                   unsigned num_tors);

} // namespace mudock
