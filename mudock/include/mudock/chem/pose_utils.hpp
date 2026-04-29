#pragma once

#include <mudock/chem/geom_ligand.hpp>
#include <mudock/cpp_implementation/mutate.hpp>
#include <mudock/molecule.hpp>
#include <string>
#include <string_view>

namespace mudock {

  inline void apply_pose(static_molecule& ligand, const chromosome& pose) {
    geom_ligand geom_lig{ligand};
    apply<cpu_vectorization::AUTO>(ligand.x(),
                                   ligand.y(),
                                   ligand.z(),
                                   pose,
                                   ligand.num_atoms(),
                                   ligand.num_rotamers(),
                                   geom_lig.fragments_masks(),
                                   geom_lig.fragmets_starts(),
                                   geom_lig.fragments_stops());
  }

  inline std::string pose_name(std::string_view base_name, const std::size_t pose_index) {
    return std::string(base_name) + "#" + std::to_string(pose_index);
  }

} // namespace mudock
