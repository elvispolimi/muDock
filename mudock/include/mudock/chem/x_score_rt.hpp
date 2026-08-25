#pragma once

#include <mudock/chem/x_score_ligand.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {

  // Rotor (RT) term of the X-Score scoring function
  //
  // The term is a property of the ligand, and it is computed once per ligand on the host (not pose dependent).

  [[nodiscard]] fp_type compute_x_score_rt(const x_score_ligand& xs_lig);

} // namespace mudock
