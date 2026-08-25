#pragma once

namespace mudock {

  // Per-atom validity used by the X-Score.
  //
  //   0 -> invalid : atom failed typing/parametrization, or is a non-polar hydrogen.
  //   1 -> valid   : atom successfully read and typed.
  //   2 -> pocket  : protein atom belonging to the binding pocket (Define_Pocket not currently implemented).

  enum class x_score_validity : int { invalid = 0, valid = 1, pocket = 2 };

} // namespace mudock
