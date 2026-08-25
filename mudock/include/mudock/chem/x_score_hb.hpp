#pragma once

#include <string_view>

namespace mudock {

  // HB Type
  //
  //   N  -> none
  //   H  -> hydrophobic
  //   P  -> polar
  //   D  -> donor
  //   A  -> acceptor
  //   DA -> donor/acceptor
  //   DH -> polar hydrogen
  //   M  -> metal

  enum class x_score_hb : int { N = 0, H = 1, P = 2, D = 3, A = 4, DA = 5, DH = 6, M = 7 };


} // namespace mudock