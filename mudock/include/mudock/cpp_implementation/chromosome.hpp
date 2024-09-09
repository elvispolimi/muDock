#pragma once

#include <array>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {

  /**
   * A Chromosome represents a geometric transformation that we apply to alter
   * the 3D displacement and shape of a given molecule. In this version we use
   * the following genes:
   *  0 -> Offset (in Angstrom) of a translation along the x axys
   *  1 -> Offset (in Angstrom) of a translation along the y axys
   *  2 -> Offset (in Angstrom) of a translation along the z axys
   *  3 -> Angle (in degree) that we rotate the whole molecule along the x axys
   *  4 -> Angle (in degree) that we rotate the whole molecule along the y axys
   *  5 -> Angle (in degree) that we rotate the whole molecule along the z axys
   *  6 -> Angle (in degree) that we rotate a fragment along the first rotable bond (if any)
   *  7 -> Angle (in degree) that we rotate a fragment along the second rotable bond (if any)
   *  8 -> Angle (in degree) that we rotate a fragment along the third rotatable bond (if any)
   *  9 -> ...
   *
   * The actual number of genes depends on the given molecule
  */
  using chromosome = std::array<fp_type, 6 + max_static_bonds()>;

  struct individual {
    chromosome genes;
    fp_type score;
  };
} // namespace mudock
