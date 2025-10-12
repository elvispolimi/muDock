#pragma once

#include <mudock/type_alias.hpp>

namespace mudock {
  /* ______________________________________________________________________________ */
  /* Nonbonded pair parameters */
  typedef struct nonbond_param {
    int a1; // ATM1
    int a2; // ATM2
    // TODO check this seems not relevant for our case
    int nonbond_type; // NBTYPE  0 = not 1_4     4 = is 1_4

    nonbond_param(): a1(0), a2(0) {}
  } non_bond_parameter;
} // namespace mudock
