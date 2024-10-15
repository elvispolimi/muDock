#pragma once

#include <mudock/grid.hpp>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>
#include <span>
#include <vector>

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

  void nonbonds(grid<uint_fast8_t, index2D>& nbmatrix,
                const std::span<const bond> ligand_bond,
                const int num_atoms);

  // weedbonds.cc for nonbondlist only for the first group
  /*___________________________________________________________________________
  |    ENDBRANCH---TORS---BRANCH---R O O T---BRANCH---ENDBRANCH                |
  |                                  /              \                          |
  |                                BRANCH            BRANCH--TORS---ENDBRANCH  |
  |                                /                  \                        |
  |                               ENDBRANCH            ENDBRANCH               |
  |____________________________________________________________________________|
  |  Eliminate all rigidly bonded atoms:                                       |
  |                                     * those atoms which are at the ROOT;   |
  |                                     * atoms between TORS and BRANCH;       |
  |                                     * atoms between BRANCH and ENDBRANCH.  |
  |  This is necessary for internal energy calculations.                       |
  |____________________________________________________________________________|
  | Weed out bonds in rigid pieces,                                            |
  |____________________________________________________________________________|
  */
  void weed_bonds(grid<uint_fast8_t, index2D>& nbmatrix,
                  std::vector<non_bond_parameter>& non_bond_list,
                  const int num_atoms,
                  const fragments<static_containers>& ligand_fragments);
} // namespace mudock
