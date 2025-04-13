#pragma once

#include <mudock/grid.hpp>
#include <mudock/grid/grid_map.hpp>
#include <mudock/molecule.hpp>
#include <mudock/molecule/fragments.hpp>
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
                  std::vector<int>& non_bond_list_a1,
                  std::vector<int>& non_bond_list_a2,
                  const int num_atoms,
                  const fragments<static_containers>& ligand_fragments);

  void non_bond_list(const static_molecule& ligand,
                     std::vector<int>& non_bond_list_a1,
                     std::vector<int>& non_bond_list_a2);

  void precompute_lennard_jones(const size_t non_bond_size,
                                std::vector<fp_type>& cA_v,
                                std::vector<fp_type>& cB_v,
                                std::vector<int>& xB_v,
                                const static_molecule& ligand,
                                const std::vector<int>& non_bond_list_a1,
                                const std::vector<int>& non_bond_list_a2);

} // namespace mudock
