#include <mudock/cpp_implementation/weed_bonds.hpp>

namespace mudock {
  // nonbonds.cc for nbmatrix required by weed_bonds
  void nonbonds(grid<uint_fast8_t, index2D>& nbmatrix,
                const std::span<const bond> ligand_bond,
                const int num_atoms) {
    //
    // in "nbmatrix", the values 1 (and 4) mean this pair of atoms will be included in the internal, non-bonded list
    //                           0                                         ignored

    // set all nonbonds in nbmatrix to 1, except "1-1 interactions" (self interaction)
    for (int i = 0; i < num_atoms; i++) {
      for (int j = 0; j < num_atoms; j++) { nbmatrix.at(i, j) = 1; } // j
      nbmatrix.at(i, i) = 0;                                         /* 2005-01-10 RH & GMM */
    }

    for (auto& bond: ligand_bond) {
      // Ignore 1-2 Interactions
      nbmatrix.at(bond.source, bond.dest) = 0;
      nbmatrix.at(bond.dest, bond.source) = 0;
    }

    for (auto& bond_1: ligand_bond)
      for (auto& bond_2: ligand_bond) { // loop over each atom "k" bonded to the current atom "j"
        int outer_1{0};
        int outer_2{0};
        if (bond_1.dest == bond_2.source) {
          outer_1 = bond_1.source;
          outer_2 = bond_2.dest;
        } else if (bond_1.dest == bond_2.dest) {
          outer_1 = bond_1.source;
          outer_2 = bond_2.source;
        } else if (bond_1.source == bond_2.source) {
          outer_1 = bond_1.dest;
          outer_2 = bond_2.dest;
        } else if (bond_1.source == bond_2.dest) {
          outer_1 = bond_1.dest;
          outer_2 = bond_2.source;
        } else
          continue;

        // Ignore "1-3 Interactions"
        nbmatrix.at(outer_2, outer_1) = 0;
        nbmatrix.at(outer_1, outer_2) = 0;

        for (auto& bond_3: ligand_bond) {
          int outer_3{0};
          int outer_4{0};
          if (outer_2 == bond_3.source) {
            outer_3 = bond_3.dest;
            outer_4 = outer_1;
          } else if (outer_2 == bond_3.dest) {
            outer_3 = bond_3.source;
            outer_4 = outer_1;
          } else if (outer_1 == bond_3.source) {
            outer_3 = bond_1.dest;
            outer_4 = outer_2;
          } else if (outer_1 == bond_3.dest) {
            outer_3 = bond_1.source;
            outer_4 = outer_2;
          } else
            continue;

          nbmatrix.at(outer_4, outer_3) = 0;
          nbmatrix.at(outer_3, outer_4) = 0;
        }
      }
  }

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
                  const fragments<static_containers>& ligand_fragments) {
    const auto& ligand_rigid_pieces = ligand_fragments.get_rigid_pieces();

    for (int j = 0; j < num_atoms; ++j) {
      for (int i = 0; i < num_atoms; ++i) {
        // Is atom "i" in the same rigid piece as atom "j"?
        if (ligand_rigid_pieces[i] == ligand_rigid_pieces[j]) {
          // Set the entry for atoms "i" and "j" in the
          //    nonbond matrix to 0 ;
          //
          // Later on, we will not calculate the interaction energy
          //    between atoms "i" and "j"
          nbmatrix.at(j, i) = 0;
          nbmatrix.at(i, j) = 0;
        }
      } // i
    }   // j
    /* 
    \   Weed out bonds across torsions,
    \______________________________________________________________
    */
    // Loop over all "ntor" torsions, "i"
    for (int i = 0; i < ligand_fragments.get_num_rotatable_bonds(); ++i) {
      const auto [atom_id1, atom_id2] = ligand_fragments.get_rotatable_atoms(i);
      // TODO check why not viceversa? weedbonds.cc:110
      nbmatrix.at(atom_id2, atom_id1) = 0;
    } // i

    /* 
    \  Weed out bonds from atoms directly connected to rigid pieces,
    \_ we think these are 1-3 interactions mp+rh, 10-2008______________________
    */
    for (int i = 0; i < ligand_fragments.get_num_rotatable_bonds(); ++i) {
      const auto [atom_id1, atom_id2] = ligand_fragments.get_rotatable_atoms(i);
      for (int j = 0; j < ligand_fragments.get_num_rotatable_bonds(); ++j) {
        const auto [atom_id3, atom_id4] = ligand_fragments.get_rotatable_atoms(j);
        if (ligand_rigid_pieces[atom_id1] == ligand_rigid_pieces[atom_id3]) {
          nbmatrix.at(atom_id4, atom_id2) = 0;
          nbmatrix.at(atom_id2, atom_id4) = 0;
        }
        if (ligand_rigid_pieces[atom_id1] == ligand_rigid_pieces[atom_id4]) {
          nbmatrix.at(atom_id3, atom_id2) = 0;
          nbmatrix.at(atom_id2, atom_id3) = 0;
        }
        if (ligand_rigid_pieces[atom_id2] == ligand_rigid_pieces[atom_id3]) {
          nbmatrix.at(atom_id4, atom_id1) = 0;
          nbmatrix.at(atom_id1, atom_id4) = 0;
        }
        if (ligand_rigid_pieces[atom_id2] == ligand_rigid_pieces[atom_id4]) {
          nbmatrix.at(atom_id3, atom_id1) = 0;
          nbmatrix.at(atom_id1, atom_id3) = 0;
        }
      }
      for (int k = 0; k < num_atoms; ++k) {
        if (ligand_rigid_pieces[atom_id1] == ligand_rigid_pieces[k]) {
          nbmatrix.at(k, atom_id2) = 0;
          nbmatrix.at(atom_id2, k) = 0;
        }
        if (ligand_rigid_pieces[atom_id2] == ligand_rigid_pieces[k]) {
          nbmatrix.at(k, atom_id1) = 0;
          nbmatrix.at(atom_id1, k) = 0;
        }
      } // k
    }

    // intramolecular non-bonds for ligand
    // TODO check what true_ligand_atoms is
    for (int i = 0; i < num_atoms; ++i) {
      for (int j = i + 1; j < num_atoms; ++j) {
        if ((nbmatrix.at(i, j) == 1 && nbmatrix.at(j, i) == 1)) {
          non_bond_list.emplace_back();
          auto& nbl        = non_bond_list.back();
          nbl.a1           = i;
          nbl.a2           = j;
          nbl.nonbond_type = nbmatrix.at(i, j);
        } else if ((nbmatrix.at(i, j) != 0 && nbmatrix.at(j, i) == 0) ||
                   (nbmatrix.at(i, j) == 0 && nbmatrix.at(j, i) != 0)) {
          std::ostringstream oss;
          // Build the formatted string
          oss << "BUG: ASSYMMETRY detected in Non-Bond Matrix at " << i << "," << j;
          error(oss.str());
        }
      } // j
    }   // i
  }
} // namespace mudock