
#include <cstdint>
#include <mudock/chem/autodock_ligand.hpp>
#include <mudock/chem/grid_const.hpp>
#include <mudock/log.hpp>
#include <mudock/molecule/containers.hpp>
#include <mudock/molecule/fragments.hpp>
#include <vector>

namespace mudock {
  // nonbonds.cc for nbmatrix required by weed_bonds
  void autodock_ligand::nonbonds(md_container<std::vector<uint_fast8_t>, 2>& nbmatrix,
                                 const std::span<const bond> ligand_bond,
                                 const int num_atoms) {
    //
    // in "nbmatrix", the values 1 (and 4) mean this pair of atoms will be included in the internal, non-bonded list
    //                           0                                         ignored

    // set all nonbonds in nbmatrix to 1, except "1-1 interactions" (self interaction)
    for (int i = 0; i < num_atoms; i++) {
      for (int j = 0; j < num_atoms; j++) { nbmatrix.get(i, j) = 1; } // j
      nbmatrix.get(i, i) = 0;                                         /* 2005-01-10 RH & GMM */
    }

    for (auto& bond: ligand_bond) {
      // Ignore 1-2 Interactions
      nbmatrix.get(bond.source, bond.dest) = 0;
      nbmatrix.get(bond.dest, bond.source) = 0;
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
        nbmatrix.get(outer_2, outer_1) = 0;
        nbmatrix.get(outer_1, outer_2) = 0;

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
          nbmatrix.get(outer_4, outer_3) = 0;
          nbmatrix.get(outer_3, outer_4) = 0;
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
  void autodock_ligand::weed_bonds(md_container<std::vector<uint_fast8_t>, 2>& nbmatrix,
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
          nbmatrix.get(j, i) = 0;
          nbmatrix.get(i, j) = 0;
        }
      } // i
    } // j
    /* 
    \   Weed out bonds across torsions,
    \______________________________________________________________
    */
    // Loop over all "ntor" torsions, "i"
    for (int i = 0; i < ligand_fragments.get_num_rotatable_bonds(); ++i) {
      const auto [atom_id1, atom_id2] = ligand_fragments.get_rotatable_atoms(i);
      // TODO check why not viceversa? weedbonds.cc:110
      nbmatrix.get(atom_id2, atom_id1) = 0;
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
          nbmatrix.get(atom_id4, atom_id2) = 0;
          nbmatrix.get(atom_id2, atom_id4) = 0;
        }
        if (ligand_rigid_pieces[atom_id1] == ligand_rigid_pieces[atom_id4]) {
          nbmatrix.get(atom_id3, atom_id2) = 0;
          nbmatrix.get(atom_id2, atom_id3) = 0;
        }
        if (ligand_rigid_pieces[atom_id2] == ligand_rigid_pieces[atom_id3]) {
          nbmatrix.get(atom_id4, atom_id1) = 0;
          nbmatrix.get(atom_id1, atom_id4) = 0;
        }
        if (ligand_rigid_pieces[atom_id2] == ligand_rigid_pieces[atom_id4]) {
          nbmatrix.get(atom_id3, atom_id1) = 0;
          nbmatrix.get(atom_id1, atom_id3) = 0;
        }
      }
      for (int k = 0; k < num_atoms; ++k) {
        if (ligand_rigid_pieces[atom_id1] == ligand_rigid_pieces[k]) {
          nbmatrix.get(k, atom_id2) = 0;
          nbmatrix.get(atom_id2, k) = 0;
        }
        if (ligand_rigid_pieces[atom_id2] == ligand_rigid_pieces[k]) {
          nbmatrix.get(k, atom_id1) = 0;
          nbmatrix.get(atom_id1, k) = 0;
        }
      } // k
    }

    // intramolecular non-bonds for ligand
    // TODO check what true_ligand_atoms is
    for (int i = 0; i < num_atoms; ++i) {
      for (int j = i + 1; j < num_atoms; ++j) {
        if ((nbmatrix.get(i, j) == 1 && nbmatrix.get(j, i) == 1)) {
          non_bond_list_a1.push_back(i);
          non_bond_list_a2.push_back(j);
        } else if ((nbmatrix.get(i, j) != 0 && nbmatrix.get(j, i) == 0) ||
                   (nbmatrix.get(i, j) == 0 && nbmatrix.get(j, i) != 0)) {
          std::ostringstream oss;
          // Build the formatted string
          oss << "BUG: ASSYMMETRY detected in Non-Bond Matrix at " << i << "," << j;
          error(oss.str());
        }
      } // j
    } // i
  }

  void autodock_ligand::non_bond_list() {
    const auto num_atoms = this->get_base_molecule().num_atoms();

    auto graph = make_graph(this->get_base_molecule().get_bonds(), num_atoms);
    const auto ligand_fragments =
        std::make_unique<mudock::fragments<mudock::static_containers>>(graph,
                                                                       this->get_base_molecule().get_bonds(),
                                                                       num_atoms);

    md_container<std::vector<uint_fast8_t>, 2> nbmatrix{num_atoms, num_atoms};
    // grid<uint_fast8_t, index2D> nbmatrix{{num_atoms, num_atoms}};
    nonbonds(nbmatrix, this->get_base_molecule().get_bonds(), num_atoms);
    weed_bonds(nbmatrix, num_atoms, *ligand_fragments);
  }

  void autodock_ligand::precompute_lennard_jones() {
    const auto non_bond_size = non_bond_list_a1.size();
    cA_v.resize(non_bond_size);
    cB_v.resize(non_bond_size);
    xB_v.resize(non_bond_size);
    for (size_t index = 0; index < non_bond_size; ++index) {
      const int& a1 = non_bond_list_a1[index];
      const int& a2 = non_bond_list_a2[index];

      const auto& hbond_i    = this->get_base_molecule().num_hbond(a1);
      const auto& hbond_j    = this->get_base_molecule().num_hbond(a2);
      const auto& Rij_hb_i   = Rij_hb(a1);
      const auto& Rij_hb_j   = Rij_hb(a2);
      const auto& Rii_i      = Rii(a1);
      const auto& Rii_j      = Rii(a2);
      const auto& epsij_hb_i = epsij_hb(a1);
      const auto& epsij_hb_j = epsij_hb(a2);
      const auto& epsii_i    = epsii(a1);
      const auto& epsii_j    = epsii(a2);

      // we need to determine the correct xA and xB exponents
      const int xA = xA_default; // for both LJ, 12-6 and HB, 12-10, xA is 12
      int xB       = xB_default; // assume we have LJ, 12-6

      fp_type Rij{(Rii_i + Rii_j) * fp_type{0.5}}, epsij{std::sqrt(epsii_i * epsii_j)};
      if ((hbond_i == 1 || hbond_i == 2) && hbond_j > 2) {
        // i is a donor and j is an acceptor.
        // i is a hydrogen, j is a heteroatom
        Rij   = Rij_hb_j;
        epsij = epsij_hb_j;
        xB    = 10;
      } else if ((hbond_i > 2) && (hbond_j == 1 || hbond_j == 2)) {
        // i is an acceptor and j is a donor.
        // i is a heteroatom, j is a hydrogen
        Rij   = Rij_hb_i;
        epsij = epsij_hb_i;
        xB    = 10;
      }
      fp_type cA{0};
      fp_type cB{0};
      if (xA != xB) {
        const fp_type tmp = epsij / static_cast<fp_type>(xA - xB);
        cA                = tmp * std::pow(Rij, static_cast<fp_type>(xA)) * static_cast<fp_type>(xB);
        cB                = tmp * std::pow(Rij, static_cast<fp_type>(xB)) * static_cast<fp_type>(xA);
      }
      cA_v[index] = cA;
      cB_v[index] = cB;
      xB_v[index] = xB;
    }
  }
} // namespace mudock
