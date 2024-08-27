
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <mudock/chem.hpp>
#include <mudock/grid.hpp>
#include <mudock/kernels/calc_energy_cpp.hpp>
#include <mudock/molecule.hpp>
#include <stdexcept>

namespace mudock {
  static constexpr fp_type RMIN_ELEC{0.5};
  static constexpr fp_type ELECSCALE{332.06363};
  static constexpr fp_type qsolpar{0.01097};
  // Non bond cutoff
  static constexpr fp_type nbc2{64 * 64};
  static constexpr fp_type ENERGYPENALTY{500};

  /* ______________________________________________________________________________ */
  /* Nonbonded pair parameters */
  typedef struct nonbond_param {
    int a1;           // ATM1
    int a2;           // ATM2
    int nonbond_type; // NBTYPE  0 = not 1_4     4 = is 1_4

    nonbond_param(): a1(0), a2(0) {}
  } non_bond_parameter;

  fp_type trilinear_interpolation(const grid_map& map, const point3D& coord) {
    const point3D p1 = map.get_index_from_coordinates(coord);
    const point3D p2 = map.get_index_from_coordinates(add(coord, point3D{grid_spacing}));

    // TODO autodock use a different method, not much readable
    const point3D xyz_d   = divide(difference(coord, p1), difference(p2, p1));
    const point3D xyz_d_o = difference(static_cast<const point3D>(point3D{1}), xyz_d);

    const fp_type c00 = map.at(p1) * xyz_d_o.x + map.at({p2.x, p1.y, p1.z}) * xyz_d.x;
    const fp_type c10 = map.at({p1.x, p2.y, p1.z}) * xyz_d_o.x + map.at({p2.x, p2.y, p1.z}) * xyz_d.x;
    const fp_type c01 = map.at({p1.x, p1.y, p2.z}) * xyz_d_o.x + map.at({p2.x, p1.y, p2.z}) * xyz_d.x;
    const fp_type c11 = map.at({p1.x, p2.y, p2.z}) * xyz_d_o.x + map.at({p1.x, p1.y, p1.z}) * xyz_d.x;

    const fp_type c0 = c00 * xyz_d_o.y + c10 * xyz_d.y;
    const fp_type c1 = c01 * xyz_d_o.y + c11 * xyz_d.y;

    return c0 * xyz_d_o.z + c1 * xyz_d.z;
  }

  // nonbonds.cc for nbmatrix required by weed_bonds
  void nonbonds(grid<uint_fast8_t, index2D>& nbmatrix, const static_molecule& ligand) {
    //
    // in "nbmatrix", the values 1 (and 4) mean this pair of atoms will be included in the internal, non-bonded list
    //                           0                                         ignored
    const size_t num_atoms = ligand.num_atoms();

    // set all nonbonds in nbmatrix to 1, except "1-1 interactions" (self interaction)
    for (size_t i = 0; i < num_atoms; i++) {
      for (size_t j = 0; j < num_atoms; j++) { nbmatrix.at(i, j) = 1; } // j
      nbmatrix.at(i, i) = 0;                                            /* 2005-01-10 RH & GMM */
    }
    for (auto& bond: ligand.get_bonds()) {
      // Ignore 1-2 Interactions
      nbmatrix.at(bond.source, bond.dest) = 0;
      nbmatrix.at(bond.dest, bond.source) = 0;
    }

    for (auto& bond_1: ligand.get_bonds())
      for (auto& bond_2: ligand.get_bonds()) { // loop over each atom "k" bonded to the current atom "j"
        if (bond_1.dest != bond_2.source)
          continue;

        // Ignore "1-3 Interactions"
        nbmatrix.at(bond_2.dest, bond_1.source) = 0;
        nbmatrix.at(bond_1.source, bond_2.dest) = 0;

        for (auto& bond_3: ligand.get_bonds()) {
          if (bond_2.dest != bond_3.source)
            continue;
          nbmatrix.at(bond_1.source, bond_3.dest) = 0;
          nbmatrix.at(bond_3.dest, bond_1.source) = 0;
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
                  const static_molecule& ligand,
                  const fragments<static_containers>& ligand_fragments) {
    const auto& ligand_rigid_pieces = ligand_fragments.get_rigid_pieces();

    for (size_t j = 0; j < ligand.num_atoms(); ++j) {
      for (size_t i = 0; i < ligand.num_atoms(); ++i) {
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
    for (size_t i = 0; i < ligand_fragments.get_num_rotatable_bonds(); ++i) {
      const auto [atom_id1, atom_id2] = ligand_fragments.get_rotatable_atoms(i);
      // TODO check why not viceversa? weedbonds.cc:110
      nbmatrix.at(atom_id2, atom_id1) = 0;
    } // i

    /* 
    \  Weed out bonds from atoms directly connected to rigid pieces,
    \_ we think these are 1-3 interactions mp+rh, 10-2008______________________
    */
    for (size_t i = 0; i < ligand_fragments.get_num_rotatable_bonds(); ++i) {
      const auto [atom_id1, atom_id2] = ligand_fragments.get_rotatable_atoms(i);
      for (size_t j = 0; j < ligand_fragments.get_num_rotatable_bonds(); ++j) {
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
      for (size_t k = 0; k < ligand.num_atoms(); ++k) {
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
    for (size_t i = 0; i < ligand.num_atoms(); ++i) {
      for (size_t j = i + 1; j < ligand.num_atoms(); ++j) {
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

  fp_type calc_energy(const dynamic_molecule& receptor,
                      const static_molecule& ligand,
                      const fragments<static_containers>& ligand_fragments,
                      const grid_atom_mapper& grid_maps,
                      const grid_map& electro_map,
                      const grid_map& desolv_map) {
    fp_type elect_total = 0;
    fp_type emap_total  = 0;
    fp_type dmap_total  = 0;
    for (size_t index = 0; index < ligand.num_atoms(); ++index) {
      const point3D coord{ligand.x(index), ligand.y(index), ligand.z(index)};
      const auto& atom_charge = ligand.charge(index);
      const auto& atom_map    = grid_maps.get_atom_map(receptor.autodock_type(index));

      if (atom_map.outside_grid(coord)) {
        const fp_type dist = distance2(coord, grid_maps.get_center());

        const fp_type epenalty = dist * ENERGYPENALTY;
        elect_total += epenalty;
        emap_total += epenalty;
      } else {
        // Trilinear Interpolation
        elect_total += trilinear_interpolation(electro_map, coord) * atom_charge;
        emap_total += trilinear_interpolation(atom_map, coord);
        dmap_total += trilinear_interpolation(desolv_map, coord) * std::fabs(atom_charge);
      }
    }

    const int n_torsions = ligand_fragments.get_num_rotatable_bonds();
    if (n_torsions > 0) {
      const size_t num_atoms = ligand.num_atoms();
      // TODO @Davide suppose that the receptor does not have Flexible residues eintcal.cc:147
      // TODO
      grid<uint_fast8_t, index2D> nbmatrix{{num_atoms, num_atoms}};
      nonbonds(nbmatrix, ligand);
      std::vector<non_bond_parameter> non_bond_list(ligand.num_bonds());
      weed_bonds(nbmatrix, non_bond_list, ligand, ligand_fragments);

      for (const auto& non_bond: non_bond_list) {
        const int& a1 = non_bond.a1;
        const int& a2 = non_bond.a2;

        fp_type distance = distance2(point3D{ligand.x(a1), ligand.y(a1), ligand.z(a1)},
                                     point3D{ligand.x(a2), ligand.y(a2), ligand.z(a2)});
        distance         = std::sqrt(std::clamp(distance, RMIN_ELEC * RMIN_ELEC, distance));

        //  Calculate  Electrostatic  Energy
        const fp_type r_dielectric = distance * calc_ddd_Mehler_Solmajer(distance);
        const fp_type e_elec       = ligand.charge(non_bond.a1) * ligand.charge(non_bond.a2) * ELECSCALE *
                               autodock_parameters::coeff_estat * r_dielectric;
        elect_total += e_elec;

        const fp_type nb_desolv =
            (ligand.vol(a2) * (ligand.solpar(a1) + qsolpar * std::fabs(ligand.charge(a1))) +
             ligand.vol(a1) * (ligand.solpar(a2) + qsolpar * std::fabs(ligand.charge(a2))));

        const fp_type e_desolv = autodock_parameters::coeff_desolv *
                                 std::exp(fp_type{-0.5} * sigma * sigma * std::sqrt(distance)) * nb_desolv;
        dmap_total += e_desolv;
        fp_type e_vdW_Hb{0};
        if (distance < nbc2) {
          //  Find internal energy parameters, i.e.  epsilon and r-equilibrium values...
          //  Lennard-Jones and Hydrogen Bond Potentials
          // This can be precomputed as in intnbtable.cc
          const auto& hbond_i    = ligand.num_hbond(a1);
          const auto& hbond_j    = ligand.num_hbond(a2);
          const auto& Rij_hb_i   = ligand.Rij_hb(a1);
          const auto& Rij_hb_j   = ligand.Rij_hb(a2);
          const auto& Rii_i      = ligand.Rii(a1);
          const auto& Rii_j      = ligand.Rii(a2);
          const auto& epsij_hb_i = ligand.epsij_hb(a1);
          const auto& epsij_hb_j = ligand.epsij_hb(a2);
          const auto& epsii_i    = ligand.epsii(a1);
          const auto& epsii_j    = ligand.epsii(a2);

          // we need to determine the correct xA and xB exponents
          size_t xA = 12; // for both LJ, 12-6 and HB, 12-10, xA is 12
          size_t xB = 6;  // assume we have LJ, 12-6

          fp_type Rij{0}, epsij{0};
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
          } else {
            // we need to calculate the arithmetic mean of Ri and Rj
            Rij = (Rii_i + Rii_j) / fp_type{2};
            // we need to calculate the geometric mean of epsi and epsj
            epsij = std::sqrt(epsii_i * epsii_j);
          }
          if (xA != xB) {
            const fp_type tmp = epsij / (xA - xB);
            const fp_type cA  = tmp * std::pow(Rij, static_cast<fp_type>(xA)) * xB;
            const fp_type cB  = tmp * std::pow(Rij, static_cast<fp_type>(xB)) * xA;

            const fp_type rA = std::pow(distance, static_cast<fp_type>(xA));
            const fp_type rB = std::pow(distance, static_cast<fp_type>(xB));

            e_vdW_Hb = std::min(EINTCLAMP, (cA / rA - cB / rB));
            // TODO check energy smoothing intnbtable.cc:215
            // with NOSQRT to False it seems to be the same as calmping to EINTCLAMP
          } else {
            throw std::runtime_error("ERROR: Exponents must be different, to avoid division by zero!");
          }
        }
        emap_total += e_vdW_Hb;
      }
    }
    return elect_total + emap_total + dmap_total;
  }

} // namespace mudock
