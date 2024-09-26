#include "mudock/chem/autodock_parameters.hpp"
#include "mudock/chem/grid_const.hpp"
#include "mudock/type_alias.hpp"

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
  static constexpr fp_type r_smooth{
      0.0}; // vdw nonbond smoothing range, not radius, Ang - default 0.5 matches AutoGrid recommendations
  // Non bond cutoff
  static constexpr fp_type nbc2{64}; // 8*8
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
    const point3D coordinates_in_grid = map.get_index_from_coordinates(coord);
    const point3D coordinates_floored = point3D{std::floor(coordinates_in_grid.x),
                                                std::floor(coordinates_in_grid.y),
                                                std::floor(coordinates_in_grid.z)};

    const fp_type delta_low_x = coordinates_in_grid.x - coordinates_floored.x;
    const fp_type delta_low_y = coordinates_in_grid.y - coordinates_floored.y;
    const fp_type delta_low_z = coordinates_in_grid.z - coordinates_floored.z;

    std::array<fp_type, 2> px{fp_type{1} - delta_low_x, delta_low_x};
    std::array<fp_type, 2> py{fp_type{1} - delta_low_y, delta_low_y};
    std::array<fp_type, 2> pz{fp_type{1} - delta_low_z, delta_low_z};
    fp_type value{0};
    for (size_t i = 0; i <= 1; ++i)
      for (size_t j = 0; j <= 1; ++j)
        for (size_t t = 0; t <= 1; ++t) {
          value += px[t] * py[j] * pz[i] *
                   map.at(coordinates_floored.x + t, coordinates_floored.y + j, coordinates_floored.z + i);
        }
    return value;
  }

  // nonbonds.cc for nbmatrix required by weed_bonds
  void nonbonds(grid<uint_fast8_t, index2D>& nbmatrix,
                const std::span<const bond> ligand_bond,
                const std::size_t num_atoms) {
    //
    // in "nbmatrix", the values 1 (and 4) mean this pair of atoms will be included in the internal, non-bonded list
    //                           0                                         ignored

    // set all nonbonds in nbmatrix to 1, except "1-1 interactions" (self interaction)
    for (size_t i = 0; i < num_atoms; i++) {
      for (size_t j = 0; j < num_atoms; j++) { nbmatrix.at(i, j) = 1; } // j
      nbmatrix.at(i, i) = 0;                                            /* 2005-01-10 RH & GMM */
    }

    for (auto& bond: ligand_bond) {
      // Ignore 1-2 Interactions
      nbmatrix.at(bond.source, bond.dest) = 0;
      nbmatrix.at(bond.dest, bond.source) = 0;
    }

    for (auto& bond_1: ligand_bond)
      for (auto& bond_2: ligand_bond) { // loop over each atom "k" bonded to the current atom "j"
        std::size_t outer_1{0};
        std::size_t outer_2{0};
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
          std::size_t outer_3{0};
          std::size_t outer_4{0};
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
                  const std::size_t num_atoms,
                  const fragments<static_containers>& ligand_fragments) {
    const auto& ligand_rigid_pieces = ligand_fragments.get_rigid_pieces();

    for (size_t j = 0; j < num_atoms; ++j) {
      for (size_t i = 0; i < num_atoms; ++i) {
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
      for (size_t k = 0; k < num_atoms; ++k) {
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
    for (size_t i = 0; i < num_atoms; ++i) {
      for (size_t j = i + 1; j < num_atoms; ++j) {
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

  fp_type calc_energy(const std::span<fp_type> ligand_x,
                      const std::span<fp_type> ligand_y,
                      const std::span<fp_type> ligand_z,
                      const std::span<fp_type> ligand_vol,
                      const std::span<fp_type> ligand_solpar,
                      const std::span<fp_type> ligand_charge,
                      const std::span<std::size_t> ligand_num_hbond,
                      const std::span<fp_type> ligand_Rij_hb,
                      const std::span<fp_type> ligand_Rii,
                      const std::span<fp_type> ligand_epsij_hb,
                      const std::span<fp_type> ligand_epsii,
                      const std::span<autodock_ff> ligand_autodock_type,
                      const std::span<const bond> ligand_bond,
                      const std::size_t num_atoms,
                      const fragments<static_containers>& ligand_fragments,
                      const grid_atom_mapper& grid_maps,
                      const grid_map& electro_map,
                      const grid_map& desolv_map) {
    fp_type elect_total_trilinear = 0;
    fp_type emap_total_trilinear  = 0;
    fp_type dmap_total_trilinear  = 0;
    for (size_t index = 0; index < num_atoms; ++index) {
      const point3D coord{ligand_x[index], ligand_y[index], ligand_z[index]};
      const auto& atom_charge = ligand_charge[index];
      const auto& atom_map    = grid_maps.get_atom_map(ligand_autodock_type[index]);

      if (atom_map.outside_grid(coord)) {
        const fp_type dist = distance2(coord, grid_maps.get_center());

        const fp_type epenalty = dist * ENERGYPENALTY;
        elect_total_trilinear += epenalty;
        emap_total_trilinear += epenalty;
      } else {
        // Trilinear Interpolation
        elect_total_trilinear += trilinear_interpolation(electro_map, coord) * atom_charge;
        emap_total_trilinear += trilinear_interpolation(atom_map, coord);
        dmap_total_trilinear += trilinear_interpolation(desolv_map, coord) * std::fabs(atom_charge);
      }
    }

    const int n_torsions = ligand_fragments.get_num_rotatable_bonds();
    fp_type elect_total_eintcal{0}, emap_total_eintcal{0}, dmap_total_eintcal{0};
    if (n_torsions > 0) {
      // TODO @Davide suppose that the receptor does not have Flexible residues eintcal.cc:147
      // TODO
      grid<uint_fast8_t, index2D> nbmatrix{{num_atoms, num_atoms}};
      nonbonds(nbmatrix, ligand_bond, num_atoms);
      std::vector<non_bond_parameter> non_bond_list;
      weed_bonds(nbmatrix, non_bond_list, num_atoms, ligand_fragments);

      for (const auto& non_bond: non_bond_list) {
        const int& a1 = non_bond.a1;
        const int& a2 = non_bond.a2;

        const fp_type distance_two       = distance2(point3D{ligand_x[a1], ligand_y[a1], ligand_z[a1]},
                                               point3D{ligand_x[a2], ligand_y[a2], ligand_z[a2]});
        const fp_type distance_two_clamp = std::clamp(distance_two, RMIN_ELEC * RMIN_ELEC, distance_two);
        const fp_type distance           = std::sqrt(distance_two_clamp);

        //  Calculate  Electrostatic  Energy
        const fp_type r_dielectric = fp_type{1} / (distance * calc_ddd_Mehler_Solmajer(distance));
        const fp_type e_elec       = ligand_charge[non_bond.a1] * ligand_charge[non_bond.a2] * ELECSCALE *
                               autodock_parameters::coeff_estat * r_dielectric;
        elect_total_eintcal += e_elec;

        const fp_type nb_desolv =
            (ligand_vol[a2] * (ligand_solpar[a1] + qsolpar * std::fabs(ligand_charge[a1])) +
             ligand_vol[a1] * (ligand_solpar[a2] + qsolpar * std::fabs(ligand_charge[a2])));

        const fp_type e_desolv = autodock_parameters::coeff_desolv *
                                 std::exp(fp_type{-0.5} / (sigma * sigma) * distance_two_clamp) * nb_desolv;
        dmap_total_eintcal += e_desolv;
        fp_type e_vdW_Hb{0};
        if (distance_two_clamp < nbc2) {
          //  Find internal energy parameters, i.e.  epsilon and r-equilibrium values...
          //  Lennard-Jones and Hydrogen Bond Potentials
          // This can be precomputed as in intnbtable.cc
          const auto& hbond_i    = ligand_num_hbond[a1];
          const auto& hbond_j    = ligand_num_hbond[a2];
          const auto& Rij_hb_i   = ligand_Rij_hb[a1];
          const auto& Rij_hb_j   = ligand_Rij_hb[a2];
          const auto& Rii_i      = ligand_Rii[a1];
          const auto& Rii_j      = ligand_Rii[a2];
          const auto& epsij_hb_i = ligand_epsij_hb[a1];
          const auto& epsij_hb_j = ligand_epsij_hb[a2];
          const auto& epsii_i    = ligand_epsii[a1];
          const auto& epsii_j    = ligand_epsii[a2];

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

            /* smooth with min function; 
              r_smooth is Angstrom range of "smoothing" */
            // TODO rework considering this precomputing with other values
            if (r_smooth > 0) {
              const fp_type rlow  = distance - r_smooth / 2;
              const fp_type rhigh = distance + r_smooth / 2;
              fp_type energy_smooth{100000};
              // ((((int)((r)*A_DIV)) > NEINT_1) ? NEINT_1 : ((int)((r)*A_DIV))
              for (std::size_t j = std::max(fp_type{0}, std::min(rlow * A_DIV, fp_type{NEINT - 1}));
                   j <= std::min(fp_type{NEINT - 1}, std::min(rhigh * A_DIV, fp_type{NEINT - 1}));
                   ++j)
                energy_smooth = std::min(energy_smooth, e_vdW_Hb);

              e_vdW_Hb = energy_smooth;
            } /* endif smoothing */

            // TODO check energy smoothing intnbtable.cc:215
            // with NOSQRT to False it seems to be the same as calmping to EINTCLAMP
          } else {
            throw std::runtime_error("ERROR: Exponents must be different, to avoid division by zero!");
          }
        }
        emap_total_eintcal += e_vdW_Hb;
      }
    }
    const fp_type tors_free_energy = n_torsions * autodock_parameters::coeff_tors;

    const fp_type total_trilnear = emap_total_trilinear + elect_total_trilinear + dmap_total_trilinear;
    const fp_type total_eintcal = emap_total_eintcal + elect_total_eintcal + dmap_total_eintcal;

    return total_trilnear + total_eintcal + tors_free_energy;
  }

} // namespace mudock
