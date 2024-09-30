#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <mudock/chem.hpp>
#include <mudock/cpp_implementation/calc_energy_cpp.hpp>
#include <mudock/cpp_implementation/trilinear_interpolation.hpp>
#include <mudock/molecule.hpp>
#include <stdexcept>

namespace mudock {
  fp_type calc_energy(const std::span<fp_type> ligand_x,
                      const std::span<fp_type> ligand_y,
                      const std::span<fp_type> ligand_z,
                      const std::span<fp_type> ligand_vol,
                      const std::span<fp_type> ligand_solpar,
                      const std::span<fp_type> ligand_charge,
                      const std::span<int> ligand_num_hbond,
                      const std::span<fp_type> ligand_Rij_hb,
                      const std::span<fp_type> ligand_Rii,
                      const std::span<fp_type> ligand_epsij_hb,
                      const std::span<fp_type> ligand_epsii,
                      const std::span<autodock_ff> ligand_autodock_type,
                      const int num_atoms,
                      const int n_torsions,
                      const std::span<non_bond_parameter> non_bond_list,
                      const grid_atom_mapper& grid_maps,
                      const grid_map& electro_map,
                      const grid_map& desolv_map) {
    fp_type elect_total_trilinear = 0;
    fp_type emap_total_trilinear  = 0;
    fp_type dmap_total_trilinear  = 0;
    for (int index = 0; index < num_atoms; ++index) {
      const point3D coord{ligand_x[index], ligand_y[index], ligand_z[index]};
      const auto& atom_charge = ligand_charge[index];
      const auto& atom_map    = grid_maps.get_atom_map(ligand_autodock_type[index]);



      if (atom_map.outside_grid(coord)) {
        // printf("Atom %d is outside\n", index);
        const fp_type dist = distance2(coord, grid_maps.get_center());

        const fp_type epenalty = dist * ENERGYPENALTY;
        elect_total_trilinear += epenalty;
        emap_total_trilinear += epenalty;
      } else {
        // Trilinear Interpolationp
        elect_total_trilinear += trilinear_interpolation(electro_map, coord) * atom_charge;
        emap_total_trilinear += trilinear_interpolation(atom_map, coord);
        dmap_total_trilinear += trilinear_interpolation(desolv_map, coord) * std::fabs(atom_charge);
      }
    }

    fp_type elect_total_eintcal{0}, emap_total_eintcal{0}, dmap_total_eintcal{0};
    if (n_torsions > 0) {
      // TODO @Davide suppose that the receptor does not have Flexible residues eintcal.cc:147
      // TODO

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

        // Calcuare desolv
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
          int xA = 12; // for both LJ, 12-6 and HB, 12-10, xA is 12
          int xB = 6;  // assume we have LJ, 12-6

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
            // if (r_smooth > 0) {
            //   const fp_type rlow  = distance - r_smooth / 2;
            //   const fp_type rhigh = distance + r_smooth / 2;
            //   fp_type energy_smooth{100000};
            //   // ((((int)((r)*A_DIV)) > NEINT_1) ? NEINT_1 : ((int)((r)*A_DIV))
            //   for (int j = std::max(fp_type{0}, std::min(rlow * A_DIV, fp_type{NEINT - 1}));
            //        j <= std::min(fp_type{NEINT - 1}, std::min(rhigh * A_DIV, fp_type{NEINT - 1}));
            //        ++j)
            //     energy_smooth = std::min(energy_smooth, e_vdW_Hb);

            //   e_vdW_Hb = energy_smooth;
            // } /* endif smoothing */

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

    const fp_type total_trilinear = emap_total_trilinear + elect_total_trilinear + dmap_total_trilinear;
    const fp_type total_eintcal   = emap_total_eintcal + elect_total_eintcal + dmap_total_eintcal;
    return total_trilinear + total_eintcal + tors_free_energy;
  }

} // namespace mudock
