#pragma once

#include <cstddef>
#include <mudock/type_alias.hpp>

namespace mudock {
  static constexpr fp_type cutoff_distance = static_cast<fp_type>(8.0);
  static constexpr fp_type grid_spacing    = static_cast<fp_type>(0.5);    //Angstrom
  static constexpr fp_type inv_spacing{1 / grid_spacing};                  // Angstrom
  static constexpr fp_type covalence_distance = static_cast<fp_type>(1.9); //Angstrom

  static constexpr fp_type unknown_distance_1 = static_cast<fp_type>(3.61); //Angstrom
  static constexpr fp_type unknown_distance_2 = static_cast<fp_type>(1.69); //Angstrom

  static constexpr fp_type unknown_distance_3 = static_cast<fp_type>(2.89); //Angstrom
  static constexpr fp_type unknown_distance_4 = static_cast<fp_type>(1.69); //Angstrom

  static constexpr int range_near_atom_receptor = 20;

  static constexpr int NEINT{2048};
  static constexpr int NDIEL{16384};

  static constexpr fp_type NBC =
      static_cast<fp_type>(8.00); /* Non-bonded cutoff for internal energy calc./Ang*/

  /* Used in distance look-up table. i.e. every 1/100-th of an Angstrom */
  static constexpr fp_type A_DIV = static_cast<fp_type>(100.00); /* Used in distance look-up table. */
  static constexpr fp_type EINTCLAMP =
      static_cast<fp_type>(100000); /* Clamp pairwise internal energies (kcal/mol )  */

  static constexpr fp_type factor =
      static_cast<fp_type>(332.0); /* Used to convert between calories and SI units */

  static constexpr fp_type solpar_q = static_cast<fp_type>(0.01097);
  static constexpr fp_type sigma    = static_cast<fp_type>(3.6);
  static constexpr fp_type sigma_square{sigma * sigma};

  static constexpr fp_type precision = static_cast<fp_type>(0.0001);

  static constexpr fp_type RMIN_ELEC = static_cast<fp_type>(0.5);
  static constexpr fp_type RMIN_ELEC_SQUARE{RMIN_ELEC * RMIN_ELEC};
  static constexpr fp_type ELECSCALE = static_cast<fp_type>(332.06363);
  static constexpr fp_type qsolpar   = static_cast<fp_type>(0.01097);
  static constexpr fp_type r_smooth  = static_cast<fp_type>(
      0.5); // vdw nonbond smoothing range, not radius, Ang - default 0.5 matches AutoGrid recommendations
  // Non bond cutoff
  static constexpr fp_type nbc2{64}; // 8*8
  static constexpr fp_type ENERGYPENALTY{500};

  // Lennard-Jones
  static constexpr int xA_default = 12;
  static constexpr int xB_default = 6;
} // namespace mudock
