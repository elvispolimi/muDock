#pragma once

#include <cstddef>
#include <mudock/type_alias.hpp>

namespace mudock {
  static constexpr fp_type cutoff_distance    = 8.0;
  static constexpr fp_type grid_spacing       = 0.5; //Angstrom
  static constexpr fp_type covalence_distance = 1.9; //Angstrom

  static constexpr fp_type unknown_distance_1 = 3.61; //Angstrom
  static constexpr fp_type unknown_distance_2 = 1.69; //Angstrom

  static constexpr fp_type unknown_distance_3 = 2.89; //Angstrom
  static constexpr fp_type unknown_distance_4 = 1.69; //Angstrom

  static constexpr std::size_t range_near_atom_receptor = 20;

  // TODO why this value?
  static constexpr size_t NEINT{
      2048}; /* Number of values in internal energy table, they are based on radius range values */
  /* Number of dielectric and desolvation values in lookup table.
  NDIEL is bigger than NEINT because electrostatic interactions are much
  longer-range than van der Waals interactions. */
  // TODO wrong comment in autogrid
  static constexpr size_t NDIEL{16384};

  /* Used in distance look-up table. i.e. every 1/100-th of an Angstrom */
  static constexpr fp_type A_DIV{100.00};     /* Used in distance look-up table. */
  static constexpr fp_type EINTCLAMP{100000}; /* Clamp pairwise internal energies (kcal/mol )  */

  static constexpr fp_type factor{332.0}; /* Used to convert between calories and SI units */

  static constexpr fp_type solpar_q{0.01097};
  static constexpr fp_type sigma{3.6};
} // namespace mudock
