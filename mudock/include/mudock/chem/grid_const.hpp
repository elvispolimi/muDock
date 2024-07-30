#pragma once

#include <cstddef>
#include <mudock/type_alias.hpp>

namespace mudock {
  static constexpr fp_type cutoff_distance    = 8.0;
  static constexpr fp_type grid_spacing       = 0.5;  //Angstrom
  static constexpr fp_type covalence_distance = 1.9;  //Angstrom
  
  static constexpr fp_type unknown_distance_1 = 3.61; //Angstrom
  static constexpr fp_type unknown_distance_2 = 1.69; //Angstrom

  static constexpr fp_type unknown_distance_3 = 2.89; //Angstrom
  static constexpr fp_type unknown_distance_4 = 1.69; //Angstrom

  static constexpr std::size_t range_near_atom_receptor = 20;
} // namespace mudock
