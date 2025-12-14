#pragma once

#include <mudock/cpp_implementation/mutate.hpp>
#include <mudock/cpp_implementation/vectorization.hpp>
#include <mudock/molecule.hpp>
#include <mudock/molecule/fragments.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {
  template<>
  void translate_molecule<cpu_vectorization::GH>(fp_type* __restrict__ x,
                                                 fp_type* __restrict__ y,
                                                 fp_type* __restrict__ z,
                                                 const int num_atoms,
                                                 const fp_type offset_x,
                                                 const fp_type offset_y,
                                                 const fp_type offset_z);

  template<>
  void rotate_molecule<cpu_vectorization::GH>(fp_type* __restrict__ x,
                                              fp_type* __restrict__ y,
                                              fp_type* __restrict__ z,
                                              const int num_atoms,
                                              const fp_type angle_x,
                                              const fp_type angle_y,
                                              const fp_type angle_z);

  template<>
  void rotate_fragment<cpu_vectorization::GH>(fp_type* __restrict__ x,
                                              fp_type* __restrict__ y,
                                              fp_type* __restrict__ z,
                                              const int num_atoms,
                                              const int* __restrict__ frag_mask,
                                              const int start_index,
                                              const int stop_index,
                                              const fp_type angle);
} // namespace mudock
