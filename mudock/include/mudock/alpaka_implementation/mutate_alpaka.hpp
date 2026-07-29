#pragma once

#include <alpaka/alpaka.hpp>
#include <mudock/chem/atom.hpp>
#include <mudock/chem/element.hpp>
#include <mudock/chem/residue.hpp>
#include <mudock/math/vector.hpp>

namespace mudock {

  template<typename TAcc>
  ALPAKA_FN_ACC inline void translate_molecule_alpaka(TAcc const& acc,
                                                      const fp_type* const __restrict__ x,
                                                      const fp_type* const __restrict__ y,
                                                      const fp_type* const __restrict__ z,
                                                      fp_type* const __restrict__ x_out,
                                                      fp_type* const __restrict__ y_out,
                                                      fp_type* const __restrict__ z_out,
                                                      const math::vector3_t<fp_type>& translation,
                                                      const int start_index,
                                                      const int num_atoms) {
    // TODO: qui andrà copiato il corpo del ciclo "for" di translate_molecule
    // attualmente presente in geom_transform_alpaka.cpp
  }

  template<typename TAcc>
  ALPAKA_FN_ACC inline void rotate_molecule_alpaka(/* TODO: Firma della funzione */) {
    // TODO: copia il corpo da geom_transform_alpaka.cpp
  }

  template<typename TAcc>
  ALPAKA_FN_ACC inline void rotate_fragment_alpaka(/* TODO: Firma della funzione */) {
    // TODO: copia il corpo da geom_transform_alpaka.cpp
  }

} // namespace mudock
