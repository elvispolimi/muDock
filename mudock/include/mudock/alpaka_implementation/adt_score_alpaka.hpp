#pragma once

/**
 * @file adt_score_alpaka.hpp
 * @brief Alpaka specialization for the AutoDock (ADT) scoring stage.
 * @details Declares template specializations of batch sizing utilities and the
 *          adt_score_kernel execution operator targeting `queue_alpaka`.
 */

#include <memory>
#include <mudock/alpaka_implementation/buffer_alpaka.hpp>
#include <mudock/alpaka_implementation/queue_alpaka.hpp>
#include <mudock/compute/adt_score.hpp>

namespace mudock {
  /**
   * @brief Computes batch alignment constraints for the ADT scoring kernel on Alpaka.
   * @param atoms Number of atoms in the candidate ligand molecule.
   * @param q_b Shared pointer to the Alpaka queue managing device execution.
   * @return batch_multiple Struct specifying block and grid dimension alignments.
   */
  template<>
  batch_multiple get_adt_score_batch_multiple<queue_alpaka>(const int atoms, std::shared_ptr<queue_alpaka> q_b);

  /**
   * @brief Kernel execution operator specializing ADT molecular scoring for Alpaka.
   * @details Computes electrostatic, van der Waals, and desolvation energies for candidate
   *          ligand poses using 3D trilinear interpolation over receptor grid maps.
   */
  template<>
  void adt_score_kernel<queue_alpaka>::operator()();
} // namespace mudock
