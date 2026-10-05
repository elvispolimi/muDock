#pragma once

/**
 * @file geom_transform_alpaka.hpp
 * @brief Alpaka specialization for the geometric transformation stage.
 * @details Declares the template specialization of `geom_kernel::operator()` targeting
 *          the Alpaka compute queue, applying translation, rigid orientation, and
 *          torsional angles to molecular candidate conformations.
 */

#include <mudock/alpaka_implementation/buffer_alpaka.hpp>
#include <mudock/alpaka_implementation/queue_alpaka.hpp>
#include <mudock/compute/geometric_transform.hpp>

namespace mudock {
  /**
   * @brief Execution operator specializing geometric transformations for Alpaka.
   * @details Dispatches the `apply_alpaka` kernel onto the accelerator device, transforming
   *          genotype genes (translations, quaternions/Euler rotations, torsion angles)
   *          into 3D atomic coordinates in the receptor binding site frame.
   */
  template<>
  void geom_kernel<queue_alpaka>::operator()();
} // namespace mudock
