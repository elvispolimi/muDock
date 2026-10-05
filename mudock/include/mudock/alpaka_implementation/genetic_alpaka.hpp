#pragma once

/**
 * @file genetic_alpaka.hpp
 * @brief Alpaka specialization for the Genetic Algorithm (GA) stage.
 * @details Declares the template specializations of `genetic_kernel` methods
 *          (initialize, operator(), finalize) targeting `queue_alpaka` to execute
 *          conformation population evolution on the accelerator device.
 */

#include <mudock/alpaka_implementation/buffer_alpaka.hpp>
#include <mudock/alpaka_implementation/queue_alpaka.hpp>
#include <mudock/compute/genetic.hpp>

namespace mudock {
  /**
   * @brief Initializes the candidate population on the accelerator device.
   * @details Dispatches the `initialize_alpaka` kernel to generate initial pseudo-random
   *          conformations (translations, rotations, and torsions) using the Alpaka PRNG engine.
   */
  template<>
  void genetic_kernel<queue_alpaka>::initialize();

  /**
   * @brief Advances the genetic algorithm population by one generation.
   * @details Dispatches the `iterate_alpaka` kernel to execute tournament selection,
   *          single-point crossover, and mutation operations on candidate chromosomes in parallel.
   */
  template<>
  void genetic_kernel<queue_alpaka>::operator()();

  /**
   * @brief Finalizes the genetic search stage on the accelerator device.
   * @details Dispatches the `finalize_alpaka` kernel to locate the best-scoring candidate
   *          pose across the evaluated population and copy final coordinates to host memory.
   */
  template<>
  void genetic_kernel<queue_alpaka>::finalize();
} // namespace mudock
