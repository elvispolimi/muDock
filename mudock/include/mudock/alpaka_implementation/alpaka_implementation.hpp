#pragma once

/**
 * @file alpaka_implementation.hpp
 * @brief Umbrella header aggregating all Alpaka backend headers.
 * @details Included by muDock core engine when `MUDOCK_USE_ALPAKA` is enabled.
 *          Brings in queue, buffer, kernel invocation, and stage specializations.
 */

#ifdef MUDOCK_USE_ALPAKA
  #include <mudock/alpaka_implementation/buffer_alpaka.hpp>
  #include <mudock/alpaka_implementation/queue_alpaka.hpp>
  #include <mudock/alpaka_implementation/adt_score_alpaka.hpp>
  #include <mudock/alpaka_implementation/genetic_alpaka.hpp>
  #include <mudock/alpaka_implementation/geom_transform_alpaka.hpp>
  #include <mudock/alpaka_implementation/invoke_kernel_alpaka.hpp>
#endif
