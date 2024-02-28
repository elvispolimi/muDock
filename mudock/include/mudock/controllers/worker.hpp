#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <mudock/controllers.hpp>
#include <mudock/controllers/squeue.hpp>
#include <mudock/devices.hpp>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {
  // This function is executed by application's thread
  // It carries out the dock and score computation, thus it requires:
  //  - an input queue
  //  - an output queue
  //  - the id of the targeting device to create the context
  // It dequeues from the input batches of ligands, and it enqueues in output the results
  template<kernel_type l_t>
  void worker(std::shared_ptr<squeue<static_molecule>> i_queue,
              std::shared_ptr<squeue<static_molecule>> o_queue,
              const std::size_t id) {
    device_context<l_t> context{
        id}; // Create and set the thread's context, based on the device's ID, and the target language
    // Allocate the scratch memory on the device
    kernel_scratchpad_impl<l_t> scratchpad{context};

    // Start the computation with the first batch
    auto molecules = i_queue->dequeue(get_max_work());
    while (!molecules.empty()) {
      // For each batch copy in, compute, and copy out
      scratchpad.copy_in(molecules);
      dock_and_score(scratchpad, context);
      scratchpad.copy_out(molecules);
      // Enqueue in output the results
      o_queue->enqueue(molecules);
      // Get a new batch of ligands if available
      molecules = i_queue->dequeue(get_max_work());
    }
  };
} // namespace mudock
