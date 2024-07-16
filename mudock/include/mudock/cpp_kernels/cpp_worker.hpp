#pragma once

#include <mudock/compute.hpp>
#include <mudock/cpp_kernels/virtual_screen.hpp>
#include <mudock/molecule.hpp>

namespace mudock {

  class cpp_worker {
    // this a reference to the work pool
    safe_stack<static_molecule>& input_stack;
    safe_stack<static_molecule>& output_stack;

    // this is the functor tha actually implement the virtual screening
    virtual_screen_cpp virtual_screen;

  public:
    // the constructor intialize the kernel and set the CPU affinity to the correct device
    cpp_worker(std::shared_ptr<dynamic_molecule>& protein,
               safe_stack<static_molecule>& input_molecules,
               safe_stack<static_molecule>& output_molecules,
               const std::size_t cpu_id);

    // this is the thread "main" loop (it will fetch ligands from the queue and compute them)
    void main();
  };

} // namespace mudock
