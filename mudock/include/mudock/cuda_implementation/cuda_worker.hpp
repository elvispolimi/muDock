#pragma once

#include <memory>
#include <mudock/compute.hpp>
#include <mudock/cuda_implementation/virtual_screen.cuh>
#include <mudock/molecule.hpp>

namespace mudock {

  class cuda_worker: public worker_interface {
    // this a reference to the input and output queues
    std::shared_ptr<safe_stack<static_molecule>> input_stack;
    std::shared_ptr<safe_stack<static_molecule>> output_stack;

    // this is a reoder buffer that we can use to fetch baches of ligands out of order
    std::shared_ptr<reorder_buffer> rob;

    // this is the functor tha actually implement the virtual screening
    virtual_screen_cuda virtual_screen;

    // utility function to process a single batch of ligands
    void process(batch b);

  public:
    // the constructor intialize the kernel and set the CPU affinity to the correct device
    cuda_worker(std::shared_ptr<dynamic_molecule> protein,
                std::shared_ptr<safe_stack<static_molecule>> input_molecules,
                std::shared_ptr<safe_stack<static_molecule>> output_molecules,
                std::shared_ptr<reorder_buffer> rb,
                const std::size_t gpu_id);

    // this is the thread "main" loop (it will fetch ligands from the queue and compute them)
    void main() override final;
  };

} // namespace mudock
