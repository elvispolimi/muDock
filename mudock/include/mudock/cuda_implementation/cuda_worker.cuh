#pragma once

#include <memory>
#include <mudock/chem/autodock_ligand.hpp>
#include <mudock/compute.hpp>
#include <mudock/cuda_implementation/virtual_screen.cuh>
#include <mudock/knobs.hpp>
#include <mudock/molecule.hpp>

namespace mudock {
  class cuda_worker: public worker_interface {
    // this a reference to the input and output queues
    std::shared_ptr<safe_stack<autodock_ligand>> input_stack;
    std::shared_ptr<safe_stack<static_molecule>> output_stack;

    // this is a reoder buffer that we can use to fetch baches of ligands out of order
    std::shared_ptr<reorder_buffer<autodock_ligand>> rob;

    // this is the functor tha actually implement the virtual screening
    virtual_screen_cuda virtual_screen;

    // utility function to process a single batch of ligands
    void process(batch<autodock_ligand>& b);

  public:
    // the constructor intialize the kernel and set the GPU affinity to the correct device
    cuda_worker(const knobs knobs,
                std::shared_ptr<safe_stack<autodock_ligand>>& input_molecules,
                std::shared_ptr<safe_stack<static_molecule>>& output_molecules,
                std::shared_ptr<reorder_buffer<autodock_ligand>> rb,
                std::shared_ptr<device> dev);

    // this is the thread "main" loop (it will fetch ligands from the queue and compute them)
    void main() override final;
  };

} // namespace mudock
