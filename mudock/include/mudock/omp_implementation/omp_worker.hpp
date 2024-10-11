#pragma once

#include <memory>
#include <mudock/compute.hpp>
#include <mudock/knobs.hpp>
#include <mudock/molecule.hpp>
#include <mudock/omp_implementation/virtual_screen.hpp>

namespace mudock {

  class omp_worker: public worker_interface {
    // this a reference to the input and output queues
    std::shared_ptr<safe_stack<static_molecule>> input_stack;
    std::shared_ptr<safe_stack<static_molecule>> output_stack;

    // this is a reoder buffer that we can use to fetch baches of ligands out of order
    std::shared_ptr<reorder_buffer> rob;

    // this is the functor tha actually implement the virtual screening
    virtual_screen_omp virtual_screen;

    // utility function to process a single batch of ligands
    void process(batch& b);

  public:
    // the constructor intialize the kernel and set the device affinity to the correct device
    omp_worker(const knobs knobs,
               std::shared_ptr<const grid_atom_mapper>& grid_atom_maps,
               std::shared_ptr<const grid_map>& electro_map,
               std::shared_ptr<const grid_map>& desolv_map,
               std::shared_ptr<safe_stack<static_molecule>>& input_molecules,
               std::shared_ptr<safe_stack<static_molecule>>& output_molecules,
               std::shared_ptr<reorder_buffer> rb,
               const std::size_t device_id);

    // this is the thread "main" loop (it will fetch ligands from the queue and compute them)
    void main() override final;
  };

} // namespace mudock
