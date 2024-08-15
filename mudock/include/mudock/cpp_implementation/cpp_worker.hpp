#pragma once

#include <memory>
#include <mudock/compute.hpp>
#include <mudock/cpp_implementation/virtual_screen.hpp>
#include <mudock/molecule.hpp>

namespace mudock {

  class cpp_worker: public worker_interface {
    // this a reference to the input and output queues
    std::shared_ptr<safe_stack<static_molecule>> input_stack;
    std::shared_ptr<safe_stack<static_molecule>> output_stack;

    // this is the functor tha actually implement the virtual screening
    virtual_screen_cpp virtual_screen;

  public:
    // the constructor intialize the kernel and set the CPU affinity to the correct device
    cpp_worker(std::shared_ptr<const dynamic_molecule> protein,
               std::shared_ptr<const grid_atom_mapper> grid_atom_maps,
               std::shared_ptr<const grid_map> electro_map,
               std::shared_ptr<const grid_map> desolv_map,
               std::shared_ptr<safe_stack<static_molecule>> input_molecules,
               std::shared_ptr<safe_stack<static_molecule>> output_molecules,
               const std::size_t cpu_id);

    // this is the thread "main" loop (it will fetch ligands from the queue and compute them)
    void main() override final;
  };

} // namespace mudock
