#include <cstdlib>
#include <iostream>
#include <mudock/log.hpp>
#include <mudock/sycl_implementation/sycl_worker.hpp>
#include <stdexcept>
#include <string>

namespace mudock {
  sycl_worker::sycl_worker(const knobs knobs,
                           std::shared_ptr<const grid_atom_mapper>& grid_atom_maps,
                           std::shared_ptr<const grid_map>& electro_map,
                           std::shared_ptr<const grid_map>& desolv_map,
                           std::shared_ptr<safe_stack<static_molecule>>& input_molecules,
                           std::shared_ptr<safe_stack<static_molecule>>& output_molecules,
                           std::shared_ptr<reorder_buffer> rb,
                           const std::size_t device_id,
                           sycl::device _dev)
      : input_stack(input_molecules),
        output_stack(output_molecules),
        rob(rb),
        device(_dev),
        queue(device, sycl::property::queue::in_order{}),
        virtual_screen(knobs, grid_atom_maps, electro_map, desolv_map, queue) {
    info("Worker SYCL on duty! Set affinity to device ", device_id);
  }

  void sycl_worker::process(batch& b) {
    try {
      virtual_screen(b);
    } catch (const std::runtime_error& e) { error("Unable to virtual screen a batch due to ", e.what()); }

    for (auto& batch_ligand: std::span(b.molecules.data(), b.num_ligands)) {
      output_stack->enqueue(std::move(batch_ligand));
    }
  }

  void sycl_worker::main() {
    // process the input ligands
    auto new_ligand = input_stack->dequeue();
    while (new_ligand) {
      auto [new_batch, is_valid] = rob->add_ligand(std::move(new_ligand));
      if (is_valid) {
        process(new_batch);
      }
      // NOTE: Clang says "warning: moving a temporary object prevents copy elision"
      // new_ligand = std::move(input_stack->dequeue());
      new_ligand = input_stack->dequeue();
    }

    // finish the half empty batches in the rob
    auto rob_is_empty = false;
    while (!rob_is_empty) {
      auto [half_batch, is_valid] = rob->flush_one();
      if (is_valid) {
        process(half_batch);
      } else {
        rob_is_empty = true;
      }
    }
  }

} // namespace mudock
