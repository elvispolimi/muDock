#include <cstdlib>
#include <cuda_runtime.h>
#include <iostream>
#include <mudock/cuda_implementation/cuda_check_error_macro.cuh>
#include <mudock/cuda_implementation/cuda_worker.hpp>
#include <mudock/log.hpp>
#include <stdexcept>
#include <string>

namespace mudock {
  cuda_worker::cuda_worker(const knobs knobs,
                           std::shared_ptr<const grid_atom_mapper>& grid_atom_maps,
                           std::shared_ptr<const grid_map>& electro_map,
                           std::shared_ptr<const grid_map>& desolv_map,
                           std::shared_ptr<safe_stack<static_molecule>>& input_molecules,
                           std::shared_ptr<safe_stack<static_molecule>>& output_molecules,
                           std::shared_ptr<reorder_buffer> rb,
                           const std::size_t gpu_id)
      : input_stack(input_molecules),
        output_stack(output_molecules),
        rob(rb),
        virtual_screen(knobs, grid_atom_maps, electro_map, desolv_map) {
    MUDOCK_CHECK(cudaSetDevice(static_cast<int>(gpu_id)));
    info("Worker CUDA on duty! Set affinity to GPU ", gpu_id);
  }

  void cuda_worker::process(batch& b) {
    try {
      virtual_screen(b);
    } catch (const std::runtime_error& e) { error("Unable to virtual screen a batch due to ", e.what()); }

    for (auto& batch_ligand: std::span(b.molecules.data(), b.num_ligands)) {
      output_stack->enqueue(std::move(batch_ligand));
    }
  }

  void cuda_worker::main() {
    // process the input ligands
    auto new_ligand = input_stack->dequeue();
    while (new_ligand) {
      auto [new_batch, is_valid] = rob->add_ligand(std::move(new_ligand));
      if (is_valid) {
        process(new_batch);
        new_ligand = std::move(input_stack->dequeue());
      }
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
