#include <cstdlib>
#include <cuda_runtime.h>
#include <iostream>
#include <mudock/cuda_kernels/cuda_check_error_macro.cuh>
#include <mudock/cuda_kernels/cuda_worker.hpp>
#include <mudock/log.hpp>
#include <stdexcept>
#include <string>

namespace mudock {
  cuda_worker::cuda_worker(std::shared_ptr<dynamic_molecule> protein,
                           std::shared_ptr<safe_stack<static_molecule>> input_molecules,
                           std::shared_ptr<safe_stack<static_molecule>> output_molecules,
                           std::shared_ptr<reorder_buffer> rb,
                           const std::size_t gpu_id)
      : input_stack(input_molecules), output_stack(output_molecules), rob(rb), virtual_screen(protein) {
    MUDOCK_CHECK(cudaSetDevice(static_cast<int>(gpu_id)));
    info("Worker CUDA on duty! Set affinity to GPU ", gpu_id);
  }

  void cuda_worker::main() {
    auto new_ligand = input_stack->dequeue();
    while (new_ligand) {
      auto [new_batch, is_valid] = rob->add_ligand(std::move(new_ligand));
      if (is_valid) {
        try {
          virtual_screen(std::span(new_batch.molecules));
        } catch (const std::runtime_error& e) {
          std::cerr << std::string{"Unable to vs a batch of ligands due to "} + e.what() + std::string{"\n"};
        }

        for (auto& batch_ligand: new_batch.molecules) { output_stack->enqueue(std::move(batch_ligand)); }
        new_ligand = std::move(input_stack->dequeue());
      }
    }
  }

} // namespace mudock
