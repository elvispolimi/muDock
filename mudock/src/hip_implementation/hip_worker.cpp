#include <cstdlib>
#include <iostream>
#include <mudock/hip_implementation/hip_check_error_macro.hpp>
#include <mudock/hip_implementation/hip_worker.hpp>
#include <mudock/log.hpp>
#include <stdexcept>
#include <string>

namespace mudock {
  hip_worker::hip_worker(const knobs knobs,
                         std::shared_ptr<safe_stack<autodock_ligand>>& input_molecules,
                         std::shared_ptr<safe_stack<static_molecule>>& output_molecules,
                         std::shared_ptr<reorder_buffer<autodock_ligand>> rb,
                         const std::shared_ptr<const device> dev)
      : input_stack(input_molecules), output_stack(output_molecules), rob(rb), virtual_screen(knobs, dev) {}

  void hip_worker::process(batch<autodock_ligand>& b) {
    try {
      virtual_screen(b);
    } catch (const std::runtime_error& e) { error("Unable to virtual screen a batch due to ", e.what()); }

    for (auto& batch_ligand: std::span(b.molecules.data(), b.num_ligands)) {
      output_stack->enqueue(std::make_unique<static_molecule>(std::move(*batch_ligand)));
    }
  }

  void hip_worker::main() {
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
