#pragma once

#include <memory>
#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/compute/safe_stack.hpp>
#include <mudock/compute/stage.hpp>
#include <mudock/compute/threadpool.hpp>
#include <mudock/cpp_implementation/vectorization.hpp>
#include <mudock/knobs.hpp>
#include <mudock/likwid_utils.hpp>
#include <mudock/log.hpp>
#include <mudock/molecule.hpp>

namespace mudock {

  template<typename stage_t>
  class worker: public worker_interface {
    // this a reference to the input and output queues
    std::shared_ptr<safe_stack<static_molecule>> input_stack;
    std::shared_ptr<safe_stack<static_molecule>> output_stack;

    // this is a reoder buffer that we can use to fetch baches of ligands out of order
    std::shared_ptr<reorder_buffer<static_molecule>> rob;

    // this is the functor tha actually implement the virtual screening
    stage_t pipeline;

    void process(batch<static_molecule>& b) {
      try {
        pipeline.prepare(b);
        pipeline();
        pipeline.teardown(b);
      } catch (const std::runtime_error& e) { error("Unable to virtual screen a batch due to ", e.what()); }

      for (auto& batch_ligand: std::span(b.molecules.data(), b.num_ligands)) {
        output_stack->enqueue(std::make_unique<static_molecule>(std::move(*batch_ligand)));
      }
    }

  public:
    worker(std::shared_ptr<safe_stack<static_molecule>>& input_molecules,
           std::shared_ptr<safe_stack<static_molecule>>& output_molecules,
           std::shared_ptr<reorder_buffer<static_molecule>> rb,
           // const std::size_t cpu_id,
           stage_t&& _pipeline)
        : input_stack(input_molecules),
          output_stack(output_molecules),
          rob(rb),
          pipeline(std::move(_pipeline)) {}

    void main() {
      LIKWID_MARKER_REGISTER("GA");

      // process the input ligands
      auto new_ligand = input_stack->dequeue();
      while (new_ligand) {
        auto [new_batch, is_valid] = rob->add_ligand(std::move(new_ligand));
        if (is_valid) {
          process(new_batch);
        }
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
  };
} // namespace mudock
