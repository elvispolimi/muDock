#pragma once

#include <memory>
#include <mudock/chem/autodock_protein.hpp>
#include <mudock/compute.hpp>
#include <mudock/cpp_implementation/vectorization.hpp>
#include <mudock/cpp_implementation/virtual_screen.hpp>
#include <mudock/knobs.hpp>
#include <mudock/log.hpp>
#include <mudock/molecule.hpp>

namespace mudock {

  template<cpu_vectorization vect>
  class cpp_worker: public worker_interface {
    // this a reference to the input and output queues
    std::shared_ptr<safe_stack<static_molecule>> input_stack;
    std::shared_ptr<safe_stack<static_molecule>> output_stack;

    // this is the functor tha actually implement the virtual screening
    virtual_screen_cpp<vect> virtual_screen;

  public:
    cpp_worker(const knobs knobs,
               const autodock_protein& adt_protein,
               std::shared_ptr<safe_stack<static_molecule>>& input_molecules,
               std::shared_ptr<safe_stack<static_molecule>>& output_molecules,
               const std::size_t cpu_id)
        : input_stack(input_molecules), output_stack(output_molecules), virtual_screen(adt_protein, knobs) {
      cpu_set_t cpuset;
      CPU_ZERO(&cpuset);
      CPU_SET(cpu_id, &cpuset); // Set affinity to the target CPU
      pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
      info("Worker CPP on duty! Set affinity to core ", cpu_id);
    }

    void main() {
      LIKWID_MARKER_REGISTER("GA");

      auto new_ligand = input_stack->dequeue();
      while (new_ligand) {
        virtual_screen(*new_ligand);
        try {
          // TODO check, probably wrong due to the previous move
          output_stack->enqueue(std::move(new_ligand));
        } catch (const std::runtime_error& e) {
          error("Unable to vs molecule ",
                new_ligand->properties.get(property_type::NAME),
                " due to ",
                e.what());
        }
        // NOTE: Clang says "warning: moving a temporary object prevents copy elision"
        // new_ligand = std::move(input_stack->dequeue());
        new_ligand = input_stack->dequeue();
      }
    }
  };
} // namespace mudock
