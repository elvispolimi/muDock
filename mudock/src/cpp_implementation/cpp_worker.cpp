#include <cstdlib>
#include <iostream>
#include <mudock/cpp_implementation/cpp_worker.hpp>
#include <mudock/log.hpp>
#include <stdexcept>
#include <string>

namespace mudock {
  cpp_worker::cpp_worker(std::shared_ptr<const dynamic_molecule> protein,
                         std::shared_ptr<const grid_atom_mapper> grid_atom_maps,
                         std::shared_ptr<const grid_map> electro_map,
                         std::shared_ptr<const grid_map> desolv_map,
                         std::shared_ptr<safe_stack<static_molecule>> input_molecules,
                         std::shared_ptr<safe_stack<static_molecule>> output_molecules,
                         const std::size_t cpu_id)
      : input_stack(input_molecules),
        output_stack(output_molecules),
        virtual_screen(protein, grid_atom_maps, electro_map, desolv_map) {
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(cpu_id, &cpuset); // Set affinity to the target CPU
    pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
    info("Worker CPP on duty! Set affinity to core ", cpu_id);
  }

  void cpp_worker::main() {
    auto new_ligand = input_stack->dequeue();
    while (new_ligand) {
      virtual_screen(*new_ligand);
      try {
        output_stack->enqueue(std::move(new_ligand));
      } catch (const std::runtime_error& e) {
        error("Unable to vs molecule ",
              new_ligand->properties.get(property_type::NAME),
              " due to ",
              e.what());
      }
      new_ligand = std::move(input_stack->dequeue());
    }
  }

} // namespace mudock
