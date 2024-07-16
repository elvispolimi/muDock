#include <cstdlib>
#include <mudock/cpp_kernels/cpp_worker.hpp>

namespace mudock {
  cpp_worker::cpp_worker(std::shared_ptr<dynamic_molecule>& protein,
                         safe_stack<static_molecule>& input_molecules,
                         safe_stack<static_molecule>& output_molecules,
                         const std::size_t cpu_id)
      : input_stack(input_molecules), output_stack(output_molecules), virtual_screen(protein) {
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(cpu_id, &cpuset); // Set affinity to the target CPU
    pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
  }

  void cpp_worker::main() {
    auto new_ligand = input_stack.dequeue();
    while (new_ligand) {
      virtual_screen(*new_ligand);
      output_stack.enqueue(std::move(new_ligand));
      new_ligand = std::move(input_stack.dequeue());
    }
  }

} // namespace mudock
