#include <mudock/cuda_implementation/virtual_screen.cuh>
#include <span>

namespace mudock {
  virtual_screen_cuda::virtual_screen_cuda(const knobs k): dist(fp_type{0}, fp_type{10}), configuration(k) {
    // NOTE: at the moment we don't care about the protein
  }

  void virtual_screen_cuda::operator()(batch& incoming_batch) {
    for (auto& ligand: std::span(incoming_batch.molecules.data(), incoming_batch.num_ligands)) {
      // Reset the random number generator to improve consistency
      generator = std::mt19937{ligand.get()->num_atoms()};
      ligand->properties.assign(property_type::SCORE, std::to_string(dist(generator)));
    }
  }
} // namespace mudock
