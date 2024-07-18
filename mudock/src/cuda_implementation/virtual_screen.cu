#include <mudock/cuda_implementation/virtual_screen.cuh>

namespace mudock {
  virtual_screen_cuda::virtual_screen_cuda([[maybe_unused]] std::shared_ptr<dynamic_molecule>& protein)
      : generator(protein->num_atoms()), dist(fp_type{0}, fp_type{10}) {
    // NOTE: at the moment we don't care about the protein
  }

  void virtual_screen_cuda::operator()(batch& incoming_batch) {
    for (auto& ligand: incoming_batch.molecules) {
      ligand->properties.assign(property_type::SCORE, std::to_string(dist(generator)));
    }
  }
} // namespace mudock
