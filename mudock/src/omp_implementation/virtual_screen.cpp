#include "mudock/molecule/properties.hpp"

#include <alloca.h>
#include <cstddef>
#include <cstring>
#include <cuda_runtime.h>
#include <mudock/cpp_implementation/center_of_mass.hpp>
#include <mudock/cpp_implementation/geometric_transformations.hpp>
#include <mudock/cpp_implementation/mutate.hpp>
#include <mudock/cpp_implementation/weed_bonds.hpp>
// #include <mudock/cuda_implementation/evaluate_fitness.cuh>
#include <mudock/grid.hpp>
#include <mudock/omp_implementation/virtual_screen.hpp>
#include <mudock/utils.hpp>
#include <omp.h>
#include <span>

namespace mudock {
  static constexpr std::size_t max_non_bonds{1024 * 10};

  virtual_screen_omp::virtual_screen_omp(const knobs k,
                                         std::shared_ptr<const grid_atom_mapper> &grid_atom_maps,
                                         std::shared_ptr<const grid_map> &electro_map,
                                         std::shared_ptr<const grid_map> &desolv_map)
      : configuration(k) {
    // TODO move this part into the cuda_worker -> once per GPU
    // Allocate grid maps
    assert(electro_map.get()->index == desolv_map.get()->index &&
           desolv_map.get()->index == grid_atom_maps.get()->get_index());
  }

  void virtual_screen_omp::operator()(batch &incoming_batch) {
#pragma omp target
    for (auto &ligand: std::span(incoming_batch.molecules.data(), incoming_batch.num_ligands)) {
      // TODO
    }
  }
} // namespace mudock
