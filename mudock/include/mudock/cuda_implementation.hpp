#pragma once

#include <mudock/compute/safe_stack.hpp>
#include <mudock/compute/threadpool.hpp>
#include <mudock/knobs.hpp>
#include <mudock/molecule.hpp>
#if defined(MUDOCK_USE_CUDA)
  #include <mudock/cuda_implementation/cuda_manager.cuh>
#else
  #include <mudock/chem/autodock_ligand.hpp>
  #include <mudock/chem/autodock_protein.hpp>
  #include <mudock/log.hpp>
namespace mudock {
  inline void manage_cuda(const std::vector<std::string>&,
                          threadpool&,
                          const knobs,
                          [[maybe_unused]] const autodock_protein& adt_protein,
                          std::shared_ptr<safe_stack<autodock_ligand> >,
                          std::shared_ptr<safe_stack<static_molecule> >) {
    warning("The CUDA implementation is disabled");
  }
} // namespace mudock
#endif
