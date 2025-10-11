#pragma once

#include <mudock/chem/autodock_ligand.hpp>
#include <mudock/chem/autodock_protein.hpp>
#include <mudock/compute/safe_stack.hpp>
#include <mudock/compute/threadpool.hpp>
#include <mudock/knobs.hpp>
#ifdef MUDOCK_USE_HIP
  #include <mudock/hip_implementation/hip_manager.hpp>
#else
  #include <mudock/log.hpp>
namespace mudock {
  inline void manage_hip(const std::vector<std::string>&,
                         threadpool&,
                         const knobs,
                         [[maybe_unused]] const autodock_protein& adt_protein,
                         std::shared_ptr<safe_stack<autodock_ligand> >,
                         std::shared_ptr<safe_stack<static_molecule> >) {
    warning("The HIP implementation is disabled");
  }
} // namespace mudock
#endif
