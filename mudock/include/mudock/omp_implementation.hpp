#pragma once

#include <mudock/chem/autodock_protein.hpp>
#include <mudock/compute/safe_stack.hpp>
#include <mudock/compute/threadpool.hpp>
#include <mudock/knobs.hpp>
#ifdef MUDOCK_USE_OMP
  #include <mudock/omp_implementation/omp_batch_sizer.hpp>
  #include <mudock/omp_implementation/omp_manager.hpp>
  #include <mudock/omp_implementation/virtual_screen.hpp>
#else
  #include <mudock/log.hpp>
namespace mudock {
  inline void manage_omp(const std::vector<std::string>&,
                         threadpool&,
                         const knobs,
                         [[maybe_unused]] const autodock_protein& adt_protein,
                         std::shared_ptr<safe_stack<static_molecule> >,
                         std::shared_ptr<safe_stack<static_molecule> >) {
    warning("The OpenMP implementation is disabled");
  }
} // namespace mudock
#endif
