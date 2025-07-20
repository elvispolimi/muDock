#pragma once

#include <mudock/chem/autodock_protein.hpp>
#include <mudock/compute/safe_stack.hpp>
#include <mudock/compute/threadpool.hpp>
#include <mudock/knobs.hpp>
#ifdef MUDOCK_USE_SYCL
  #include <mudock/sycl_implementation/sycl_manager.hpp>
#else
  #include <mudock/log.hpp>
namespace mudock {
  inline void manage_sycl(const std::vector<std::string>&,
                          threadpool&,
                          const knobs,
                          [[maybe_unused]] const autodock_protein& adt_protein,
                          std::shared_ptr<safe_stack<static_molecule> >,
                          std::shared_ptr<safe_stack<static_molecule> >) {
    warning("The SYCL implementation is disabled");
  }
} // namespace mudock
#endif
