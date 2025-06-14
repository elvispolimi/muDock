#pragma once

#include <mudock/chem/autodock_protein.hpp>
#include <mudock/compute/safe_stack.hpp>
#include <mudock/compute/threadpool.hpp>
#include <mudock/knobs.hpp>
#ifdef MUDOCK_USE_HIP
  #include <mudock/hip_implementation/hip_batch_sizer.hpp>
  #include <mudock/hip_implementation/hip_check_error_macro.hpp>
  #include <mudock/hip_implementation/hip_manager.hpp>
  #include <mudock/hip_implementation/virtual_screen.hpp>
#else
  #include <mudock/log.hpp>
namespace mudock {
  inline void manage_hip(const std::vector<std::string>&,
                         threadpool&,
                         const knobs,
                         [[maybe_unused]] const autodock_protein& adt_protein,
                         std::shared_ptr<safe_stack<static_molecule> >,
                         std::shared_ptr<safe_stack<static_molecule> >) {
    warning("The HIP implementation is disabled");
  }
} // namespace mudock
#endif
