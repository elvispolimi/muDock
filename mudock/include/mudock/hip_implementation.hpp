#pragma once

#include <mudock/knobs.hpp>
#ifdef MUDOCK_USE_HIP
  #include <mudock/hip_implementation/hip_batch_sizer.hpp>
  #include <mudock/hip_implementation/hip_check_error_macro.hpp>
  #include <mudock/hip_implementation/hip_manager.hpp>
  #include <mudock/hip_implementation/virtual_screen.hpp>
#else
  #include <mudock/log.hpp>
namespace mudock {
  inline void manage_hip(std::string_view,
                         threadpool&,
                         const knobs,
                         [[maybe_unused]] std::shared_ptr<const grid_atom_mapper> grid_atom_maps,
                         [[maybe_unused]] std::shared_ptr<const grid_map> electro_map,
                         [[maybe_unused]] std::shared_ptr<const grid_map> desolv_map,
                         std::shared_ptr<safe_stack<static_molecule> >,
                         std::shared_ptr<safe_stack<static_molecule> >) {
    warning("The HIP implementation is disabled");
  }
} // namespace mudock
#endif
