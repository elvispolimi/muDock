#pragma once

#include <mudock/knobs.hpp>
#ifdef MUDOCK_USE_SYCL
  #include <mudock/sycl_implementation/sycl_batch_sizer.hpp>
  #include <mudock/sycl_implementation/sycl_manager.hpp>
  #include <mudock/sycl_implementation/sycl_object.hpp>
  #include <mudock/sycl_implementation/sycl_wrapper.hpp>
  #include <mudock/sycl_implementation/virtual_screen.hpp>
#else
  #include <mudock/log.hpp>
namespace mudock {
  inline void manage_sycl(std::string_view,
                          threadpool&,
                          const knobs,
                          [[maybe_unused]] std::shared_ptr<const grid_atom_mapper> grid_atom_maps,
                          [[maybe_unused]] std::shared_ptr<const grid_map> electro_map,
                          [[maybe_unused]] std::shared_ptr<const grid_map> desolv_map,
                          std::shared_ptr<safe_stack<static_molecule> >,
                          std::shared_ptr<safe_stack<static_molecule> >) {
    warning("The SYCL implementation is disabled");
  }
} // namespace mudock
#endif
