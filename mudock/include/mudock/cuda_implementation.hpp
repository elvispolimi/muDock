#pragma once

#ifdef MUDOCK_ENABLE_CUDA
  #include <mudock/cuda_implementation/cuda_batch_sizer.cuh>
  #include <mudock/cuda_implementation/cuda_check_error_macro.cuh>
  #include <mudock/cuda_implementation/cuda_manager.hpp>
  #include <mudock/cuda_implementation/cuda_object.cuh>
  #include <mudock/cuda_implementation/cuda_wrapper.cuh>
  #include <mudock/cuda_implementation/virtual_screen.cuh>
#else
  #include <mudock/log.hpp>
namespace mudock {
  inline void manage_cuda(std::string_view,
                          threadpool&,
                          std::shared_ptr<dynamic_molecule>,
                          std::shared_ptr<safe_stack<static_molecule> >,
                          std::shared_ptr<safe_stack<static_molecule> >) {
    warning("The CUDA implementation is disabled");
  }
} // namespace mudock
#endif
