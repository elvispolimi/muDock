#pragma once

#ifdef MUDOCK_ENABLE_CUDA
  #include <mudock/cuda_kernels/cuda_batch_sizer.cuh>
  #include <mudock/cuda_kernels/cuda_check_error_macro.cuh>
  #include <mudock/cuda_kernels/cuda_object.cuh>
  #include <mudock/cuda_kernels/cuda_wrapper.cuh>
  #include <mudock/cuda_kernels/virtual_screen.cuh>
#endif
