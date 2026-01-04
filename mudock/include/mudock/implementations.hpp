#pragma once

#include <mudock/cpp_implementation/cpp_implementation.hpp>
#include <mudock/cuda_implementation/cuda_implementation.hpp>
#include <mudock/gh_implementation/gh_implementation.hpp>
#include <mudock/hip_implementation/hip_implementation.hpp>
#include <mudock/implementation_types.hpp>
#include <mudock/xsimd_implementation/xsimd_implementation.hpp>

namespace mudock {

  struct device_type_desc {
    static constexpr auto cpu_token = "CPU";
    static constexpr auto gpu_token = "GPU";
  };

  enum class device_type { CPU = 0, GPU };

  template<implementation_type impl_t, typename queue_t>
    requires std::derived_from<queue_t, queue>
  struct kernel_type_traits_impl {
    using type = queue_t;
  };

  template<implementation_type impl_t>
  struct kernel_type_traits {};

  template<>
  struct kernel_type_traits<implementation_type::CPP> {
    using type = kernel_type_traits_impl<implementation_type::CPP, queue_cpp>::type;
  };
#ifdef MUDOCK_USE_GH
  template<>
  struct kernel_type_traits<implementation_type::GH> {
    using type = kernel_type_traits_impl<implementation_type::GH, queue_gh>::type;
  };
#endif
#ifdef MUDOCK_USE_XSIMD
  template<>
  struct kernel_type_traits<implementation_type::XSIMD> {
    using type = kernel_type_traits_impl<implementation_type::XSIMD, queue_xsimd>::type;
  };
#endif
#ifdef MUDOCK_USE_CUDA
  template<>
  struct kernel_type_traits<implementation_type::CUDA> {
    using type = kernel_type_traits_impl<implementation_type::CUDA, queue_cuda>::type;
  };
#endif
#ifdef MUDOCK_USE_HIP
  template<>
  struct kernel_type_traits<implementation_type::HIP> {
    using type = kernel_type_traits_impl<implementation_type::HIP, queue_hip>::type;
  };
#endif

  static constexpr int num_cpu_kernel_type() {
    return 1
#ifdef MUDOCK_USE_GH
           + 1
#endif
#ifdef MUDOCK_USE_XSIMD
           + 1
#endif
        ;
  }
  static constexpr std::array<implementation_type, num_cpu_kernel_type()> cpu_kernel_type{
      implementation_type::CPP
#ifdef MUDOCK_USE_GH
      ,
      implementation_type::GH
#endif
#ifdef MUDOCK_USE_XSIMD
      ,
      implementation_type::XSIMD
#endif
  };

  static constexpr int num_gpu_kernel_type() {
    return 0
#ifdef MUDOCK_USE_CUDA
           + 1
#endif
#ifdef MUDOCK_USE_HIP
           + 1
#endif
        ;
  }
  static constexpr std::array<implementation_type, num_gpu_kernel_type()> gpu_kernel_type{
#ifdef MUDOCK_USE_CUDA
      implementation_type::CUDA,
#endif
#ifdef MUDOCK_USE_HIP
      implementation_type::HIP,
#endif
  };

  inline device_type get_device_type(const std::string_view& impl) {
    if (impl == device_type_desc::cpu_token)
      return device_type::CPU;
    if (impl == device_type_desc::gpu_token)
      return device_type::GPU;
    throw std::runtime_error("Requested device not available");
  };
} // namespace mudock
