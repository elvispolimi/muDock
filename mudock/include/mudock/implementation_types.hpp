#pragma once

#include <string_view>
namespace mudock {
  struct implementation_type_desc {
    static constexpr auto cpp_token = "CPP";
#ifdef MUDOCK_USE_GH
    static constexpr auto gh_token = "GH";
#endif
#ifdef MUDOCK_USE_XSIMD
    static constexpr auto xsimd_token = "XSIMD";
#endif
#ifdef MUDOCK_USE_CUDA
    static constexpr auto cuda_token = "CUDA";
#endif
  };

  enum class implementation_type {
    CPP = 0
#ifdef MUDOCK_USE_GH
    ,
    GH
#endif
#ifdef MUDOCK_USE_XSIMD
    ,
    XSIMD
#endif
#ifdef MUDOCK_USE_CUDA
    ,
    CUDA
#endif
  };

  inline implementation_type get_impl_type(const std::string_view& impl) {
    if (impl == implementation_type_desc::cpp_token)
      return implementation_type::CPP;
#ifdef MUDOCK_USE_GH
    if (impl == implementation_type_desc::gh_token)
      return implementation_type::GH;
#endif
#ifdef MUDOCK_USE_XSIMD
    if (impl == implementation_type_desc::xsimd_token)
      return implementation_type::XSIMD;
#endif
#ifdef MUDOCK_USE_CUDA
    if (impl == implementation_type_desc::cuda_token)
      return implementation_type::CUDA;
#endif
    throw std::runtime_error("Requested implementation not available");
  };

} // namespace mudock
