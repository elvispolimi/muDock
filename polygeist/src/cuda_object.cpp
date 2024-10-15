#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/cuda_implementation/cuda_object.cuh>
#include <mudock/type_alias.hpp>
#include <polygeist/cuda_random.cuh>

namespace mudock {

  template<>
  cuda_object<int>::~cuda_object() noexcept(false) {}

  template<>
  cuda_object<fp_type>::~cuda_object() noexcept(false) {}

  template<>
  cuda_object<fp_type*>::~cuda_object() noexcept(false) {}

  template<>
  cuda_object<XORWOWState>::~cuda_object() noexcept(false) {}

  template<>
  cuda_object<chromosome>::~cuda_object() noexcept(false) {}
} // namespace mudock
