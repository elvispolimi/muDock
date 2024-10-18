#include <mudock/cpp_implementation/chromosome.hpp>
#include <polygeist/cuda_object.cuh>
#include <mudock/type_alias.hpp>
#include <polygeist/cuda_random.cuh>

namespace mudock {

  template<>
  cuda_object<int>::~cuda_object() {}

  template<>
  cuda_object<fp_type>::~cuda_object() {}

  template<>
  cuda_object<fp_type*>::~cuda_object() {}

  template<>
  cuda_object<XORWOWState>::~cuda_object() {}

  template<>
  cuda_object<chromosome>::~cuda_object() {}
} // namespace mudock
