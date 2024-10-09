#pragma once

#include <cassert>
#include <cuda_runtime.h>
#include <mudock/cuda_implementation/cuda_object.cuh>

namespace mudock {

  template<template<class...> class container_type, typename T, class... args>
  class cuda_wrapper: cuda_object<T> {
  public:
    cuda_wrapper()                                = default;
    cuda_wrapper(const cuda_wrapper &)            = delete;
    cuda_wrapper(cuda_wrapper &&)                 = default;
    cuda_wrapper &operator=(const cuda_wrapper &) = delete;
    cuda_wrapper &operator=(cuda_wrapper &&)      = default;

    // this is the container host side
    container_type<T, args...> host;

    inline void alloc(const std::size_t num_elements) {
      host.resize(num_elements);
      cuda_object<T>::alloc(num_elements);
    };
    // CUDA set const byte value
    inline void alloc(const std::size_t num_elements, const int value) {
      host.resize(num_elements, value);
      cuda_object<T>::alloc(host.size());
      cuda_object<T>::set_to_value(value);
    };

    inline void copy_host2device() {
      cuda_object<T>::alloc(host.size());
      cuda_object<T>::copy_host2device(host.data());
    };
    inline void copy_device2host() {
      host.resize(cuda_object<T>::num_elements());
      cuda_object<T>::copy_device2host(host.data());
    };

    [[nodiscard]] inline auto dev_pointer() const { return cuda_object<T>::dev_pointer(); }
    [[nodiscard]] inline auto host_pointer() const { return host.data(); }
    [[nodiscard]] inline auto host_pointer() { return host.data(); }
    [[nodiscard]] inline auto num_elements() const {
      assert(cuda_object<T>::num_elements() == host.size());
      return cuda_object<T>::num_elements();
    };
  };

} // namespace mudock
