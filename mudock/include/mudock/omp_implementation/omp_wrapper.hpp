#pragma once

#include <cassert>
#include <mudock/omp_implementation/omp_object.hpp>

namespace mudock {

  template<template<class...> class container_type, typename T, class... args>
  class omp_wrapper: omp_object<T> {
  public:
    omp_wrapper()                               = default;
    omp_wrapper(const omp_wrapper &)            = delete;
    omp_wrapper(omp_wrapper &&)                 = default;
    omp_wrapper &operator=(const omp_wrapper &) = delete;
    omp_wrapper &operator=(omp_wrapper &&)      = default;

    // this is the container host side
    container_type<T, args...> host;

    inline void alloc(const std::size_t num_elements) {
      host.resize(num_elements);
      omp_object<T>::alloc(num_elements);
    };
    // CUDA set const byte value
    inline void alloc(const std::size_t num_elements, const int value) {
      host.resize(num_elements, value);
      omp_object<T>::alloc(host.size());
      omp_object<T>::set_to_value(value);
    };

    inline void copy_host2device() {
      omp_object<T>::alloc(host.size());
      omp_object<T>::copy_host2device(host.data());
    };
    inline void copy_device2host() {
      host.resize(omp_object<T>::num_elements());
      omp_object<T>::copy_device2host(host.data());
    };

    [[nodiscard]] inline auto dev_pointer() const { return omp_object<T>::dev_pointer(); }
    [[nodiscard]] inline auto host_pointer() const { return host.data(); }
    [[nodiscard]] inline auto host_pointer() { return host.data(); }
    [[nodiscard]] inline auto num_elements() const {
      assert(omp_object<T>::num_elements() == host.size());
      return omp_object<T>::num_elements();
    };
  };

} // namespace mudock
