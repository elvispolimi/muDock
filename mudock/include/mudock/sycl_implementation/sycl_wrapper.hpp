#pragma once

#include <cassert>
#include <mudock/sycl_implementation/sycl_object.hpp>
#include <sycl/sycl.hpp>

namespace mudock {

  template<template<class...> class container_type, typename T, class... args>
  class sycl_wrapper: sycl_object<T> {
  public:
    sycl_wrapper(sycl::queue &_queue): sycl_object<T>(_queue){};
    sycl_wrapper(const sycl_wrapper &)            = delete;
    sycl_wrapper(sycl_wrapper &&)                 = default;
    sycl_wrapper &operator=(const sycl_wrapper &) = delete;
    sycl_wrapper &operator=(sycl_wrapper &&)      = default;

    // this is the container host side
    container_type<T, args...> host;

    inline void alloc(const std::size_t num_elements) {
      host.resize(num_elements);
      sycl_object<T>::alloc(num_elements);
    };
    inline void alloc(const std::size_t num_elements, const T value) {
      host.resize(num_elements, value);
      sycl_object<T>::alloc(host.size());
      sycl_object<T>::set_to_value(value);
    };

    inline void copy_host2device() {
      sycl_object<T>::alloc(host.size());
      sycl_object<T>::copy_host2device(host.data());
    };
    inline void copy_device2host() {
      host.resize(sycl_object<T>::num_elements());
      sycl_object<T>::copy_device2host(host.data());
    };

    [[nodiscard]] inline auto dev_pointer() const { return sycl_object<T>::dev_pointer(); }
    [[nodiscard]] inline auto host_pointer() const { return host.data(); }
    [[nodiscard]] inline auto host_pointer() { return host.data(); }
    [[nodiscard]] inline auto num_elements() const {
      assert(sycl_object<T>::num_elements() == host.size());
      return sycl_object<T>::num_elements();
    };
  };

} // namespace mudock
