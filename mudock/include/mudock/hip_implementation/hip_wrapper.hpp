#pragma once

#include <cassert>
#include <mudock/hip_implementation/hip_check_error_macro.hpp>
#include <mudock/hip_implementation/hip_object.hpp>

namespace mudock {

  template<template<class...> class container_type, typename T, class... args>
  class hip_wrapper: hip_object<T> {
  public:
    hip_wrapper(const hipStream_t &_stream): hip_object<T>::hip_object(_stream){};
    hip_wrapper(const hip_wrapper &)            = delete;
    hip_wrapper(hip_wrapper &&)                 = default;
    hip_wrapper &operator=(const hip_wrapper &) = delete;
    hip_wrapper &operator=(hip_wrapper &&)      = default;

    // this is the container host side
    container_type<T, args...> host;

    inline void alloc(const std::size_t num_elements) {
      host.resize(num_elements);
      hip_object<T>::alloc(num_elements);
    };
    // Set ONLY const byte value
    inline void alloc(const std::size_t num_elements, const int value) {
      host.resize(num_elements, value);
      hip_object<T>::alloc(host.size());
      hip_object<T>::set_to_value(value);
    };

    inline void copy_host2device() {
      hip_object<T>::alloc(host.size());
      hip_object<T>::copy_host2device(host.data());
    };
    inline void copy_device2host() {
      host.resize(hip_object<T>::num_elements());
      hip_object<T>::copy_device2host(host.data());
    };

    [[nodiscard]] inline auto dev_pointer() const { return hip_object<T>::dev_pointer(); }
    [[nodiscard]] inline auto host_pointer() const { return host.data(); }
    [[nodiscard]] inline auto host_pointer() { return host.data(); }
    [[nodiscard]] inline auto num_elements() const {
      assert(hip_object<T>::num_elements() == host.size());
      return hip_object<T>::num_elements();
    };
  };

} // namespace mudock
