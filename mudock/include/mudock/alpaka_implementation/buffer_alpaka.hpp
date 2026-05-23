#pragma once

#include <alpaka/alpaka.hpp>

#include <cstdint>
#include <memory>
#include <mudock/alpaka_implementation/queue_alpaka.hpp>
#include <mudock/compute/buffer.hpp>
#include <optional>
#include <vector>

namespace mudock {
  // Specialization for the generic buffer_impl declared in compute/buffer.hpp.
  template<template<class...> class container_type, typename T, class... args>
  struct buffer_impl<container_type, T, queue_alpaka, args...> {
  private:
    using dim         = queue_alpaka::dim;
    using idx         = queue_alpaka::idx;
    using dev_acc     = queue_alpaka::dev_acc;
    using buffer_type = alpaka::Buf<dev_acc, T, dim, idx>;
    using byte_type   = std::uint8_t;

    container_type<T, args...> host;
    std::shared_ptr<queue_alpaka> q;
    std::optional<buffer_type> buffer;
    T* ptr                 = nullptr;
    std::size_t size       = 0;
    std::size_t alloc_size = 0;
    bool valid             = false;

    static auto extent(const std::size_t num_elements) {
      return alpaka::Vec<dim, idx>{static_cast<idx>(num_elements)};
    }

    static auto host_device() { return alpaka::getDevByIdx(alpaka::PlatformCpu{}, 0); }

    void alloc_device(const std::size_t num_elements) {
      if (num_elements > alloc_size) {
        q->synchronize();
        buffer.emplace(alpaka::allocBuf<T, idx>(q->native_device(), extent(num_elements)));
        ptr        = alpaka::getPtrNative(*buffer);
        alloc_size = num_elements;
      }
      size = num_elements;
    }

  public:
    buffer_impl(std::shared_ptr<queue_alpaka> _queue, const std::size_t num_elements = 0): q(_queue) {
      if (num_elements)
        this->alloc(num_elements);
    }
    buffer_impl(const buffer_impl& other) = delete;

    buffer_impl& operator=(const buffer_impl& other) = delete;

    buffer_impl(buffer_impl&&) = default;

    buffer_impl& operator=(buffer_impl&&) = default;

    bool is_valid() { return valid; }
    void set_valid() { valid = true; }
    void set_not_valid() { valid = false; }

    ~buffer_impl() {
      if (buffer) {
        q->synchronize();
      }
    }
    [[nodiscard]] inline auto host_pointer() const { return host.data(); }
    [[nodiscard]] inline auto host_pointer() { return host.data(); }
    [[nodiscard]] inline auto num_elements() const { return host.size(); };

    auto* operator()() { return host.data(); }

    inline void alloc(const std::size_t num_elements) {
      host.resize(num_elements);
      alloc_device(num_elements);
    };
    inline void alloc(const std::size_t num_elements, const int value) {
      host.resize(num_elements, value);
      alloc_device(num_elements);

      const auto bytes = size * sizeof(T);
      auto bytes_host  = std::vector<byte_type>(bytes, static_cast<byte_type>(value));
      auto host_view   = alpaka::createView(host_device(),
                                          bytes_host.data(),
                                          alpaka::Vec<dim, idx>{static_cast<idx>(bytes)});
      auto device_view = alpaka::createView(q->native_device(),
                                            reinterpret_cast<byte_type*>(ptr),
                                            alpaka::Vec<dim, idx>{static_cast<idx>(bytes)});
      alpaka::memcpy(q->native_queue(), device_view, host_view);
    };

    inline void copy_host2device(const std::size_t copy_size = 0) {
      valid        = true;
      const auto n = copy_size ? copy_size : host.size();
      alloc_device(n);
      auto host_view   = alpaka::createView(host_device(), host.data(), extent(n));
      auto device_view = alpaka::createView(q->native_device(), ptr, extent(n));
      alpaka::memcpy(q->native_queue(), device_view, host_view);
    };
    inline void copy_device2host() {
      valid            = false;
      auto device_view = alpaka::createView(q->native_device(), ptr, extent(size));
      auto host_view   = alpaka::createView(host_device(), host.data(), extent(size));
      alpaka::memcpy(q->native_queue(), host_view, device_view);
    };
    inline void copy_device2device(const buffer_impl<container_type, T, queue_alpaka, args...>& other,
                                   const int n = 0) {
      valid = true;
      alloc_device(other.num_elements());
      host.resize(other.num_elements());

      const auto copy_size = n ? static_cast<std::size_t>(n) : size;
      auto src_view        = alpaka::createView(q->native_device(), other.ptr, extent(copy_size));
      auto dest_view       = alpaka::createView(q->native_device(), ptr, extent(copy_size));
      alpaka::memcpy(q->native_queue(), dest_view, src_view);
    }

    [[nodiscard]] inline auto dev_pointer() { return ptr; }

    [[nodiscard]] T** dev_pointer_ref() { return &ptr; }
  };
} // namespace mudock
