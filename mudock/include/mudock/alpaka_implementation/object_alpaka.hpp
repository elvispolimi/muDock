#pragma once

#include <alpaka/alpaka.hpp>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <mudock/alpaka_implementation/queue_alpaka.hpp>
#include <mudock/compute/object.hpp>
#include <optional>
#include <vector>

namespace mudock {
  template<typename T>
  struct object<T, queue_alpaka> {
  private:
    using dim         = queue_alpaka::dim;
    using idx         = queue_alpaka::idx;
    using dev_acc     = queue_alpaka::dev_acc;
    using buffer_type = alpaka::Buf<dev_acc, T, dim, idx>;
    using byte_type   = std::uint8_t;

    T* ptr                 = nullptr;
    std::size_t size       = 0;
    std::size_t alloc_size = 0;
    std::shared_ptr<queue_alpaka> q;
    std::optional<buffer_type> buffer;

    static auto extent(const std::size_t num_elements) {
      return alpaka::Vec<dim, idx>{static_cast<idx>(num_elements)};
    }

  public:
    object(std::shared_ptr<queue_alpaka> _q): q(_q) {}
    object(const object&)            = delete;
    object(object&&)                 = delete;
    object& operator=(const object&) = delete;
    object& operator=(object&&)      = delete;
    ~object() {
      if (buffer) {
        q->synchronize();
      }
    }

    void change_queue(std::shared_ptr<queue_alpaka> _q) { q = _q; }

    void alloc(const std::size_t num_elements) {
      if (num_elements > alloc_size) {
        q->synchronize();
        buffer.emplace(alpaka::allocBuf<T, idx>(q->native_device(), extent(num_elements)));
        ptr        = alpaka::getPtrNative(*buffer);
        alloc_size = num_elements;
      }
      size = num_elements;
    }

    void set_to_value(const int value) {
      const auto bytes = size * sizeof(T);
      auto host        = std::vector<byte_type>(bytes, static_cast<byte_type>(value));
      auto host_view   = alpaka::createView(alpaka::getDevByIdx(alpaka::PlatformCpu{}, 0),
                                          host.data(),
                                          alpaka::Vec<dim, idx>{static_cast<idx>(bytes)});
      auto device_view = alpaka::createView(q->native_device(),
                                            reinterpret_cast<byte_type*>(ptr),
                                            alpaka::Vec<dim, idx>{static_cast<idx>(bytes)});
      alpaka::memcpy(q->native_queue(), device_view, host_view);
    }

    void copy_host2device(const T* host, const std::size_t copy_size = 0) {
      const auto n    = copy_size ? copy_size : size;
      auto host_view  = alpaka::createView(alpaka::getDevByIdx(alpaka::PlatformCpu{}, 0), host, extent(n));
      auto device_view = alpaka::createView(q->native_device(), ptr, extent(n));
      alpaka::memcpy(q->native_queue(), device_view, host_view);
    }

    void copy_device2host(T* const host) const {
      auto device_view = alpaka::createView(q->native_device(), ptr, extent(size));
      auto host_view   = alpaka::createView(alpaka::getDevByIdx(alpaka::PlatformCpu{}, 0), host, extent(size));
      alpaka::memcpy(q->native_queue(), host_view, device_view);
    }

    void copy_device2device(object<T, queue_alpaka>& other, const int copy_size = 0) {
      const auto n   = copy_size ? static_cast<std::size_t>(copy_size) : size;
      auto src_view  = alpaka::createView(q->native_device(), other.ptr, extent(n));
      auto dest_view = alpaka::createView(q->native_device(), ptr, extent(n));
      alpaka::memcpy(q->native_queue(), dest_view, src_view);
    }

    [[nodiscard]] T* dev_pointer() const { return ptr; }
    [[nodiscard]] T** dev_pointer_ref() { return &ptr; }
    [[nodiscard]] std::size_t num_elements() const { return size; }
  };
} // namespace mudock
