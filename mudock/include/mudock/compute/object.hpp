#pragma once

#include <cstddef>
#include <memory>
#include <mudock/compute/queue.hpp>

namespace mudock {
  template<typename T>
  struct object {
  protected:
    T* ptr                 = nullptr;
    std::size_t size       = 0;
    std::size_t alloc_size = 0;
    std::shared_ptr<queue> q;

  public:
    object(std::shared_ptr<queue> _q): q(_q) {};
    object(const object&)  = delete;
    object(object&& other) = delete;
    // TODO fix me the noexcept, change the mudock check
    ~object() noexcept(false) {
      if (ptr && q)
        q->free((void**) &ptr);
    };
    object& operator=(const object&) = delete;
    object& operator=(object&&)      = delete;

    void change_queue(std::shared_ptr<queue> _q) { q = _q; };

    void alloc(const size_t num_elements) {
      if (num_elements > alloc_size) {
        if (ptr)
          q->free((void**) &ptr);
        q->alloc((void**) &ptr, num_elements * sizeof(T));
        alloc_size = num_elements;
      }
      size = num_elements;
    };

    // It works only if we considers the first 8 bit of the value, as for major implementations it sets bytes
    void set_to_value(const int value) { q->set_to_value(ptr, size, value); }

    void copy_host2device(const T* host, const std::size_t copy_size = 0) {
      q->copy_host2device(host, ptr, (copy_size ? copy_size : size) * sizeof(T));
    };
    void copy_device2host(T* const host) const { q->copy_device2host(ptr, host, size * sizeof(T)); };
    void copy_device2device(object<T>& other, const int copy_size = -1) {
      const auto n = (copy_size ? copy_size : size);
      // TODO check if necessary
      // alloc(n);
      q->copy_device2device(ptr, other.ptr, n * sizeof(T));
    };

    [[nodiscard]] T* dev_pointer() const { return ptr; };
    [[nodiscard]] std::size_t num_elements() const { return size; };
  };
} // namespace mudock
