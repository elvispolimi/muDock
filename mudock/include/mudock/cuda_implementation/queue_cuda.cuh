#pragma once

#include <mudock/compute/queue.hpp>
#include <mudock/likwid_utils.hpp>

#define BLOCK_SIZE 32

namespace mudock {
  struct queue_cuda: queue {
    queue_cuda(const int _id);
    ~queue_cuda();

    // non-copyable, but movable (optional)
    queue_cuda(const queue_cuda&)            = delete;
    queue_cuda& operator=(const queue_cuda&) = delete;

    queue_cuda(queue_cuda&&) noexcept;
    queue_cuda& operator=(queue_cuda&&) noexcept;

    void launch_kernel(void* f, void*[], const int);

    void alloc(void**, const size_t);
    void free(void**);
    void set_to_value(void*, const size_t, const char);
    void copy_host2device(const void*, void*, const size_t);
    void copy_device2host(const void*, void*, const size_t);
    void copy_device2device(const void*, void*, const size_t);
    bool obj_required() { return true; }

    void operator()();

    void synchronize();

  private:
    struct impl;
    std::unique_ptr<impl> impl_;
  };

} // namespace mudock
