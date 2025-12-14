#pragma once

#include <functional>
#include <memory>

namespace mudock {

  struct queue {
    queue(const int _id): id(_id) {};
    virtual void launch_kernel(std::function<void()>)                                                 = 0;
    virtual void alloc(void* ptr, const size_t num_bytes)                                             = 0;
    virtual void free(void* ptr)                                                                      = 0;
    virtual void set_to_value(void* ptr, const size_t num_bytes, const int value)                     = 0;
    virtual void copy_host2device(const void* const host, const void* device, const size_t num_bytes) = 0;
    virtual void copy_device2host(const void* const device, const void* host, const size_t num_bytes) = 0;
    virtual void
        copy_device2device(const void* const device1, const void* device2, const size_t num_bytes) = 0;

    virtual void operator()() = 0;

    virtual bool obj_required() = 0;

  protected:
    const int id;
  };
} // namespace mudock
