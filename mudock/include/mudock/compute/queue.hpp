#pragma once

#include <functional>
#include <memory>

namespace mudock {

  struct queue {
    queue(const int _id): id(_id) {};
    virtual ~queue() = default;

    // non-copyable, but movable (optional)
    queue(const queue&)            = delete;
    queue& operator=(const queue&) = delete;

    queue(queue&&) noexcept;
    queue& operator=(queue&&) noexcept;

    virtual void launch_kernel(void*, const int, void*[])             = 0;
    virtual void alloc(void**, const size_t)                          = 0;
    virtual void free(void**)                                         = 0;
    virtual void set_to_value(void*, const size_t, const char)        = 0;
    virtual void copy_host2device(const void*, void*, const size_t)   = 0;
    virtual void copy_device2host(const void*, void*, const size_t)   = 0;
    virtual void copy_device2device(const void*, void*, const size_t) = 0;

    virtual void operator()() = 0;

    virtual bool obj_required() = 0;

    virtual void synchronize();

    int get_id() { return id; };

  protected:
    const int id;
  };
} // namespace mudock
