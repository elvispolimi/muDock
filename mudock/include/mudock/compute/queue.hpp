#pragma once

#include <cstddef>
#include <mudock/devices.hpp>

namespace mudock {

  struct queue {
    queue(const int _id, const device_type d_t): id(_id), dev_type(d_t) {};
    virtual ~queue() = default;

    queue(const queue&)            = delete;
    queue& operator=(const queue&) = delete;

    queue(queue&&) noexcept            = default; // OK
    queue& operator=(queue&&) noexcept = delete;  // const id prevents assignment

    // TODO would be nice to unify launch kernel interface yet again
    // virtual void launch_kernel(void*, void*[] = nullptr, const int = 0) = 0;
    // template<const char* Region>
    // void launch_kernel(void* f, void* args[] = nullptr, int batch = 0) {
    //   this->launch_kernel_impl(f, args, std::string_view{Region}, batch);
    // }
    virtual void alloc(void**, const size_t)                          = 0;
    virtual void free(void**)                                         = 0;
    virtual void set_to_value(void*, const size_t, const char)        = 0;
    virtual void copy_host2device(const void*, void*, const size_t)   = 0;
    virtual void copy_device2host(const void*, void*, const size_t)   = 0;
    virtual void copy_device2device(const void*, void*, const size_t) = 0;

    virtual void operator()() = 0;

    virtual bool obj_required() = 0;

    virtual bool honors_stage_bucket_policy() const = 0;

    virtual void synchronize() = 0;

    int get_id() { return id; };

    device_type get_dev_type() { return dev_type; };

  protected:
    const int id;
    const device_type dev_type;
  };
} // namespace mudock
