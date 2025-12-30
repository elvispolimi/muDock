#pragma once

#include <cstddef>

namespace mudock {

  struct queue {
    queue(const int _id): id(_id) {};
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

    virtual void synchronize() = 0;

    int get_id() { return id; };

  protected:
    const int id;

    // virtual void launch_kernel_impl(void*, void*[] = nullptr, const std::string_view = {}, const int = 0) = 0;
  };
} // namespace mudock
