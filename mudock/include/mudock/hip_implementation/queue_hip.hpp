#pragma once

#include <memory>
#include <mudock/compute/queue.hpp>
#include <mudock/grid/mdindex.hpp>

#ifdef __HIP_PLATFORM_AMD__
  #define BLOCK_SIZE 64
#else
  #define BLOCK_SIZE 32
#endif

namespace mudock {
  struct queue_hip: queue {
    queue_hip(const int _id, const device_type dev_t);
    ~queue_hip();

    // non-copyable, but movable (optional)
    queue_hip(const queue_hip&)            = delete;
    queue_hip& operator=(const queue_hip&) = delete;

    queue_hip(queue_hip&&) noexcept;
    queue_hip& operator=(queue_hip&&) noexcept;

    void launch_kernel(void*, void*[], const index3D, const index3D);
    void launch_kernel(void*, void*[], const int, const int);

    void alloc(void**, const size_t);
    void free(void**);
    void set_to_value(void*, const size_t, const char);
    void copy_host2device(const void*, void*, const size_t);
    void copy_device2host(const void*, void*, const size_t);
    void copy_device2device(const void*, void*, const size_t);
    bool obj_required() { return true; }
    bool honors_stage_bucket_policy() const override { return true; }

    void operator()();

    void synchronize();

  private:
    struct impl;
    std::unique_ptr<impl> impl_;
  };

} // namespace mudock
