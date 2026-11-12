#pragma once

#include <memory>
#include <mudock/compute/batch_multiple.hpp>
#include <mudock/compute/queue.hpp>
#include <mudock/grid/mdindex.hpp>

namespace mudock {
  struct queue_sycl: queue {
    queue_sycl(const int _id, const device_type dev_t);
    ~queue_sycl();

    // non-copyable, but movable (optional)
    queue_sycl(const queue_sycl&)            = delete;
    queue_sycl& operator=(const queue_sycl&) = delete;

    queue_sycl(queue_sycl&&) noexcept;
    queue_sycl& operator=(queue_sycl&&) noexcept;

    // Implementation is missing, please inculde also invoke_kernel_header.hpp for TU which call this method
    template<class F, class... Args>
    inline void invoke_kernel(const index3D gridDim, const index3D blockDim, Args&&... args);

    template<class F, class... Args>
    inline void invoke_kernel(const int gridDim, const int blockDim, Args&&... args);

    void alloc(void**, const size_t) override;
    void free(void**) override;
    void set_to_value(void*, const size_t, const char) override;
    void copy_host2device(const void*, void*, const size_t) override;
    void copy_device2host(const void*, void*, const size_t) override;
    void copy_device2device(const void*, void*, const size_t) override;
    bool obj_required() override { return true; }
    bool honors_stage_bucket_policy() const override { return dev_type == device_type::GPU; }

    void operator()() override;

    void synchronize() override;

    template<class kernel_name>
    batch_multiple get_batch_multiple();

  private:
    struct impl;
    std::unique_ptr<impl> impl_;
  };

} // namespace mudock
