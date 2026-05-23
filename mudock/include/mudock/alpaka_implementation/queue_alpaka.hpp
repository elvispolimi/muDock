#pragma once

#include <cstddef>
#include <memory>
#include <mudock/alpaka_implementation/alpaka_types.hpp>
#include <mudock/compute/queue.hpp>
#include <mudock/grid/mdindex.hpp>

namespace mudock {
  struct queue_alpaka: queue {
    using dim       = alpaka_backend::dim;
    using idx       = alpaka_backend::idx;
    using acc       = alpaka_backend::acc;
    using dev_acc   = alpaka_backend::dev_acc;
    using queue_acc = alpaka_backend::queue_acc;

    queue_alpaka(const int _id, const device_type _dev_type);
    ~queue_alpaka();

    queue_alpaka(const queue_alpaka&)            = delete;
    queue_alpaka& operator=(const queue_alpaka&) = delete;

    queue_alpaka(queue_alpaka&&) noexcept;
    queue_alpaka& operator=(queue_alpaka&&) noexcept = delete;

    template<class F, class... Args>
    inline void invoke_kernel(const index3D gridDim, const index3D blockDim, Args&&... args);

    template<class F, class... Args>
    inline void invoke_kernel(const int gridDim, const int blockDim, Args&&... args);

    static constexpr int block_threads(const int requested_threads);

    queue_acc& native_queue();
    const dev_acc& native_device() const;

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

  private:
    struct impl;
    std::unique_ptr<impl> impl_;
  };
} // namespace mudock
