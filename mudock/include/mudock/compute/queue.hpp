#pragma once

#include <cstddef>
#include <atomic>
#include <memory>
#include <mudock/devices.hpp>

namespace mudock {

  struct device_memory_tracker {
    std::atomic<std::size_t> current{0};
    std::atomic<std::size_t> peak{0};

    void allocate(const std::size_t bytes) {
      const auto value = current.fetch_add(bytes, std::memory_order_relaxed) + bytes;
      auto previous    = peak.load(std::memory_order_relaxed);
      while (previous < value &&
             !peak.compare_exchange_weak(previous, value, std::memory_order_relaxed)) {}
    }

    void release(const std::size_t bytes) { current.fetch_sub(bytes, std::memory_order_relaxed); }
  };

  struct queue {
    queue(const int _id,
          const device_type d_t,
          std::shared_ptr<device_memory_tracker> tracker = {})
        : id(_id), dev_type(d_t), memory_tracker(tracker ? std::move(tracker)
                                                          : std::make_shared<device_memory_tracker>()) {};
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
    virtual void free(void**, const size_t bytes)                     = 0;
    virtual void set_to_value(void*, const size_t, const char)        = 0;
    virtual void copy_host2device(const void*, void*, const size_t)   = 0;
    virtual void copy_device2host(const void*, void*, const size_t)   = 0;
    virtual void copy_device2device(const void*, void*, const size_t) = 0;

    virtual void operator()() = 0;

    virtual bool obj_required() = 0;

    virtual bool honors_stage_bucket_policy() const = 0;

    virtual void synchronize() = 0;

    virtual std::size_t allocated_bytes() const       = 0;
    virtual std::size_t peak_allocated_bytes() const  = 0;

    std::size_t device_allocated_bytes() const { return memory_tracker->current.load(); }
    std::size_t device_peak_allocated_bytes() const { return memory_tracker->peak.load(); }
    std::shared_ptr<device_memory_tracker> get_memory_tracker() const { return memory_tracker; }

    int get_id() { return id; };

    device_type get_dev_type() { return dev_type; };

  protected:
    const int id;
    const device_type dev_type;
    std::shared_ptr<device_memory_tracker> memory_tracker;
  };
} // namespace mudock
