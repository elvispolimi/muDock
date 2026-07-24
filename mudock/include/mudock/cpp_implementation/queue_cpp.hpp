#pragma once

#include <cassert>
#include <functional>
#include <mudock/compute/queue.hpp>
#include <mudock/likwid_utils.hpp>
#include <mudock/log.hpp>
#include <mutex>
#include <stdexcept>

namespace mudock {
  struct queue_cpp: queue {
    queue_cpp(const int _id, [[maybe_unused]] const device_type _dev_type): queue(_id, _dev_type) {
      assert(dev_type == device_type::CPU && "CPP version supports only CPUs devices");
    };

    // TODO add kernel verion on CPU which as templates arguments for args
    template<const char* region_name, class F, class... Args>
    inline void invoke_kernel(F&& f, Args&&... args) {
      static std::once_flag flag;
      static std::string full; // stable storage
      std::call_once(flag, [] {
        full = std::string("CPP ") + region_name;
        MUDOCK_CPP_MARKER_REGISTER(full.c_str());
      });
      MUDOCK_CPP_MARKER_START(full.c_str());
      std::invoke(std::forward<F>(f), std::forward<Args>(args)...);
      MUDOCK_CPP_MARKER_STOP(full.c_str());
    }

    template<class F, class... Args>
    inline void invoke_kernel(F&& f, Args&&... args) {
      std::invoke(std::forward<F>(f), std::forward<Args>(args)...);
    }

    void alloc(void**, const size_t) {}
    void free(void**) {}
    void set_to_value(void*, const size_t, const char) {}
    void copy_host2device(const void*, void*, const size_t) {}
    void copy_device2host(const void*, void*, const size_t) {}
    void copy_device2device(const void*, void*, const size_t) {}
    bool obj_required() { return false; }
    bool honors_stage_bucket_policy() const override { return false; }

    void operator()() {
      cpu_set_t cpuset;
      CPU_ZERO(&cpuset);
      CPU_SET(id, &cpuset); // Set affinity to the target CPU
      pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
      mudock::info("Worker CPP on duty! Set affinity to core ", id);
    }

    void synchronize() {};
  };
} // namespace mudock
