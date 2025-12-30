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
    queue_cpp(const int _id): queue(_id) {};
    // void launch_kernel(void* f, void* args[], const int s) {
    //   assert(args != nullptr && "Wrong args to invoke the kernels");
    //   assert(s == 0 && "Invoked CPU kernels with dimensions different than 0");
    //   auto fp = reinterpret_cast<cpu_invoker_t>(f);
    //   fp(args[0]); // calls the function
    // }
    // void launch_kernel(void* f, void*[], const std::string_view region_name, const int) {
    //   assert(!region_name.empty() && "Requeste kernel launch with empty region name");
    //   auto fp = reinterpret_cast<void (*)()>(f);
    //   if (!region_name.empty()) {
    //     LIKWID_MARKER_REGISTER(region_name);
    //     LIKWID_MARKER_START(region_name);
    //   }
    //   fp(); // calls the function
    //   if (!region_name.empty())
    //     LIKWID_MARKER_STOP(region_name);
    // }
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

    // void launch_kernel(void*, void*[], const std::string_view, const int) {
    //   throw std::runtime_error("Launch kernel function not defined for CPU code");
    // }
    //
    // void launch_kernel(void*, void*[], const int) {
    //   throw std::runtime_error("Launch kernel function not defined for CPU code");
    // }

    void alloc(void**, const size_t) {
      // throw std::runtime_error("Invoked alloc on CPP queue");
    }
    void free(void**) {}
    void set_to_value(void*, const size_t, const char) {
      // throw std::runtime_error("Invoked set_to_value on CPP queue");
    }
    void copy_host2device(const void*, void*, const size_t) {
      // throw std::runtime_error("Invoked copy_host2device on CPP queue");
    }
    void copy_device2host(const void*, void*, const size_t) {
      // throw std::runtime_error("Invoked device2host on CPP queue");
    }
    void copy_device2device(const void*, void*, const size_t) {
      // throw std::runtime_error("Invoked device2device mem on CPP queue");
    }
    bool obj_required() { return false; }

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
