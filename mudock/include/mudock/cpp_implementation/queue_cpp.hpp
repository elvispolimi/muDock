#pragma once

#include <mudock/compute/queue.hpp>
#include <mudock/log.hpp>

namespace mudock {
  struct queue_cpp: queue {
    queue_cpp(const int _id): queue(_id) {};
    void launch_kernel(void* f, const int, void*[]) {
      auto fp = reinterpret_cast<void (*)()>(f);
      fp(); // calls the function
    }

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
