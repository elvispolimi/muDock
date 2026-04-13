#pragma once

#include <cassert>
#include <mudock/compute/queue.hpp>
#include <mudock/log.hpp>

namespace mudock {
  struct queue_alpaka: queue {
    queue_alpaka(const int _id, const device_type _dev_type): queue(_id, _dev_type) {
      assert((dev_type == device_type::CPU || dev_type == device_type::GPU) &&
             "Alpaka skeleton supports only CPU or GPU device tags");
    }

    void alloc(void**, const size_t) override {}
    void free(void**) override {}
    void set_to_value(void*, const size_t, const char) override {}
    void copy_host2device(const void*, void*, const size_t) override {}
    void copy_device2host(const void*, void*, const size_t) override {}
    void copy_device2device(const void*, void*, const size_t) override {}
    bool obj_required() override { return false; }

    void operator()() override {
      mudock::info("Worker ALPAKA skeleton on duty for device ", id);
    }

    void synchronize() override {}
  };
} // namespace mudock
