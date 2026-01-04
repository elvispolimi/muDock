#pragma once

#include <hip/hip_runtime.h>
#include <mudock/hip_implementation/hip_utils.hpp>
#include <mudock/hip_implementation/queue_hip.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {
  struct hip_texture_devices {
    fp_type* tex_dev;
    const int dev_id = 0;

    hip_texture_devices(const int id, const int xyz, const int num_tex, const fp_type* src): dev_id(id) {
      const int num_elements = xyz * num_tex;

      MUDOCK_CHECK(hipSetDevice(dev_id));
      MUDOCK_CHECK(hipMalloc(&tex_dev, num_elements * sizeof(fp_type)));
      MUDOCK_CHECK(hipMemcpy(tex_dev, src, num_elements * sizeof(fp_type), hipMemcpyHostToDevice));
      MUDOCK_CHECK(hipDeviceSynchronize());
    };

    // no copy
    hip_texture_devices(const hip_texture_devices&)            = delete;
    hip_texture_devices& operator=(const hip_texture_devices&) = delete;

    // no move
    hip_texture_devices(hip_texture_devices&& other)            = delete;
    hip_texture_devices& operator=(hip_texture_devices&& other) = delete;

    ~hip_texture_devices() noexcept {
      if (tex_dev) {
        (void) hipSetDevice(dev_id);
        (void) hipFree(tex_dev);
      }
    }
  };
} // namespace mudock
