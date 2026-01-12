#pragma once

#include <mudock/sycl_implementation/queue_sycl.hpp>
#include <mudock/sycl_implementation/sycl_utils.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {
  struct sycl_texture_devices {
    fp_type* tex_dev = nullptr;
    const int dev_id = 0;
    queue_sycl q;

    sycl_texture_devices(const int id,
                         const device_type dev_type,
                         const int xyz,
                         const int num_tex,
                         const fp_type* src)
        : dev_id(id), q(id, dev_type) {
      const int num_elements  = xyz * num_tex;
      const std::size_t bytes = num_elements * sizeof(fp_type);

      q.alloc(reinterpret_cast<void**>(&tex_dev), bytes);
      q.copy_host2device(src, tex_dev, bytes);

      q.synchronize();
    };

    // no copy
    sycl_texture_devices(const sycl_texture_devices&)            = delete;
    sycl_texture_devices& operator=(const sycl_texture_devices&) = delete;

    // no move (matches your HIP version)
    sycl_texture_devices(sycl_texture_devices&&)            = delete;
    sycl_texture_devices& operator=(sycl_texture_devices&&) = delete;

    ~sycl_texture_devices() noexcept {
      if (tex_dev) {
        try {
          q.free(reinterpret_cast<void**>(&tex_dev));
        } catch (const sycl::exception& e) {}
        tex_dev = nullptr;
      }
    }
  };
} // namespace mudock
