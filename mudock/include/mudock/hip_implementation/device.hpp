#pragma once

#include <hip/hip_runtime.h>
#include <mudock/chem/autodock_protein.hpp>
#include <mudock/grid.hpp>
#include <mudock/hip_implementation/hip_check_error_macro.hpp>
#include <mudock/hip_implementation/hip_texture.hpp>
#include <mudock/hip_implementation/hip_wrapper.hpp>

namespace mudock {
  struct hipStream_wrapper {
    hipStream_wrapper(const int id) {
      MUDOCK_CHECK(hipSetDevice(static_cast<int>(id)));
      MUDOCK_CHECK(hipStreamCreate(&stream););
    };
    ~hipStream_wrapper() noexcept(false) { MUDOCK_CHECK(hipStreamDestroy(stream)) };
    hipStream_t& operator()() { return stream; };

  private:
    hipStream_t stream;
  };

  struct device {
    // Device ID
    const std::size_t id;
    // Grid Maps
    const point<fp_type, 3> center_maps;
    // FIXME remove from here and add to get_stream() moeve to private

    const autodock_protein& adt_protein;

    device(const std::size_t gpu_id, const autodock_protein& adt_protein);

    hipStream_wrapper get_stream() const;

    int get_wavefront() const { return wavefront_size; };
    const hipTexture_wrapper* get_tex_dev_pointer() const { return atom_tex.dev_pointer(); };

  private:
    int wavefront_size{0};
    hipStream_wrapper stream;
    hip_wrapper<std::vector, hipTexture_wrapper> atom_tex;
  };
} // namespace mudock
