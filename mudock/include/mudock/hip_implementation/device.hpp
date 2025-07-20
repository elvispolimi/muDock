#pragma once

#include <hip/hip_runtime.h>
#include <mudock/chem/autodock_protein.hpp>
#include <mudock/grid.hpp>
#include <mudock/hip_implementation/hip_check_error_macro.hpp>
#include <mudock/hip_implementation/hip_texture.hpp>
#include <mudock/hip_implementation/hip_wrapper.hpp>

namespace mudock {
  struct device {
    // Device ID
    const std::size_t id;
    // Grid Maps
    const point<fp_type, 3> center_maps;
    // FIXME remove from here and add to get_stream() moeve to private
    struct hipStream_wrapper {
      hipStream_t stream;
      hipStream_wrapper(const hipStream_t stream): stream(stream) {};
      ~hipStream_wrapper() noexcept(false) { MUDOCK_CHECK(hipStreamDestroy(stream)) };
      hipStream_t& operator()() { return stream; };
    } stream;

    hip_wrapper<std::vector, hipTexture_wrapper> atom_tex;
    const autodock_protein& adt_protein;

    device(const std::size_t gpu_id, const autodock_protein& adt_protein);

    hipStream_t get_stream() const;

    int get_wavefront() const { return wavefront_size; };

  private:
    int wavefront_size{0};
  };
} // namespace mudock
