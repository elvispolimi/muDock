#pragma once

#include <mudock/chem/autodock_protein.hpp>
#include <mudock/sycl_implementation/sycl_texture.hpp>
#include <mudock/sycl_implementation/sycl_wrapper.hpp>
#include <sycl/sycl.hpp>

namespace mudock {
  struct device {
    const point<fp_type, 3> center, minimum, maximum;

    const autodock_protein& adt_protein;

    device(const sycl::device& dev, const autodock_protein& adt_protein);

    sycl::queue get_queue() const;
    int get_sub_group_size() const;

    const fp_type* get_tex_dev_pointer() const { return atom_tex.tex.dev_pointer(); };

  private:
    // Device ID
    const sycl::device dev;
    sycl::queue queue;
    syclTexture_wrapper atom_tex;
  };
} // namespace mudock
