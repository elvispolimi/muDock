#pragma once

#include <memory>
#include <mudock/devices/context.hpp>
#include <mudock/devices/kernel_types.hpp>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>
#include <span>

namespace mudock {
  // The scratchpad is used to abstract the memory management on the device and host
  // Each language manage device memory differently
  // Thus, each worker will have a different scratchpad, which depends on the
  // language and device(though device_context) types

  // The scratchpad interfaces define which functionalities should each language expose to manager
  // interactions between host and device code 
  template<kernel_type l_t>
  class kernel_scratchpad_interface {
  public:
    // The scratchpad requires the device context in which move in and out the data
    kernel_scratchpad_interface(const device_context<l_t>& d_c): d_context(d_c){};
    virtual ~kernel_scratchpad_interface() {}

    // A scratchpad cannot be copied or moved
    kernel_scratchpad_interface(const kernel_scratchpad_interface&)            = delete;
    kernel_scratchpad_interface(kernel_scratchpad_interface&&)                 = delete;
    kernel_scratchpad_interface& operator=(const kernel_scratchpad_interface&) = delete;
    kernel_scratchpad_interface& operator=(kernel_scratchpad_interface&&)      = delete;

    // Given a span of input molecule on the host, copy ON the device memory pointed by device_context
    virtual void copy_in(const std::span<std::unique_ptr<static_molecule>> i_molecules) = 0;
    // Given a span of output molecule on the host, copy FROM the device memory pointed by device_context
    virtual void copy_out(std::span<std::unique_ptr<static_molecule>> o_molecules) = 0;

  protected:
    // Is the device's context holding this scratchpad's memory
    const device_context<l_t>& d_context;
  };

  // The scratchpad implementation class is used to provide full class specialization
  // NOTE: this enables a useful mechanism of implementation selection based on the kernel type
  template<kernel_type l_t>
  class kernel_scratchpad_impl {};
} // namespace mudock
