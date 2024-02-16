#pragma once

#include <memory>
#include <mudock/devices/device_contexts.hpp>
#include <mudock/devices/language_types.hpp>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>
#include <span>

namespace mudock {
  template<language_types l_t>
  class molecules_scratchpad {
  public:
    molecules_scratchpad();
    ~molecules_scratchpad();

    molecules_scratchpad(const molecules_scratchpad&)            = delete;
    molecules_scratchpad(molecules_scratchpad&&)                 = delete;
    molecules_scratchpad& operator=(const molecules_scratchpad&) = delete;
    molecules_scratchpad& operator=(molecules_scratchpad&&)      = delete;

    void copy_in(const std::span<std::unique_ptr<static_molecule>>, const device_context<l_t>&);
    void copy_out(std::span<std::unique_ptr<static_molecule>>, const device_context<l_t>&);
  };
} // namespace mudock
