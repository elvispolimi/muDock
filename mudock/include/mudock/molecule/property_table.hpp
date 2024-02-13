#pragma once

#include <array>
#include <mudock/molecule/properties.hpp>
#include <string>

namespace mudock {

  class property_map {
    std::array<std::string, 1> map;

    inline constexpr auto to_index(const property_type type) const { return static_cast<std::size_t>(type); }

  public:
    void assign(const property_type name, std::string value);
    void initialize(const property_type name, std::string value);
    const std::string& get(const property_type name) const;
  };

} // namespace mudock
