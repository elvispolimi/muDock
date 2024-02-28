#pragma once

#include <array>
#include <cassert>
#include <mudock/molecule/properties.hpp>
#include <string>

namespace mudock {

  class property_map {
    std::array<std::string, property_type_size()> map;

    inline constexpr auto to_index(const property_type type) const {
      assert(static_cast<std::size_t>(type) < property_type_size());
      return static_cast<std::size_t>(type);
    }

  public:
    property_map();                                               // fill all the properties with "N/A"
    void assign(const property_type name, std::string value);     // overwrite previous values
    void initialize(const property_type name, std::string value); // overwrite "N/A" values only
    [[nodiscard]] const std::string& get(const property_type name) const; // get the value
  };

} // namespace mudock
