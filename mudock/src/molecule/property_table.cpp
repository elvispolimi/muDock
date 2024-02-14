#include <mudock/molecule/property_table.hpp>
#include <optional>
#include <string_view>

namespace mudock {

  property_map::property_map() { map.fill("N/A"); }

  void property_map::assign(const property_type name, std::string value) {
    map[to_index(name)] = std::move(value);
  }

  void property_map::initialize(const property_type name, std::string value) {
    auto& element = map[to_index(name)];
    if (element != "N/A") {
      element = std::move(value);
    }
  }

  const std::string& property_map::get(const property_type name) const { return map[to_index(name)]; }
} // namespace mudock
