#include <algorithm>
#include <mudock/chem/autodock_grid_types.hpp>
#include <stdexcept>

namespace mudock {

  autodock_grid_type parse_map_symbol(const std::string_view symbol) {
    const auto element_it = std::find_if(std::begin(MAP_DICTIONARY),
                                         std::end(MAP_DICTIONARY),
                                         [&symbol](const auto& e) { return e.name == symbol; });
    if (element_it != std::end(MAP_DICTIONARY))
      return element_it->value;
    else
      throw std::runtime_error("Missing map type");
  }

  const std::array<map_description, num_autodock_grids()> MAP_DICTIONARY = {
      {{autodock_grid_type::A, "A"},
       {autodock_grid_type::C, "C"},
       {autodock_grid_type::H, "H"},
       {autodock_grid_type::HD, "HD"},
       {autodock_grid_type::N, "N"},
       {autodock_grid_type::NA, "NA"},
       {autodock_grid_type::OA, "OA"},
       {autodock_grid_type::SA, "SA"},
       {autodock_grid_type::Cl, "Cl"},
       {autodock_grid_type::F, "F"},
       {autodock_grid_type::S, "S"},
       {autodock_grid_type::Si, "Si"},
       {autodock_grid_type::Br, "Br"},
       {autodock_grid_type::P, "P"},
       {autodock_grid_type::I, "I"},
       {autodock_grid_type::ELEC, "Elec"},
       {autodock_grid_type::DESOLV, "Desolv"}}};
} // namespace mudock
