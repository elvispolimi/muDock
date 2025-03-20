#include <algorithm>
#include <mudock/chem/ligand_maps.hpp>
#include <stdexcept>

namespace mudock {

  ligand_map_types parse_map_symbol(const std::string_view symbol) {
    const auto element_it = std::find_if(std::begin(MAP_DICTIONARY),
                                         std::end(MAP_DICTIONARY),
                                         [&symbol](const auto& e) { return e.name == symbol; });
    if (element_it != std::end(MAP_DICTIONARY))
      return element_it->value;
    else
      throw std::runtime_error("Missing map type");
  }

  const std::array<map_description, 15> MAP_DICTIONARY = {{{ligand_map_types::A, "A"},
                                                           {ligand_map_types::C, "C"},
                                                           {ligand_map_types::H, "H"},
                                                           {ligand_map_types::HD, "HD"},
                                                           {ligand_map_types::N, "N"},
                                                           {ligand_map_types::NA, "NA"},
                                                           {ligand_map_types::OA, "OA"},
                                                           {ligand_map_types::SA, "SA"},
                                                           {ligand_map_types::Cl, "Cl"},
                                                           {ligand_map_types::F, "F"},
                                                           {ligand_map_types::S, "S"},
                                                           {ligand_map_types::Br, "Br"},
                                                           {ligand_map_types::P, "P"},
                                                           {ligand_map_types::I, "I"},
                                                           {ligand_map_types::Si, "Si"}}};
} // namespace mudock
