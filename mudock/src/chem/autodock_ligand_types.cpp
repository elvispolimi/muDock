#include <algorithm>
#include <mudock/chem/autodock_ligand_types.hpp>
#include <stdexcept>

namespace mudock {

  autodock_ligand_ff parse_map_symbol(const std::string_view symbol) {
    const auto element_it = std::find_if(std::begin(MAP_DICTIONARY),
                                         std::end(MAP_DICTIONARY),
                                         [&symbol](const auto& e) { return e.name == symbol; });
    if (element_it != std::end(MAP_DICTIONARY))
      return element_it->value;
    else
      throw std::runtime_error("Missing map type");
  }

  const std::array<map_description, num_autodock_ligand_types()> MAP_DICTIONARY = {
      {{autodock_ligand_ff::A, "A"},
       {autodock_ligand_ff::C, "C"},
       {autodock_ligand_ff::H, "H"},
       {autodock_ligand_ff::HD, "HD"},
       {autodock_ligand_ff::N, "N"},
       {autodock_ligand_ff::NA, "NA"},
       {autodock_ligand_ff::OA, "OA"},
       {autodock_ligand_ff::SA, "SA"},
       {autodock_ligand_ff::Cl, "Cl"},
       {autodock_ligand_ff::F, "F"},
       {autodock_ligand_ff::S, "S"},
       {autodock_ligand_ff::Br, "Br"},
       {autodock_ligand_ff::P, "P"},
       {autodock_ligand_ff::I, "I"}}};
} // namespace mudock
