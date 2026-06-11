#include <algorithm>
#include <mudock/chem/bond_types.hpp>
#include <stdexcept>

//===------------------------------------------------------------------------------------------------------
// WARNING: This file has been automatically generated from chem/periodic_table.json
//===------------------------------------------------------------------------------------------------------

namespace mudock {

  const std::array<bond_type_description, 5> BOND_DICTIONARY = {{
      {bond_type::SINGLE, "Single"},
      {bond_type::DOUBLE, "Double"},
      {bond_type::TRIPLE, "Triple"},
      {bond_type::AMIDE, "Amide"},
      {bond_type::AROMATIC, "Aromatic"},
  }};

  bond_type parse_bond_type(const std::string_view symbol) {
    const auto element_it = std::find_if(std::begin(BOND_DICTIONARY),
                                         std::end(BOND_DICTIONARY),
                                         [&symbol](const auto& e) { return e.name == symbol; });
    if (element_it != std::end(BOND_DICTIONARY))
      return element_it->value;
    else
      throw std::runtime_error("Missing bond type");
  }
} // namespace mudock
