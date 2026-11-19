#pragma once

#include <string_view>

namespace mudock {
  struct pdbqt {
    static constexpr std::string_view molecule_name_token = "REMARK  Name =";
    std::string_view::size_type next_molecule_start_index(std::string_view text) const;
  };
} // namespace mudock
