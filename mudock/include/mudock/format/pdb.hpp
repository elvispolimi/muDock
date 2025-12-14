#pragma once

#include <string_view>

namespace mudock {

  class pdb {
  public:
    static constexpr auto PDB_START_TOKEN = "HEADER";
    std::string_view::size_type next_molecule_start_index(std::string_view text) const;
  };
} // namespace mudock
