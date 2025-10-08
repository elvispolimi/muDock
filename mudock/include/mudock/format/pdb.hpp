#pragma once

#include <string_view>

namespace mudock {

  class pdb {
  public:
    std::string_view::size_type next_molecule_start_index(std::string_view text) const;
  };
} // namespace mudock
