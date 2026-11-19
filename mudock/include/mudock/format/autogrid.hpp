#pragma once

#include <filesystem>
#include <mudock/chem/autodock_grid.hpp>

namespace mudock {
  [[nodiscard]] autodock_grid load_autogrid_maps(const std::filesystem::path& dpf_path);
}
