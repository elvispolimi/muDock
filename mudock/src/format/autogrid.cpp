#include <array>
#include <cmath>
#include <fstream>
#include <mudock/chem/autodock_grid.hpp>
#include <mudock/format/autogrid.hpp>
#include <mudock/utils.hpp>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>

namespace mudock {
  namespace {
    std::filesystem::path resolve_reference(const std::filesystem::path& owner,
                                            const std::string& reference) {
      const std::filesystem::path path{reference};
      if (path.is_absolute()) return path;
      return owner.parent_path() / path;
    }

    std::string read_text(const std::filesystem::path& path) {
      std::ifstream file(path, std::ios::binary);
      if (!file) throw std::runtime_error("Cannot open AutoGrid input: " + path.string());
      return read_from_stream(file);
    }

    std::filesystem::path find_dpf_fld(const std::filesystem::path& dpf_path) {
      const auto dpf = read_text(dpf_path);
      std::istringstream input(dpf);
      std::string line;
      while (std::getline(input, line)) {
        std::istringstream fields(line);
        std::string directive, path;
        if (fields >> directive >> path; directive == "fld")
          return resolve_reference(dpf_path, path);
      }
      throw std::runtime_error("DPF has no fld directive: " + dpf_path.string());
    }

    autodock_grid load_fld(const std::filesystem::path& fld_path) {
      const auto fld = read_text(fld_path);
      std::istringstream input(fld);
      std::string line;
      int nx = 0, ny = 0, nz = 0;
      point3D min, max;
      bool bounds_seen = false;
      std::array<autodock_grid_type, num_autodock_grids()> types{};
      std::array<bool, num_autodock_grids()> seen_types{};
      std::array<std::filesystem::path, num_autodock_grids()> files{};
      std::size_t labels = 0;

      while (std::getline(input, line)) {
        if (line.starts_with("dim1=")) nx = std::stoi(line.substr(5));
        else if (line.starts_with("dim2=")) ny = std::stoi(line.substr(5));
        else if (line.starts_with("dim3=")) nz = std::stoi(line.substr(5));
        else if (line.starts_with("label=")) {
          if (labels >= types.size()) throw std::runtime_error("Too many channels in AutoGrid FLD");
          auto label = line.substr(6);
          if (const auto end = label.find_first_of(" \t#"); end != std::string::npos) label.resize(end);
          const auto dash = label.find('-');
          const auto symbol = label.substr(0, dash);
          if (symbol == "Electrostatics") types[labels] = autodock_grid_type::ELEC;
          else if (symbol == "Desolvation") types[labels] = autodock_grid_type::DESOLV;
          else types[labels] = parse_map_symbol(symbol);
          const auto type_index = static_cast<std::size_t>(types[labels]);
          if (seen_types[type_index]) throw std::runtime_error("Duplicate channel in AutoGrid FLD: " + symbol);
          seen_types[type_index] = true;
          ++labels;
        } else if (line.starts_with("coord 1 ")) {
          if (line.find("filetype=ascii") == std::string::npos)
            throw std::runtime_error("Only ASCII AutoGrid coordinate files are supported");
          const std::regex file_pattern{R"(file=([^\s]+))"};
          std::smatch match;
          if (!std::regex_search(line, match, file_pattern))
            throw std::runtime_error("AutoGrid FLD coord 1 lacks file reference");
          const auto xyz_path = resolve_reference(fld_path, match[1]);
          std::istringstream xyz(read_text(xyz_path));
          std::string xyz_line;
          std::array<fp_type, 3> mins{}, maxs{};
          for (int axis = 0; axis < 3; ++axis) {
            if (!std::getline(xyz, xyz_line)) throw std::runtime_error("Malformed AutoGrid XYZ file");
            std::istringstream values(xyz_line);
            if (!(values >> mins[static_cast<std::size_t>(axis)] >> maxs[static_cast<std::size_t>(axis)]))
              throw std::runtime_error("Malformed AutoGrid XYZ bounds");
          }
          min = point3D{mins[0], mins[1], mins[2]};
          max = point3D{maxs[0], maxs[1], maxs[2]};
          bounds_seen = true;
        } else if (line.starts_with("variable ")) {
          if (line.find("filetype=ascii") == std::string::npos)
            throw std::runtime_error("Only ASCII AutoGrid map files are supported");
          std::istringstream fields(line);
          std::string variable, file_token;
          std::size_t id = 0;
          fields >> variable >> id >> file_token;
          if (id == 0 || id > labels) throw std::runtime_error("AutoGrid FLD variable id is out of range");
          constexpr std::string_view prefix{"file="};
          if (!file_token.starts_with(prefix)) throw std::runtime_error("Malformed AutoGrid FLD file entry");
          files[static_cast<std::size_t>(types[id - 1])] =
              resolve_reference(fld_path, file_token.substr(prefix.size()));
        }
      }

      if (nx < 2 || ny < 2 || nz < 2) throw std::runtime_error("AutoGrid FLD has invalid dimensions");
      if (!bounds_seen) throw std::runtime_error("AutoGrid FLD has no coord 1 bounds");
      autodock_grid grid{md_index<3>{nx, ny, nz}, min, max};
      if (!seen_types[static_cast<std::size_t>(autodock_grid_type::ELEC)] ||
          !seen_types[static_cast<std::size_t>(autodock_grid_type::DESOLV)])
        throw std::runtime_error("AutoGrid FLD must contain electrostatic and desolvation maps");
      for (std::size_t channel = 0; channel < labels; ++channel) {
        const auto type = types[channel];
        const auto path = files[static_cast<std::size_t>(type)];
        if (path.empty()) throw std::runtime_error("Missing AutoGrid map for channel " +
                                                   std::string(get_description(type).name));
        std::istringstream map(read_text(path));
        for (int header = 0; header < 6; ++header)
          if (!std::getline(map, line)) throw std::runtime_error("Truncated AutoGrid map header: " + path.string());
        auto view = grid.get_atom_map(type);
        for (std::size_t z = 0; z < grid.index.size_z(); ++z)
          for (std::size_t y = 0; y < grid.index.size_y(); ++y)
            for (std::size_t x = 0; x < grid.index.size_x(); ++x) {
              fp_type value{};
              if (!(map >> value) || !std::isfinite(value))
                throw std::runtime_error("Truncated or malformed AutoGrid map: " + path.string());
              view.get(x, y, z) = value;
            }
        grid.mark_map_present(type);
      }
      return grid;
    }
  } // namespace

  autodock_grid load_autogrid_maps(const std::filesystem::path& dpf_path) {
    return load_fld(find_dpf_fld(dpf_path));
  }
} // namespace mudock
