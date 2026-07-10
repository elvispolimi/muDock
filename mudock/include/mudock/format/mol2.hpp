#pragma once

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mudock/chem/bond_types.hpp>
#include <mudock/chem/elements.hpp>
#include <mudock/molecule.hpp>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace mudock {
  struct mol2_tokens {
    static constexpr auto MOLECULE_TOKEN = "@<TRIPOS>MOLECULE";
  };

  namespace detail {
    inline std::string trim_copy(std::string value) {
      const auto begin = value.find_first_not_of(" \t\r\n");
      if (begin == std::string::npos) {
        return {};
      }
      const auto end = value.find_last_not_of(" \t\r\n");
      return value.substr(begin, end - begin + 1);
    }

    inline std::string mol2_element_token(const std::string_view atom_type,
                                          const std::string_view atom_name) {
      auto normalize = [](std::string token) {
        token                = trim_copy(std::move(token));
        const auto dot_index = token.find('.');
        if (dot_index != std::string::npos) {
          token.resize(dot_index);
        }
        token.erase(
            std::remove_if(token.begin(), token.end(), [](unsigned char c) { return !std::isalpha(c); }),
            token.end());
        if (token.empty()) {
          return token;
        }
        token[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(token[0])));
        for (std::size_t i = 1; i < token.size(); ++i) {
          token[i] = static_cast<char>(std::tolower(static_cast<unsigned char>(token[i])));
        }
        if (token.size() > 2) {
          token.resize(2);
        }
        return token;
      };

      auto token = normalize(std::string(atom_type));
      if (!token.empty()) {
        return token;
      }
      return normalize(std::string(atom_name));
    }

    inline bond_type parse_mol2_bond_type(std::string token) {
      token = trim_copy(std::move(token));
      std::transform(token.begin(), token.end(), token.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
      });
      if (token == "1" || token == "un") {
        return bond_type::SINGLE;
      }
      if (token == "2") {
        return bond_type::DOUBLE;
      }
      if (token == "3") {
        return bond_type::TRIPLE;
      }
      if (token == "am") {
        return bond_type::AMIDE;
      }
      if (token == "ar") {
        return bond_type::AROMATIC;
      }
      throw std::runtime_error("Unsupported MOL2 bond type '" + token + "'");
    }

    inline const char* print_mol2_bond_type(const bond_type type) {
      switch (type) {
        case bond_type::SINGLE: return "1";
        case bond_type::DOUBLE: return "2";
        case bond_type::TRIPLE: return "3";
        case bond_type::AMIDE: return "am";
        case bond_type::AROMATIC: return "ar";
        default: throw std::runtime_error("Unsupported internal bond type for MOL2 output");
      }
    }

    inline bool is_mol2_aromatic_atom(const std::string& atom_type) {
      return atom_type.find(".ar") != std::string::npos || atom_type.find(".AR") != std::string::npos;
    }
  } // namespace detail

  class mol2 {
  public:
    static constexpr auto MOLECULE_TOKEN = mol2_tokens::MOLECULE_TOKEN;

    std::string_view::size_type next_molecule_start_index(std::string_view text) const;

    template<class molecule_type>
      requires is_molecule<molecule_type>
    static void print(const molecule_type& molecule, std::ostream& out_s) {
      // Header
      out_s << "@<TRIPOS>MOLECULE" << std::endl;
      out_s << molecule.properties.get(property_type::NAME) << std::endl;
      out_s << molecule.num_atoms() << " " << molecule.num_bonds() << " 0 0 0"
            << std::endl; // Atom count, bond count, etc.
      out_s << "SMALL" << std::endl;
      out_s << "USER_CHARGES" << std::endl;
      out_s << std::endl;

      // Atoms
      out_s << "@<TRIPOS>ATOM" << std::endl;
      for (int atom_index = 0; atom_index < molecule.num_atoms(); ++atom_index) {
        const auto element_symbol = std::string_view{get_description(molecule.elements(atom_index)).symbol};

        const std::string generated_atom_name = std::string(element_symbol) + std::to_string(atom_index + 1);

        const auto& stored_atom_name = molecule.atom_name(atom_index);

        const auto atom_name = stored_atom_name.empty() ? std::string_view{generated_atom_name}
                                                        : std::string_view{stored_atom_name};

        const auto sybyl_type = molecule.sybyl_type(atom_index);

        const auto atom_type =
            sybyl_type == sybyl_atom_type::UNKNOWN ? element_symbol : to_string(sybyl_type);

        const auto residue_id    = molecule.residue_id(atom_index);
        const auto& residue_name = molecule.residue_name(atom_index);
        out_s << std::setw(5) << atom_index + 1 << " " << std::setw(8) << atom_name << " " << std::setw(10)
              << std::fixed << std::setprecision(4) << molecule.x(atom_index) << " " << std::setw(10)
              << std::fixed << std::setprecision(4) << molecule.y(atom_index) << " " << std::setw(10)
              << std::fixed << std::setprecision(4) << molecule.z(atom_index) << " " << std::setw(8)
              << atom_type << " " << std::setw(5) << residue_id << " " << std::setw(8) << residue_name << " "
              << std::setw(10) << std::fixed << std::setprecision(4) << molecule.charge(atom_index)
              << std::endl;
      }

      // Bonds
      out_s << "@<TRIPOS>BOND" << std::endl;
      for (int bond_index = 0; bond_index < molecule.num_bonds(); ++bond_index) {
        const auto bond = molecule.bonds(bond_index);
        out_s << std::setw(5) << bond_index + 1 << " " << std::setw(5) << bond.source + 1 << " "
              << std::setw(5) << bond.dest + 1 << " " << std::setw(2)
              << detail::print_mol2_bond_type(bond.type) << std::endl;
      }
      out_s << std::endl;
    }

    template<class molecule_type>
      requires is_molecule<molecule_type>
    static void parse(molecule_type& molecule, const std::string_view description) {
      // NOTE:
      // This is the native MOL2 path used for internal molecules. It preserves the raw MOL2 atom typing and
      // computes can_rotate with a local heuristic instead of forwarding OpenBabel rotor callbacks.
      enum class state { NONE, MOLECULE, ATOM, BOND };

      state current_state = state::NONE;
      std::string line;
      std::istringstream desc{std::string(description)};
      int expected_atoms = 0;
      int expected_bonds = 0;
      int atom_index     = 0;
      int bond_index     = 0;

      while (std::getline(desc, line)) {
        line = detail::trim_copy(std::move(line));
        if (line.empty()) {
          continue;
        }
        if (line == "@<TRIPOS>MOLECULE") {
          current_state = state::MOLECULE;
          atom_index    = 0;
          bond_index    = 0;
          continue;
        }
        if (line == "@<TRIPOS>ATOM") {
          current_state = state::ATOM;
          continue;
        }
        if (line == "@<TRIPOS>BOND") {
          current_state = state::BOND;
          continue;
        }
        if (line.rfind("@<TRIPOS>", 0) == 0) {
          current_state = state::NONE;
          continue;
        }

        switch (current_state) {
          case state::NONE: break;
          case state::MOLECULE: {
            if (expected_atoms == 0 && molecule.properties.get(property_type::NAME) == "N/A") {
              molecule.properties.assign(property_type::NAME, detail::trim_copy(line));
              continue;
            }

            if (expected_atoms == 0 && expected_bonds == 0) {
              std::istringstream stream(line);
              stream >> expected_atoms >> expected_bonds;
              if (expected_atoms <= 0 || expected_bonds < 0) {
                throw std::runtime_error("Invalid MOL2 counts line");
              }
              molecule.resize(expected_atoms, expected_bonds);
            }
            break;
          }
          case state::ATOM: {
            if (atom_index >= expected_atoms) {
              throw std::runtime_error("MOL2 atom section exceeds declared atom count");
            }
            try {
              std::istringstream stream(line);
              std::vector<std::string> tokens;
              for (std::string token; stream >> token;) {
                tokens.push_back(std::move(token));
              }

              if (tokens.size() < 6) {
                throw std::runtime_error("Invalid MOL2 atom record: too few fields");
              }

              const int atom_id = std::stoi(tokens[0]);
              if (atom_id != atom_index + 1) {
                throw std::runtime_error("Invalid MOL2 atom record: unexpected atom id " +
                                         std::to_string(atom_id) + ", expected " +
                                         std::to_string(atom_index + 1));
              }

              const auto& atom_name  = tokens[1];
              const fp_type x        = static_cast<fp_type>(std::stod(tokens[2]));
              const fp_type y        = static_cast<fp_type>(std::stod(tokens[3]));
              const fp_type z        = static_cast<fp_type>(std::stod(tokens[4]));
              const auto& atom_type  = tokens[5];

              molecule.x(atom_index)           = x;
              molecule.y(atom_index)           = y;
              molecule.z(atom_index)           = z;
              molecule.is_aromatic(atom_index) = detail::is_mol2_aromatic_atom(atom_type) ? 1 : 0;
              molecule.atom_name(atom_index)   = atom_name;
              molecule.sybyl_type(atom_index)  = parse_sybyl_atom_type(atom_type);
              molecule.elements(atom_index) =
                  parse_element_symbol(detail::mol2_element_token(atom_type, atom_name));

              if (tokens.size() >= 7) {
                molecule.residue_id(atom_index) = std::stoi(tokens[6]);
              }
              if (tokens.size() >= 8) {
                molecule.residue_name(atom_index)      = tokens[7];
                molecule.atom_residue_type(atom_index) = parse_residue_type(tokens[7]);
              }
              if (tokens.size() >= 9) {
                molecule.charge(atom_index) = static_cast<fp_type>(std::stod(tokens[8]));
              }

              ++atom_index;
            } catch (const std::exception& e) {
              throw std::runtime_error("Invalid MOL2 atom record: " + line + " (" + e.what() + ")");
            }
            break;
          }
          case state::BOND: {
            if (bond_index >= expected_bonds) {
              throw std::runtime_error("MOL2 bond section exceeds declared bond count");
            }
            std::istringstream stream(line);
            int bond_id = 0;
            int atom_1  = 0;
            int atom_2  = 0;
            std::string bond_type_token;

            stream >> bond_id >> atom_1 >> atom_2 >> bond_type_token;
            if (!stream) {
              throw std::runtime_error("Invalid MOL2 bond record");
            }

            molecule.bonds(
                bond_index) = {atom_1 - 1, atom_2 - 1, detail::parse_mol2_bond_type(bond_type_token), false};
            ++bond_index;
            break;
          }
        }
      }

      if (expected_atoms == 0) {
        throw std::runtime_error("MOL2 molecule is missing a valid counts line");
      }
      if (atom_index != expected_atoms || bond_index != expected_bonds) {
        throw std::runtime_error("MOL2 molecule size does not match declared counts");
      }

      std::vector<int> degrees(static_cast<std::size_t>(molecule.num_atoms()), 0);
      for (int i = 0; i < molecule.num_bonds(); ++i) {
        const auto& bond = molecule.bonds(i);
        ++degrees[static_cast<std::size_t>(bond.source)];
        ++degrees[static_cast<std::size_t>(bond.dest)];
      }
      for (int i = 0; i < molecule.num_bonds(); ++i) {
        auto& bond      = molecule.bonds(i);
        bond.can_rotate = (bond.type == bond_type::SINGLE) &&
                          (degrees[static_cast<std::size_t>(bond.source)] > 1) &&
                          (degrees[static_cast<std::size_t>(bond.dest)] > 1);
      }
    }
  };

} // namespace mudock
