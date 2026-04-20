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

    inline std::string mol2_element_token(const std::string_view atom_type, const std::string_view atom_name) {
      auto normalize = [](std::string token) {
        token = trim_copy(std::move(token));
        const auto dot_index = token.find('.');
        if (dot_index != std::string::npos) {
          token.resize(dot_index);
        }
        token.erase(std::remove_if(token.begin(), token.end(), [](unsigned char c) {
                      return !std::isalpha(c);
                    }),
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
        out_s << std::setw(5) << atom_index + 1 << " "                                        // Atom ID
              << std::setw(8) << get_description(molecule.elements(atom_index)).symbol << " " // Atom name
              << std::setw(10) << std::fixed << std::setprecision(4) << molecule.x(atom_index) << " "
              << std::setw(10) << std::fixed << std::setprecision(4) << molecule.y(atom_index) << " "
              << std::setw(10) << std::fixed << std::setprecision(4) << molecule.z(atom_index) << " "
              << std::setw(8) << get_description(molecule.elements(atom_index)).symbol << " " // Atom type
              << std::setw(5) << 1 << " "
              << std::setw(8) << "LIG" << " "
              << std::setw(10) << std::fixed << std::setprecision(4) << molecule.charge(atom_index)
              << std::endl;
      }

      // Bonds
      out_s << "@<TRIPOS>BOND" << std::endl;
      for (int bond_index = 0; bond_index < molecule.num_bonds(); ++bond_index) {
        const auto bond = molecule.bonds(bond_index);
        out_s << std::setw(5) << bond_index + 1 << " " << std::setw(5) << bond.source + 1 << " " << std::setw(5)
              << bond.dest + 1 << " " << std::setw(2) << detail::print_mol2_bond_type(bond.type) << std::endl;
      }
      out_s << std::endl;
    }

    template<class molecule_type>
      requires is_molecule<molecule_type>
    static void parse(molecule_type& molecule, const std::string_view description) {
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
            std::istringstream stream(line);
            int atom_id = 0;
            std::string atom_name;
            fp_type x = 0, y = 0, z = 0;
            std::string atom_type;
            int subst_id = 0;
            std::string subst_name;
            fp_type charge = 0;

            stream >> atom_id >> atom_name >> x >> y >> z >> atom_type;
            if (!stream) {
              throw std::runtime_error("Invalid MOL2 atom record");
            }
            stream >> subst_id >> subst_name >> charge;

            molecule.x(atom_index)      = x;
            molecule.y(atom_index)      = y;
            molecule.z(atom_index)      = z;
            molecule.charge(atom_index) = charge;
            molecule.is_aromatic(atom_index) =
                detail::is_mol2_aromatic_atom(atom_type) ? 1 : 0;
            molecule.elements(atom_index) =
                parse_element_symbol(detail::mol2_element_token(atom_type, atom_name));
            ++atom_index;
            break;
          }
          case state::BOND: {
            if (bond_index >= expected_bonds) {
              throw std::runtime_error("MOL2 bond section exceeds declared bond count");
            }
            std::istringstream stream(line);
            int bond_id = 0;
            int atom_1 = 0;
            int atom_2 = 0;
            std::string bond_type_token;

            stream >> bond_id >> atom_1 >> atom_2 >> bond_type_token;
            if (!stream) {
              throw std::runtime_error("Invalid MOL2 bond record");
            }

            molecule.bonds(bond_index) =
                {atom_1 - 1, atom_2 - 1, detail::parse_mol2_bond_type(bond_type_token), false};
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
        auto& bond = molecule.bonds(i);
        bond.can_rotate = (bond.type == bond_type::SINGLE) &&
                          (degrees[static_cast<std::size_t>(bond.source)] > 1) &&
                          (degrees[static_cast<std::size_t>(bond.dest)] > 1);
      }
    }
  };

} // namespace mudock
