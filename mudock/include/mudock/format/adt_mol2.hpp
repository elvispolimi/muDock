#pragma once

#include <fstream>
#include <iomanip>
#include <iostream>
#include <mudock/chem/autodock_types.hpp>
#include <mudock/molecule.hpp>
#include <string>
#include <string_view>

namespace mudock {
  struct adt_mol2_tokens {
    static constexpr auto ATOM_TOKEN     = "@<TRIPOS>ATOM";
    static constexpr auto BOND_TOKEN     = "@<TRIPOS>BOND";
    static constexpr auto MOLECULE_TOKEN = "@<TRIPOS>MOLECULE";
  };

  enum class adt_mol2_state { NONE = 0, MOLECULE, ATOM, BOND };

  class adt_mol2 {
  public:
    static constexpr auto MOLECULE_TOKEN = adt_mol2_tokens::MOLECULE_TOKEN;

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
      out_s << "GASTEIGER" << std::endl;
      out_s << std::endl;

      // Atoms
      out_s << "@<TRIPOS>ATOM" << std::endl;
      for (int atom_index = 0; atom_index < molecule.num_atoms(); ++atom_index) {
        out_s
            << std::setw(5) << atom_index << " "                                            // Atom ID
            << std::setw(8) << get_description(molecule.elements(atom_index)).symbol << " " // Atom name
            << std::setw(10) << std::fixed << std::setprecision(4) << molecule.x(atom_index) << " "
            << std::setw(10) << std::fixed << std::setprecision(4) << molecule.y(atom_index) << " "
            << std::setw(10) << std::fixed << std::setprecision(4) << molecule.z(atom_index)
            << " "
            // << std::setw(8) << get_description(molecule.elements(atom_index)).symbol << " " // Atom element
            << std::setw(8) << get_description(molecule.autodock_type(atom_index)).name << " " // Atom ADT
            << std::setw(10) << std::fixed << std::setprecision(4) << molecule.charge(atom_index) << " "
            << std::setw(10) << std::fixed << std::setprecision(0) << molecule.is_aromatic(atom_index) << " "
            << std::endl;
      }

      // Bonds
      out_s << "@<TRIPOS>BOND" << std::endl;
      for (int bond_index = 0; bond_index < molecule.num_bonds(); ++bond_index) {
        const auto bond = molecule.bonds(bond_index);
        out_s << std::setw(5) << bond_index << " " << std::setw(5) << bond.source + 1 << " " << std::setw(5)
              << bond.dest + 1 << " " << std::setw(2) << get_description(bond.type).name << " "
              << std::setw(2) << bond.can_rotate << std::endl;
      }

      out_s << std::endl;
    }

    template<class molecule_type>
      requires is_molecule<molecule_type>
    static void parse(molecule_type& molecule, const std::string_view description) {
      std::string line;
      std::istringstream desc{std::string(description)};

      adt_mol2_state state = adt_mol2_state::NONE;
      int atom_index{0}, bond_index{0};
      while (std::getline(desc, line)) {
        switch (state) {
          case adt_mol2_state::NONE:
            if (line.find(adt_mol2_tokens::MOLECULE_TOKEN) != std::string::npos) {
              state = adt_mol2_state::MOLECULE;
              std::getline(desc, line);
              molecule.properties.assign(property_type::NAME, line);
            }
            break;
          case adt_mol2_state::MOLECULE:
            if (line.find(adt_mol2_tokens::ATOM_TOKEN) != std::string::npos) {
              state      = adt_mol2_state::ATOM;
              atom_index = 0;
              bond_index = 0;
            }
            break;
          case adt_mol2_state::ATOM:
            if (line.find(adt_mol2_tokens::BOND_TOKEN) != std::string::npos) {
              state = adt_mol2_state::BOND;
            } else {
              std::istringstream stream(line);
              fp_type x, y, z, charge;
              bool is_aromatic;
              std::string _, element, adt;

              stream >> _ >> element >> x >> y >> z >> adt >> charge >> is_aromatic;

              molecule.x(atom_index)             = x;
              molecule.y(atom_index)             = y;
              molecule.z(atom_index)             = z;
              molecule.elements(atom_index)      = parse_element_symbol(element);
              molecule.autodock_type(atom_index) = parse_autodock_type(adt);
              molecule.charge(atom_index)        = charge;
              molecule.is_aromatic(atom_index)   = is_aromatic;

              atom_index += 1;
            }
            break;
          case adt_mol2_state::BOND:
            // if (line.find(adt_mol2_tokens::MOLECULE_TOKEN) != std::string::npos) {
            //   state = adt_mol2_state::MOLECULE;
            // } else
            if (line.empty()) {
              state = adt_mol2_state::NONE;
              molecule.resize(atom_index, bond_index);
            } else {
              std::istringstream stream{line};
              std::string _, bond_type;
              int atom_1, atom_2;
              bool can_rotate;

              stream >> _ >> atom_1 >> atom_2 >> bond_type >> can_rotate;

              molecule.bonds(bond_index) = {atom_1 - 1, atom_2 - 1, parse_bond_type(bond_type), can_rotate};

              bond_index += 1;
            }
          default: break;
        }
      }
      assert(state == adt_mol2_state::NONE);
    }
  };
} // namespace mudock
