#pragma once

#include "mudock/chem/residue_types.hpp"
#include "mudock/chem/sybyl_atom_types.hpp"

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
        const auto element_symbol = std::string_view{get_description(molecule.elements(atom_index)).symbol};

        const std::string generated_atom_name = std::string(element_symbol) + std::to_string(atom_index + 1);

        const auto& stored_atom_name = molecule.atom_name(atom_index);

        const auto atom_name = stored_atom_name.empty() ? std::string_view{generated_atom_name}
                                                        : std::string_view{stored_atom_name};

        const auto sybyl_type = molecule.sybyl_type(atom_index);

        const auto atom_type =
            sybyl_type == sybyl_atom_type::UNKNOWN ? element_symbol : to_string(sybyl_type);

        const auto residue_id        = molecule.residue_id(atom_index);
        const auto& residue_name     = molecule.residue_name(atom_index);
        const auto residue_type_name = to_string(molecule.atom_residue_type(atom_index));
        out_s << std::setw(5) << atom_index + 1 << " " << std::setw(8) << atom_name << " " << std::fixed
              << std::setprecision(4) << std::setw(10) << molecule.x(atom_index) << " " << std::setw(10)
              << molecule.y(atom_index) << " " << std::setw(10) << molecule.z(atom_index) << " "
              << std::setw(8) << mudock::to_string(molecule.sybyl_type(atom_index)) << " " << std::setw(5)
              << molecule.residue_id(atom_index) << " " << std::setw(8) << molecule.residue_name(atom_index)
              << " " << std::setw(8) << get_description(molecule.autodock_type(atom_index)).name << " "
              << std::setw(10) << molecule.charge(atom_index) << " " << std::setw(10)
              << static_cast<int>(molecule.is_aromatic(atom_index)) << std::endl;
      }

      // Bonds
      out_s << "@<TRIPOS>BOND" << std::endl;
      for (int bond_index = 0; bond_index < molecule.num_bonds(); ++bond_index) {
        const auto bond = molecule.bonds(bond_index);
        out_s << std::setw(5) << bond_index + 1 << " " << std::setw(5) << bond.source + 1 << " "
              << std::setw(5) << bond.dest + 1 << " " << std::setw(2) << get_description(bond.type).name
              << " " << std::setw(2) << bond.can_rotate << std::endl;
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

              // Molecule name line
              if (!std::getline(desc, line)) {
                throw std::runtime_error("Invalid ADT-MOL2 file: missing molecule name");
              }

              molecule.properties.assign(property_type::NAME, line);

              // Counts line: num_atoms num_bonds ...
              if (!std::getline(desc, line)) {
                throw std::runtime_error("Invalid ADT-MOL2 file: missing molecule counts line");
              }

              std::istringstream counts_stream(line);

              int num_atoms = 0;
              int num_bonds = 0;

              counts_stream >> num_atoms >> num_bonds;

              if (!counts_stream || num_atoms < 0 || num_bonds < 0) {
                throw std::runtime_error("Invalid ADT-MOL2 file: invalid molecule counts line: " + line);
              }

              molecule.resize(num_atoms, num_bonds);

              atom_index = 0;
              bond_index = 0;
            }
            break;

          case adt_mol2_state::MOLECULE:
            if (line.find(adt_mol2_tokens::ATOM_TOKEN) != std::string::npos) {
              state      = adt_mol2_state::ATOM;
              atom_index = 0;
            }
            break;

          case adt_mol2_state::ATOM:
            if (line.find(adt_mol2_tokens::BOND_TOKEN) != std::string::npos) {
              if (atom_index != molecule.num_atoms()) {
                throw std::runtime_error("Invalid ADT-MOL2 file: expected " +
                                         std::to_string(molecule.num_atoms()) + " atoms, parsed " +
                                         std::to_string(atom_index));
              }

              state      = adt_mol2_state::BOND;
              bond_index = 0;
            } else if (!line.empty()) {
              if (atom_index >= molecule.num_atoms()) {
                throw std::runtime_error("Invalid ADT-MOL2 file: too many atom records");
              }

              std::istringstream stream(line);

              int atom_id    = 0;
              int residue_id = 0;

              fp_type x      = 0;
              fp_type y      = 0;
              fp_type z      = 0;
              fp_type charge = 0;

              bool is_aromatic = false;

              std::string atom_name;
              std::string element;
              std::string adt;
              std::string residue_name;
              std::string sybyl_type;
              std::string residue_type_name;

              stream >> atom_id >> atom_name >> x >> y >> z >> sybyl_type >> residue_id >> residue_name >>
                  residue_type_name >> adt >> charge >> is_aromatic;

              if (!stream) {
                throw std::runtime_error("Invalid ADT-MOL2 atom record: " + line);
              }

              // Optional sanity check: atom IDs in the file should normally be 1-based.
              if (atom_id != atom_index + 1) {
                throw std::runtime_error("Invalid ADT-MOL2 atom record: unexpected atom id " +
                                         std::to_string(atom_id) + ", expected " +
                                         std::to_string(atom_index + 1));
              }

              molecule.x(atom_index) = x;
              molecule.y(atom_index) = y;
              molecule.z(atom_index) = z;

              const auto parsed_sybyl_type       = parse_sybyl_atom_type(sybyl_type);
              molecule.elements(atom_index)      = get_element(parsed_sybyl_type);
              molecule.autodock_type(atom_index) = parse_autodock_type(adt);

              molecule.residue_id(atom_index)        = residue_id;
              molecule.residue_name(atom_index)      = residue_name;
              molecule.atom_residue_type(atom_index) = parse_residue_type(residue_type_name);

              molecule.atom_name(atom_index)  = atom_name;
              molecule.sybyl_type(atom_index) = parsed_sybyl_type;

              molecule.charge(atom_index)      = charge;
              molecule.is_aromatic(atom_index) = is_aromatic;

              atom_index += 1;
            }
            break;

          case adt_mol2_state::BOND:
            if (line.empty()) {
              if (bond_index != molecule.num_bonds()) {
                throw std::runtime_error("Invalid ADT-MOL2 file: expected " +
                                         std::to_string(molecule.num_bonds()) + " bonds, parsed " +
                                         std::to_string(bond_index));
              }

              state = adt_mol2_state::NONE;
            } else {
              if (bond_index >= molecule.num_bonds()) {
                throw std::runtime_error("Invalid ADT-MOL2 file: too many bond records");
              }

              std::istringstream stream{line};

              int bond_id     = 0;
              int atom_1      = 0;
              int atom_2      = 0;
              bool can_rotate = false;
              std::string bond_type;

              stream >> bond_id >> atom_1 >> atom_2 >> bond_type >> can_rotate;

              if (!stream) {
                throw std::runtime_error("Invalid ADT-MOL2 bond record: " + line);
              }

              // Optional sanity check: bond IDs should normally be 1-based.
              if (bond_id != bond_index + 1) {
                throw std::runtime_error("Invalid ADT-MOL2 bond record: unexpected bond id " +
                                         std::to_string(bond_id) + ", expected " +
                                         std::to_string(bond_index + 1));
              }

              if (atom_1 <= 0 || atom_1 > molecule.num_atoms() || atom_2 <= 0 ||
                  atom_2 > molecule.num_atoms()) {
                throw std::runtime_error("Invalid ADT-MOL2 bond record: atom index out of range: " + line);
              }

              molecule.bonds(bond_index) = {atom_1 - 1, atom_2 - 1, parse_bond_type(bond_type), can_rotate};

              bond_index += 1;
            }
            break;
        }
      }
      if (state == adt_mol2_state::BOND) {
        if (bond_index != molecule.num_bonds()) {
          throw std::runtime_error("Invalid ADT-MOL2 file: expected " + std::to_string(molecule.num_bonds()) +
                                   " bonds, parsed " + std::to_string(bond_index));
        }

        state = adt_mol2_state::NONE;
      }

      if (state != adt_mol2_state::NONE) {
        throw std::runtime_error("Invalid ADT-MOL2 file: unexpected end of file");
      }
    }
  };
} // namespace mudock
