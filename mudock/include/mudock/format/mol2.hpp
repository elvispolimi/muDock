#pragma once

#include <fstream>
#include <iomanip>
#include <iostream>
#include <mudock/molecule.hpp>
#include <string_view>

namespace mudock {

  class mol2 {
  public:
    std::string_view::size_type next_molecule_start_index(std::string_view text) const;

    static void print(const static_molecule& molecule, std::ostream& out_s) {
      // Header
      out_s << "@<TRIPOS>MOLECULE" << std::endl;
      out_s << molecule.properties.get(property_type::NAME) << std::endl;
      out_s << molecule.num_atoms() << " " << molecule.num_bonds() << " 0 0 0"
            << std::endl; // Atom count, bond count, etc.
      out_s << "SMALL" << std::endl;
      out_s << "NO_CHARGES" << std::endl;
      out_s << std::endl;

      // Atoms
      out_s << "@<TRIPOS>ATOM" << std::endl;
      for (int atom_index = 0; atom_index < molecule.num_atoms(); ++atom_index) {
        out_s << std::setw(5) << atom_index << " "                                            // Atom ID
              << std::setw(8) << get_description(molecule.elements(atom_index)).symbol << " " // Atom name
              << std::setw(10) << std::fixed << std::setprecision(4) << molecule.x(atom_index) << " "
              << std::setw(10) << std::fixed << std::setprecision(4) << molecule.y(atom_index) << " "
              << std::setw(10) << std::fixed << std::setprecision(4) << molecule.z(atom_index) << " "
              << std::setw(8) << get_description(molecule.elements(atom_index)).symbol
              << " " // Atom type in MOL2 format
              // << std::setw(5) << "1"
              // << " "                                   // Assume molecule number is 1
              // << molecule.name << std::endl;
              << std::endl;
      }

      // Bonds
      out_s << "@<TRIPOS>BOND" << std::endl;
      for (int bond_index = 0; bond_index < molecule.num_bonds(); ++bond_index) {
        const auto bond = molecule.bonds(bond_index);
        out_s << std::setw(5) << bond_index << " " << std::setw(5) << bond.source + 1 << " " << std::setw(5)
              << bond.dest + 1 << " " << std::setw(2) << get_description(bond.type).name << std::endl;
      }
    };
  };

} // namespace mudock
