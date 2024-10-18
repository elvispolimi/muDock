#pragma once

#include <fstream>
#include <mudock/format/rdkit_wrapper.hpp>
#include <mudock/molecule.hpp>
#include <string_view>

namespace mudock {

  class mol2 {
  public:
    std::string_view::size_type next_molecule_start_index(std::string_view text) const;

    template<class molecule_type>
      requires is_molecule<molecule_type>
    void parse(molecule_type&& molecule, std::string_view molecule_description) const;

    // TODO rework print function
    // Some data are getting lost during the conversion to rdkit
    // TODO find a way to use rdkit function
    static void print(const static_molecule& molecule, std::ofstream& output) {
      // Header
      output << "@<TRIPOS>MOLECULE" << std::endl;
      output << molecule.properties.get(property_type::NAME) << std::endl;
      output << molecule.num_atoms() << " " << molecule.num_bonds() << " 0 0 0"
             << std::endl; // Atom count, bond count, etc.
      output << "SMALL" << std::endl;
      output << "NO_CHARGES" << std::endl;
      output << std::endl;

      // Atoms
      output << "@<TRIPOS>ATOM" << std::endl;
      for (int atom_index = 0; atom_index < molecule.num_atoms(); ++atom_index) {
        output << std::setw(5) << atom_index << " "                                            // Atom ID
               << std::setw(8) << get_description(molecule.elements(atom_index)).symbol << " " // Atom name
               << std::setw(10) << std::fixed << std::setprecision(4) << molecule.x(atom_index) << " "
               << std::setw(10) << std::fixed << std::setprecision(4) << molecule.y(atom_index) << " "
               << std::setw(10) << std::fixed << std::setprecision(4) << molecule.z(atom_index) << " "
               << std::setw(8) << get_description(molecule.elements(atom_index)).symbol << std::endl;
      }

      // Bonds
      auto get_bond_type_string = [](const bond_type& b_t) {
        switch (b_t) {
          case bond_type::SINGLE: return "1";    // Single bond
          case bond_type::DOUBLE: return "2";    // Double bond
          case bond_type::TRIPLE: return "3";    // Triple bond
          case bond_type::AROMATIC: return "ar"; // Aromatic bond
          case bond_type::AMIDE: return "am";    // Amide bond (if applicable)
          default: return "un";                  // Unknown bond type (MOL2 often uses "un" for unknown)
        }
      };
      output << "@<TRIPOS>BOND" << std::endl;
      for (int bond_index = 0; bond_index < molecule.num_bonds(); ++bond_index) {
        const auto bond = molecule.bonds(bond_index);
        output << std::setw(5) << bond_index << " " << std::setw(5) << bond.source + 1 << " " << std::setw(5)
               << bond.dest + 1 << " " << std::setw(2) << get_bond_type_string(bond.type) << std::endl;
      }
    };
  };

  //===------------------------------------------------------------------------------------------------------
  // Out-of-class method definitions
  //===------------------------------------------------------------------------------------------------------

  template<class molecule_type>
    requires is_molecule<molecule_type>
  void mol2::parse(molecule_type&& molecule, std::string_view molecule_description) const {
    // NOTE: we use RDKit to parse everything and deal with all the chemistry complexity
    convert(molecule, parse_mol2(molecule_description));
  }

} // namespace mudock
