#pragma once

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

    static void print(const static_molecule& molecule) {
      // Header
      std::cout << "@<TRIPOS>MOLECULE" << std::endl;
      std::cout << molecule.properties.get(property_type::NAME) << std::endl;
      std::cout << molecule.num_atoms() << " " << molecule.num_bonds() << " 0 0 0"
                << std::endl; // Atom count, bond count, etc.
      std::cout << "SMALL" << std::endl;
      std::cout << "NO_CHARGES" << std::endl;
      std::cout << std::endl;

      // Atoms
      std::cout << "@<TRIPOS>ATOM" << std::endl;
      for (int atom_index = 0; atom_index < molecule.num_atoms(); ++atom_index) {
        std::cout << std::setw(5) << atom_index << " "                                            // Atom ID
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
      std::cout << "@<TRIPOS>BOND" << std::endl;
      for (int bond_index = 0; bond_index < molecule.num_bonds(); ++bond_index) {
        const auto bond = molecule.bonds(bond_index);
        std::cout << std::setw(5) << bond_index << " " << std::setw(5) << bond.source+1 << " " << std::setw(5)
                  << bond.dest+1 << " " << std::setw(2) << get_description(bond.type).name << std::endl;
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
