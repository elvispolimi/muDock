#pragma once

#include <cassert>
#include <memory>
#include <mudock/chem.hpp>
#include <mudock/chem/autodock_types.hpp>
#include <mudock/chem/elements.hpp>
#include <mudock/format/pdbqt.hpp>
#include <mudock/format/supported_format.hpp>
#include <mudock/log.hpp>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>
#include <mudock/utils.hpp>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/data.h>
#include <openbabel/elements.h>
#include <openbabel/generic.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/parsmart.h>
#include <sys/types.h>

namespace mudock{

  template<class molecule_type>
  void foo(molecule_type&& dest, const ob_mol_wrapper& source, size_t num_atoms) {
  /// Calculate interactive pairs if the molecule is a ligand
    if constexpr (std::same_as<std::remove_cvref_t<molecule_type>, static_molecule>) {
      for (size_t atom_id = 0; atom_id < num_atoms; ++atom_id) {
        const std::vector<int> atom_neighbors = calc_neighbors(source, atom_id);
        assert(atom_neighbors.size() <= max_static_neighbors());

        for(size_t i = 0; i < max_static_neighbors(); ++i){
          if(i < atom_neighbors.size()){
            dest.neighbors(atom_id, i) = atom_neighbors[i];
          } else {
            dest.neighbors(atom_id, i) = -1;
          }
        }
      }
    }

  }


  bool isHydrophobicAtom(const ob_mol_wrapper &mol, const OpenBabel::OBAtom *atom) {
    // Define a simple hydrophobic SMARTS: non-polar carbon, sp3
    OpenBabel::OBSmartsPattern smarts;
    smarts.Init(
        "[c,s,F,Cl,Br,I,S&H0&v2,$([D3,D4;#6])&!$([#6]~[#7,#8,#9])&!$([#6X4H0]);+0]"); // aliphatic carbon not bonded to N/O/S

    if (!smarts.Match(*mol.get()))
      return false;

    // Check if the atom is part of any match
    for (const auto &match: smarts.GetMapList()) {
      for (uint idx: match) {
        if (idx == atom->GetIdx())
          return true;
      }
    }

    return false;
  }

  std::unordered_map<int, std::vector<int>> get_atoms_in_frag(
      const std::span<const bond>& bonds, 
      const std::size_t num_atom
      ){
    auto graph = make_graph(bonds, num_atom);
    const auto ligand_fragments =
      std::make_unique<mudock::fragments<mudock::static_containers>>(graph,
          bonds,
          num_atom);

    auto rigid_pieces = ligand_fragments.get()->get_rigid_pieces();

    std::unordered_map<int, std::vector<int>> atoms_in_fragment;

    for (size_t i = 0; i < num_atom; ++i) {
      atoms_in_fragment[rigid_pieces[i]].emplace_back(i);
    }

    return atoms_in_fragment;
  }

  /// The neighbors of an atom are defined as the atoms that are connected to it a number of bonds <= 3
  std::vector<int> calc_neighbors(const ob_mol_wrapper& mol, int atomIdx) {

    std::set<int> neighbors;

    OpenBabel::OBAtom* oba0 = mol->GetAtom(atomIdx + 1); // OpenBabel uses 1-based indexing
    neighbors.insert(oba0->GetIndex());

    OpenBabel::OBBondIterator it0 = oba0->BeginBonds();
    for (OpenBabel::OBAtom* oba1 = oba0->BeginNbrAtom(it0); oba1 != nullptr; oba1 = oba0->NextNbrAtom(it0)) {
      neighbors.insert(oba1->GetIndex());
      OpenBabel::OBBondIterator it1 = oba0->BeginBonds();
      for (OpenBabel::OBAtom* oba2 = oba1->BeginNbrAtom(it1); oba2 != nullptr; oba2 = oba1->NextNbrAtom(it1)) {
        neighbors.insert(oba2->GetIndex());
        OpenBabel::OBBondIterator it2 = oba1->BeginBonds();
        for (OpenBabel::OBAtom* oba3 = oba2->BeginNbrAtom(it2); oba3 != nullptr; oba3 = oba2->NextNbrAtom(it2)) {
          neighbors.insert(oba3->GetIndex());
        }
      }
    }

    std::vector<int> out(neighbors.begin(), neighbors.end());
    return out;
  }
}
