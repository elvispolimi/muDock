#pragma once

#include <mudock/chem/autodock_parameters.hpp>
#include <mudock/chem/autodock_types.hpp>
#include <mudock/chem/residue.hpp>
#include <mudock/chem/x_score_xtool_types.hpp>
#include <mudock/chem/residue_types.hpp>
#include <mudock/chem/sybyl_atom_types.hpp>
#include <mudock/format/ob_helper.hpp>
#include <mudock/molecule.hpp>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/data.h>
#include <openbabel/elements.h>
#include <openbabel/generic.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/oberror.h>
#include <openbabel/residue.h>

namespace mudock {

  //===------------------------------------------------------------------------------------------------------
  // Useful type alias to work with OpenBabel data structures
  //===------------------------------------------------------------------------------------------------------

  using ob_mol_wrapper = std::unique_ptr<OpenBabel::OBMol>;

  [[nodiscard]] bond_type parse_ob_bond_type(const OpenBabel::OBBond& bond_type);

  [[nodiscard]] bool ob_rotate_check(::OpenBabel::OBBond&);
  [[nodiscard]] bool pdbqt_rotate_check(::OpenBabel::OBBond&);

  [[nodiscard]] autodock_ff convert_ob_elements(OpenBabel::OBAtom&);

  template<class molecule_type>
    requires is_molecule<molecule_type>
  // TODO use this for the default case of the molecule
  void apply_autodock_forcefield_ob(molecule_type&& molecule, const ob_mol_wrapper& ob_mol) {
    int index = 0;
    for (auto atom_it = ob_mol->BeginAtoms(); atom_it < ob_mol->EndAtoms(); ++atom_it) {
      const auto atom = *atom_it;

      autodock_ff adt               = convert_ob_elements(*atom);
      molecule.autodock_type(index) = adt;
      const auto& ff_entry          = get_description(adt);
      molecule.autodock_type(index) = ff_entry.value;
      molecule.Rii(index)           = ff_entry.Rii;
      molecule.epsii(index)         = ff_entry.epsii * autodock_parameters::coeff_vdW;
      molecule.vol(index)           = ff_entry.vol;
      molecule.solpar(index)        = ff_entry.solpar;
      molecule.Rij_hb(index)        = ff_entry.Rij_hb;
      molecule.epsij_hb(index)      = ff_entry.epsij_hb * autodock_parameters::coeff_hbond;
      molecule.num_hbond(index)     = ff_entry.hbond;
      ++index;
    }
  } // namespace mudock

  //===------------------------------------------------------------------------------------------------------
  // Translate an OpenBabel molecule to our internal format
  //===------------------------------------------------------------------------------------------------------

  template<class molecule_type>
    requires derived_from_molecule<molecule_type>
  void convert(molecule_type& dest,
               const ob_mol_wrapper& source,
               std::function<bool(OpenBabel::OBBond&)> check_rotor_bond = ob_rotate_check) {
    const size_t num_atoms = source->NumAtoms();
    const size_t num_bonds = source->NumBonds();
    // set the molecule geometry
    // NOTE: for static molecules we need to enforce the constraint on the maximum number
    //       ot atoms or bonds by throwing an exception
    if constexpr (std::same_as<std::remove_cvref_t<molecule_type>, static_molecule>) {
      if (num_atoms > max_static_atoms() || num_bonds > max_static_bonds()) {
        throw std::runtime_error("Number of atoms or bonds exceeding static storage");
      }
    }
    dest.resize(static_cast<int>(num_atoms), static_cast<int>(num_bonds));

    // TODO check charges
    int mudock_atom_index{0};
    std::unordered_map<unsigned int, std::size_t> index_translator;

    OpenBabel::OBTypeTable ttab;
    ttab.SetFromType("INT");
    // TODO add flags to customize this value
    ttab.SetToType("XYZ");

    OpenBabel::OBTypeTable sybyl_table;
    sybyl_table.SetFromType("INT");
    sybyl_table.SetToType("SYB");

    [[maybe_unused]] unsigned long max_atom_index{0};
    for (auto atom_it = source->BeginAtoms(); atom_it < source->EndAtoms(); ++atom_it) {
      const auto atom      = *atom_it;
      const auto atom_id   = atom->GetId();
      const auto atom_type = atom->GetType();
      std::string ob_type  = atom->GetType();

      const auto atom_element = parse_element_symbol(ttab.Translate(atom_type));
      const auto sybyl_type   = parse_sybyl_atom_type(sybyl_table.Translate(ob_type));

      dest.atom_name(mudock_atom_index)  = get_ob_atom_name(atom);
      dest.sybyl_type(mudock_atom_index) = sybyl_type;

      dest.residue_id(mudock_atom_index)        = get_ob_residue_id(atom);
      dest.residue_name(mudock_atom_index)      = get_ob_residue_name(atom);
      dest.atom_residue_type(mudock_atom_index) = parse_residue_type(dest.residue_name(mudock_atom_index));

      // Get which atoms are aromatic
      dest.elements(mudock_atom_index) = atom_element;
      dest.x(mudock_atom_index)        = static_cast<fp_type>(atom->GetX());
      dest.y(mudock_atom_index)        = static_cast<fp_type>(atom->GetY());
      dest.z(mudock_atom_index)        = static_cast<fp_type>(atom->GetZ());
      // if constexpr (std::derived_from<molecule_type, autodock_static_layer> ||
      //               std::derived_from<molecule_type, autodock_dynamic_layer>) {
      dest.charge(mudock_atom_index)      = static_cast<fp_type>(atom->GetPartialCharge());
      dest.is_aromatic(mudock_atom_index) = atom->IsAromatic();
      // }

			// parsing residue for proteins
      if constexpr (std::same_as<std::remove_cvref_t<molecule_type>, dynamic_molecule>) {
        if (OpenBabel::OBResidue* ob_res = atom->GetResidue()) {
          std::string res_name = ob_res->GetName();
          std::string atom_name = ob_res->GetAtomID(atom);
          dest.residue_types(mudock_atom_index) = parse_residue_name(res_name);
          // remove spaces from atom name, e.g. " CA " -> "CA"
          atom_name.erase(std::remove(atom_name.begin(), atom_name.end(), ' '), atom_name.end());
          dest.atom_name(mudock_atom_index) = atom_name;
        } else {
          dest.residue_types(mudock_atom_index) = residue::UNKNOWN;
          dest.atom_name(mudock_atom_index) = "UNKNOWN";
        }
      }

      index_translator.emplace(atom_id, mudock_atom_index);
      ++mudock_atom_index;
      max_atom_index = std::max(max_atom_index, atom_id);
    }
    // Verify that there are no gaps in molecule indexes
    assert(((max_atom_index + 1) == num_atoms) && (static_cast<int>(num_atoms) == mudock_atom_index));

    // fill the bond information
    int mudock_bond_index{0};
    for (auto bond_it = source->BeginBonds(); bond_it < source->EndBonds(); ++bond_it) {
      const auto bond                  = *bond_it;
      const std::size_t atom_id_source = bond->GetBeginAtomIdx() - 1;
      const std::size_t atom_id_dest   = bond->GetEndAtomIdx() - 1;
      auto& mudock_bond                = dest.bonds(mudock_bond_index);
      mudock_bond.source     = static_cast<int>(index_translator.at(static_cast<int>(atom_id_source)));
      mudock_bond.dest       = static_cast<int>(index_translator.at(static_cast<int>(atom_id_dest)));
      mudock_bond.type       = parse_ob_bond_type(*bond);
      mudock_bond.can_rotate = check_rotor_bond(*bond);
      ++mudock_bond_index;
    }

    // store the molecule name
    auto name = source->GetTitle();
    dest.properties.assign(property_type::NAME, name);
  }

  static inline int element_to_atomic_num(const element e) {
    // Replace with your actual mapping
    return static_cast<int>(e); // placeholder if your enum already matches Z
  }

  // map your bond type -> OpenBabel bond order and aromatic flag
  static inline std::pair<int, bool> to_ob_bond_order_and_aromatic(bond_type t) {
    // Adjust to your enum; common choices shown:
    switch (t) {
      case bond_type::SINGLE: return {1, false};
      case bond_type::DOUBLE: return {2, false};
      case bond_type::TRIPLE: return {3, false};
      case bond_type::AROMATIC: return {5, true}; // OB uses 5 for aromatic
      default: return {1, false};
    }
  }

  template<class molecule_type>
    requires derived_from_molecule<molecule_type>
  void convert(ob_mol_wrapper& dest, const molecule_type& src) {
    OpenBabel::OBMol& mol = *dest;
    mol.Clear(); // start fresh if reusing wrapper
    mol.BeginModify();

    const size_t natoms = src.num_atoms();
    const size_t nbonds = src.num_bonds();

    // --- Atoms ---
    std::vector<OpenBabel::OBAtom*> ob_atoms;
    ob_atoms.reserve(natoms);

    for (int i = 0; i < static_cast<int>(natoms); ++i) {
      OpenBabel::OBAtom* a = mol.NewAtom(); // creates atom with new index (1-based)
      a->SetAtomicNum(element_to_atomic_num(src.elements(i)));
      a->SetVector(src.x(i), src.y(i), src.z(i));
      a->SetPartialCharge(src.charge(i));
      // Atom aromaticity (mostly derived from bonds in OB, but setting doesn’t hurt)
      a->SetAromatic(src.is_aromatic(i));
      a->SetId(static_cast<unsigned int>(i)); // keep our 0-based id if you use it elsewhere
      ob_atoms.push_back(a);
    }

    // --- Bonds ---
    for (size_t k = 0; k < nbonds; ++k) {
      const auto& b = src.bonds(static_cast<int>(k));
      // Your indices appear 0-based; OB API expects 1-based atom indices
      const int a1 = static_cast<int>(b.source) + 1;
      const int a2 = static_cast<int>(b.dest) + 1;

      const auto [order, arom] = to_ob_bond_order_and_aromatic(b.type);
      // AddBond(beginIdx, endIdx, order)
      mol.AddBond(a1, a2, order);

      // Optionally flag aromatic explicitly
      if (arom) {
        if (OpenBabel::OBBond* obb = mol.GetBond(mol.NumBonds() - 1)) {
          obb->SetAromatic(true);
        }
      }
    }
    // --- Title / name ---
    mol.SetTitle(src.properties.get(property_type::NAME).c_str());

    // If your coordinates are already in Å, leave as is.
    // If you need explicit dimensionality/flags, set them here:
    // mol.SetDimension(3);

    // OpenBabel often (re)perceives aromaticity/valence; do it if you want OB’s perception:
    // mol.PerceiveBondOrders();
    // mol.PerceiveAromaticity(); // requires kekulized structure/bond orders

    mol.EndModify();

    // If you rely on OB’s internal properties (e.g., implicit hydrogens), consider:
    // OBMolAtomIter ai(mol); for (OBAtom* a = ai.Begin(); a; a = ai.Next()) a->ImplicitHydrogenCount(); // etc.
  }
} // namespace mudock
