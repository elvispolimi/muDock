#pragma once

#include <array>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <memory>
#include <mudock/chem.hpp>
#include <mudock/chem/autodock_types.hpp>
#include <mudock/chem/elements.hpp>
#include <mudock/format/pdbqt.hpp>
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
#include <stdexcept>
#include <string>
#include <string_view>
#include <sys/types.h>

namespace mudock {

  //===------------------------------------------------------------------------------------------------------
  // Useful type alias to work with OpenBabel data structures
  //===------------------------------------------------------------------------------------------------------

  using ob_mol_wrapper = std::unique_ptr<OpenBabel::OBMol>;

  enum class supported_format : int { MOL2 = 0, PDBQT, PDB };

  struct format_description {
    supported_format format;
    std::string_view extension;
  };

  static constexpr std::array<format_description, 3> FORMAT_EXTENSIONS = {
      {{supported_format::MOL2, "mol2"}, {supported_format::PDBQT, "pdbqt"}, {supported_format::PDB, "pdb"}}};

  [[nodiscard]] supported_format parse_supported_format(const std::string_view extension);
  [[nodiscard]] inline supported_format parse_supported_format(const std::filesystem::path path) {
    assert(path.has_extension());
    const auto extension = path.extension();
    assert(!extension.empty());
    return parse_supported_format(std::string_view(extension.string().substr(1)));
  }
  [[nodiscard]] std::string_view parse_supported_format(const supported_format format);

  [[nodiscard]] ob_mol_wrapper parser(const std::filesystem::path file_path);
  void writer(const ob_mol_wrapper& mol, const std::filesystem::path out_path);

  //===------------------------------------------------------------------------------------------------------
  // Utility functions to parse convert simple data from OpenBabel to our simple data structure
  //===------------------------------------------------------------------------------------------------------

  // the list of functions that we use to parse a molecule
  template<supported_format format>
  [[nodiscard]] ob_mol_wrapper format_parser(const std::string_view description);

  template<>
  ob_mol_wrapper format_parser<supported_format::PDBQT>(const std::string_view description);
  template<>
  ob_mol_wrapper format_parser<supported_format::MOL2>(const std::string_view description);
  template<>
  ob_mol_wrapper format_parser<supported_format::PDB>(const std::string_view description);

  template<supported_format format>
  void format_writer(const ob_mol_wrapper& mol, const std::filesystem::path out_path);

  template<>
  void format_writer<supported_format::PDBQT>(const ob_mol_wrapper& mol,
                                              const std::filesystem::path out_path);
  template<>
  void format_writer<supported_format::MOL2>(const ob_mol_wrapper& mol, const std::filesystem::path out_path);
  template<>
  void format_writer<supported_format::PDB>(const ob_mol_wrapper& mol, const std::filesystem::path out_path);

  [[nodiscard]] bond_type parse_ob_bond_type(const OpenBabel::OBBond& bond_type);

  //===------------------------------------------------------------------------------------------------------
  // Translate an OpenBabel molecule to our internal format
  //===------------------------------------------------------------------------------------------------------

  template<class molecule_type>
    requires is_molecule<molecule_type>
  void convert(molecule_type&& dest, const ob_mol_wrapper& source) {
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
    dest.resize(num_atoms, num_bonds);

    // TODO check charges
    size_t mudock_atom_index{0};
    std::unordered_map<unsigned int, int> index_translator;
    OpenBabel::OBTypeTable ttab;
    ttab.SetFromType("INT");
    ttab.SetToType("XYZ");
    [[maybe_unused]] unsigned long max_atom_index{0};
    for (auto atom_it = source->BeginAtoms(); atom_it < source->EndAtoms(); ++atom_it) {
      const auto atom         = *atom_it;
      const auto atom_id      = atom->GetId();
      const auto atom_type    = atom->GetType();
      const auto atom_element = parse_element_symbol(ttab.Translate(atom_type));
      // Get which atoms are aromatic
      // TODO check the cast between bool and uinfast8_t
      dest.elements(mudock_atom_index) = atom_element;
      // TODO apparently is not working on PDBQT files, manually sets aromaticity
      // Ask Gadio how can I compute aromaticity!!
      dest.is_aromatic(mudock_atom_index) = atom->IsAromatic();
      dest.x(mudock_atom_index)           = static_cast<fp_type>(atom->GetX());
      dest.y(mudock_atom_index)           = static_cast<fp_type>(atom->GetY());
      dest.z(mudock_atom_index)           = static_cast<fp_type>(atom->GetZ());
      dest.charge(mudock_atom_index)      = atom->GetPartialCharge();
      index_translator.emplace(atom_id, mudock_atom_index);
      ++mudock_atom_index;
      max_atom_index = std::max(max_atom_index, atom_id);
    }
    // Verify that there are no gaps in molecule indexes
    assert(((max_atom_index + 1) == num_atoms) && (num_atoms == mudock_atom_index));

    // FIXME
    auto IsImide = [](OpenBabel::OBBond* querybond) {
      if (querybond->GetBondOrder() != 2)
        return (false);

      OpenBabel::OBAtom* bgn = querybond->GetBeginAtom();
      OpenBabel::OBAtom* end = querybond->GetEndAtom();
      if ((bgn->GetAtomicNum() == 6 && end->GetAtomicNum() == 7) ||
          (bgn->GetAtomicNum() == 7 && end->GetAtomicNum() == 6))
        return (true);

      return (false);
    };

    auto IsAmidine = [IsImide](OpenBabel::OBBond* querybond) {
      OpenBabel::OBAtom *c, *n;
      c = n = nullptr;

      // Look for C-N bond
      OpenBabel::OBAtom* bgn = querybond->GetBeginAtom();
      OpenBabel::OBAtom* end = querybond->GetEndAtom();
      if (bgn->GetAtomicNum() == 6 && end->GetAtomicNum() == 7) {
        c = bgn;
        n = end;
      }
      if (bgn->GetAtomicNum() == 7 && end->GetAtomicNum() == 6) {
        c = end;
        n = bgn;
      }
      if (!c || !n)
        return (false);
      if (querybond->GetBondOrder() != 1)
        return (false);
      if (n->GetTotalDegree() != 3)
        return false; // must be a degree 3 nitrogen

      // Make sure C is attached to =N
      OpenBabel::OBBond* bond;
      std::vector<OpenBabel::OBBond*>::iterator i;
      for (bond = c->BeginBond(i); bond; bond = c->NextBond(i)) {
        if (IsImide(bond))
          return (true);
      }

      // Return
      return (false);
    };

    // fill the bond information
    auto mudock_bond_index = int{0};
    for (auto bond_it = source->BeginBonds(); bond_it < source->EndBonds(); ++bond_it) {
      const auto bond          = *bond_it;
      const int atom_id_source = bond->GetBeginAtomIdx() - 1;
      const int atom_id_dest   = bond->GetEndAtomIdx() - 1;
      auto& mudock_bond        = dest.bonds(mudock_bond_index);
      mudock_bond.source       = index_translator.at(atom_id_source);
      mudock_bond.dest         = index_translator.at(atom_id_dest);
      mudock_bond.type         = parse_ob_bond_type(*bond);
      //mudock_bond.can_rotate   = bond->IsRotor(true);
      mudock_bond.can_rotate = true;
      if ((bond->GetBondOrder() != 1 || bond->IsAromatic() || bond->IsAmide() || IsAmidine(bond) ||
           bond->IsInRing()) ||
          (((bond->GetBeginAtom())->GetExplicitDegree() == 1) ||
           ((bond->GetEndAtom())->GetExplicitDegree() == 1))) {
        mudock_bond.can_rotate = false;
      }

      ++mudock_bond_index;
    }

    // store the molecule name
    auto name = source->GetTitle();
    dest.properties.assign(property_type::NAME, name);
  }

  // template<class molecule_type>
  //   requires is_molecule<molecule_type>
  // void convert(const ob_mol_wrapper& dest, molecule_type&& source) {
  //   // TODO
  //   throw std::runtime_error("Converted not implemented yet!");
  // }
  namespace openbabel {
    // The following function relies on the same order of atom loading and file
    template<class molecule_type>
      requires is_molecule<molecule_type>
    void apply_autodock_forcefield(molecule_type&& molecule, const std::filesystem::path input_path) {
      const auto format = parse_supported_format(input_path);
      assert(format == supported_format::PDBQT);

      const auto desc = read_from_stream(std::ifstream(input_path));
      std::stringstream desc_s{desc};

      const std::size_t num_atoms = molecule.num_atoms();

      std::size_t index = 0;
      std::string line;
      while (std::getline(desc_s, line)) {
        if (line.find(pdbqt::PDBQT_ATOM_TOKEN) != std::string::npos ||
            line.find(pdbqt::PDBQT_HETATOM_TOKEN) != std::string::npos) {
          assert(index < num_atoms);
          // FIXMED Really bad, at the moment we rely on OpenBabel PDBQT structure
          // What if PDBQT is standardized..
          if (line.size() < 79)
            line += " ";
          std::string adt_value = line.substr(77, 2);
          assert(adt_value.size() == 2);
          if (adt_value[1] == ' ')
            adt_value.pop_back();
          // FIX ME add check that the order of atoms is the same

          const auto adt                = parse_autodock_type(adt_value);
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
      }
    }
    template<class molecule_type>
      requires is_molecule<molecule_type>
    void apply_autodock_forcefield(molecule_type&& molecule, const ob_mol_wrapper& ob_mol) {
      const auto convert_obelem = [](const unsigned int ob_elem) {
        switch (ob_elem) {
          case (OpenBabel::OBElements::H):;
          case (OpenBabel::OBElements::He): return autodock_ff::He;
          case (OpenBabel::OBElements::Li): return autodock_ff::Li;
          case (OpenBabel::OBElements::Be): return autodock_ff::Be;
          case (OpenBabel::OBElements::B): return autodock_ff::B; // the original code is dead code
          case (OpenBabel::OBElements::C): return autodock_ff::C;
          case (OpenBabel::OBElements::N): return autodock_ff::N;
          case (OpenBabel::OBElements::F): return autodock_ff::F;
          case (OpenBabel::OBElements::Ne): return autodock_ff::Ne;
          case (OpenBabel::OBElements::Na): return autodock_ff::Na;
          case (OpenBabel::OBElements::Mg): return autodock_ff::Mg;
          case (OpenBabel::OBElements::Al): return autodock_ff::Al;
          case (OpenBabel::OBElements::Si): return autodock_ff::Si;
          case (OpenBabel::OBElements::P): return autodock_ff::P;
          case (OpenBabel::OBElements::S): return autodock_ff::S;
          case (OpenBabel::OBElements::Cl): return autodock_ff::Cl;
          case (OpenBabel::OBElements::Ar): return autodock_ff::Ar;
          case (OpenBabel::OBElements::K): return autodock_ff::K;
          case (OpenBabel::OBElements::Ca): return autodock_ff::Ca;
          case (OpenBabel::OBElements::Sc): return autodock_ff::Sc;
          case (OpenBabel::OBElements::Ti): return autodock_ff::Ti;
          case (OpenBabel::OBElements::V): return autodock_ff::V;
          case (OpenBabel::OBElements::Cr): return autodock_ff::Cr;
          case (OpenBabel::OBElements::Mn): return autodock_ff::Mn;
          case (OpenBabel::OBElements::Fe): return autodock_ff::Fe;
          case (OpenBabel::OBElements::Co): return autodock_ff::Co;
          case (OpenBabel::OBElements::Ni): return autodock_ff::Ni;
          case (OpenBabel::OBElements::Cu): return autodock_ff::Cu;
          case (OpenBabel::OBElements::Zn): return autodock_ff::Zn;
          case (OpenBabel::OBElements::Ga): return autodock_ff::Ga;
          case (OpenBabel::OBElements::Ge): return autodock_ff::Ge;
          case (OpenBabel::OBElements::As): return autodock_ff::As;
          case (OpenBabel::OBElements::Se): return autodock_ff::Se;
          case (OpenBabel::OBElements::Br): return autodock_ff::Br;
          case (OpenBabel::OBElements::Kr): return autodock_ff::Kr;
          case (OpenBabel::OBElements::Rb): return autodock_ff::Rb;
          case (OpenBabel::OBElements::Sr): return autodock_ff::Sr;
          case (OpenBabel::OBElements::Y): return autodock_ff::Y;
          case (OpenBabel::OBElements::Zr): return autodock_ff::Zr;
          case (OpenBabel::OBElements::Nb): return autodock_ff::Nb;
          case (OpenBabel::OBElements::Mo): return autodock_ff::Mo;
          case (OpenBabel::OBElements::Tc): return autodock_ff::Tc;
          case (OpenBabel::OBElements::Ru): return autodock_ff::Ru;
          case (OpenBabel::OBElements::Rh): return autodock_ff::Rh;
          case (OpenBabel::OBElements::Pd): return autodock_ff::Pd;
          case (OpenBabel::OBElements::Ag): return autodock_ff::Ag;
          case (OpenBabel::OBElements::Cd): return autodock_ff::Cd;
          case (OpenBabel::OBElements::In): return autodock_ff::In;
          case (OpenBabel::OBElements::Sn): return autodock_ff::Sn;
          case (OpenBabel::OBElements::Sb): return autodock_ff::Sb;
          case (OpenBabel::OBElements::Te): return autodock_ff::Te;
          case (OpenBabel::OBElements::I): return autodock_ff::I;
          case (OpenBabel::OBElements::Xe): return autodock_ff::Xe;
          case (OpenBabel::OBElements::Cs): return autodock_ff::Cs;
          case (OpenBabel::OBElements::Ba): return autodock_ff::Ba;
          case (OpenBabel::OBElements::La): return autodock_ff::La;
          case (OpenBabel::OBElements::Ce): return autodock_ff::Ce;
          case (OpenBabel::OBElements::Pr): return autodock_ff::Pr;
          case (OpenBabel::OBElements::Nd): return autodock_ff::Nd;
          case (OpenBabel::OBElements::Pm): return autodock_ff::Pm;
          case (OpenBabel::OBElements::Sm): return autodock_ff::Sm;
          case (OpenBabel::OBElements::Eu): return autodock_ff::Eu;
          case (OpenBabel::OBElements::Gd): return autodock_ff::Gd;
          case (OpenBabel::OBElements::Tb): return autodock_ff::Tb;
          case (OpenBabel::OBElements::Dy): return autodock_ff::Dy;
          case (OpenBabel::OBElements::Ho): return autodock_ff::Ho;
          case (OpenBabel::OBElements::Er): return autodock_ff::Er;
          case (OpenBabel::OBElements::Tm): return autodock_ff::Tm;
          case (OpenBabel::OBElements::Yb): return autodock_ff::Yb;
          case (OpenBabel::OBElements::Lu): return autodock_ff::Lu;
          case (OpenBabel::OBElements::Hf): return autodock_ff::Hf;
          case (OpenBabel::OBElements::Ta): return autodock_ff::Ta;
          case (OpenBabel::OBElements::W): return autodock_ff::W;
          case (OpenBabel::OBElements::Re): return autodock_ff::Re;
          case (OpenBabel::OBElements::Os): return autodock_ff::Os;
          case (OpenBabel::OBElements::Ir): return autodock_ff::Ir;
          case (OpenBabel::OBElements::Pt): return autodock_ff::Pt;
          case (OpenBabel::OBElements::Au): return autodock_ff::Au;
          case (OpenBabel::OBElements::Hg): return autodock_ff::Hg;
          case (OpenBabel::OBElements::Tl): return autodock_ff::Tl;
          case (OpenBabel::OBElements::Pb): return autodock_ff::Pb;
          case (OpenBabel::OBElements::Bi): return autodock_ff::Bi;
          case (OpenBabel::OBElements::Po): return autodock_ff::Po;
          case (OpenBabel::OBElements::At): return autodock_ff::At;
          case (OpenBabel::OBElements::Rn): return autodock_ff::Rn;
          case (OpenBabel::OBElements::Fr): return autodock_ff::Fr;
          case (OpenBabel::OBElements::Ra): return autodock_ff::Ra;
          case (OpenBabel::OBElements::Ac): return autodock_ff::Ac;
          case (OpenBabel::OBElements::Th): return autodock_ff::Th;
          case (OpenBabel::OBElements::Pa): return autodock_ff::Pa;
          case (OpenBabel::OBElements::U): return autodock_ff::U;
          case (OpenBabel::OBElements::Np): return autodock_ff::Np;
          case (OpenBabel::OBElements::Pu): return autodock_ff::Pu;
          case (OpenBabel::OBElements::Am): return autodock_ff::Am;
          case (OpenBabel::OBElements::Cm): return autodock_ff::Cm;
          case (OpenBabel::OBElements::Bk): return autodock_ff::Bk;
          case (OpenBabel::OBElements::Cf): return autodock_ff::Cf;
          case (OpenBabel::OBElements::Es): return autodock_ff::Es;
          case (OpenBabel::OBElements::Fm): return autodock_ff::Fm;
          case (OpenBabel::OBElements::Md): return autodock_ff::Md;
          case (OpenBabel::OBElements::No): return autodock_ff::No;
          case (OpenBabel::OBElements::Lr): return autodock_ff::Lr;
          case (OpenBabel::OBElements::Rf): return autodock_ff::Rf;
          case (OpenBabel::OBElements::Db): return autodock_ff::Db;
          case (OpenBabel::OBElements::Sg): return autodock_ff::Sg;
          case (OpenBabel::OBElements::Bh): return autodock_ff::Bh;
          case (OpenBabel::OBElements::Hs): return autodock_ff::Hs;
          case (OpenBabel::OBElements::Mt): return autodock_ff::Mt;
          case (OpenBabel::OBElements::Ds): return autodock_ff::Ds;
          case (OpenBabel::OBElements::Rg): return autodock_ff::Rg;
          case (OpenBabel::OBElements::Cn): return autodock_ff::Cn;
          case (OpenBabel::OBElements::Nh): return autodock_ff::Nh;
          case (OpenBabel::OBElements::Fl): return autodock_ff::Fl;
          case (OpenBabel::OBElements::Mc): return autodock_ff::Mc;
          case (OpenBabel::OBElements::Lv): return autodock_ff::Lv;
          case (OpenBabel::OBElements::Ts): return autodock_ff::Ts;
          case (OpenBabel::OBElements::Og): return autodock_ff::Og;
          default: throw std::runtime_error("Error OpenBabel typing of autodock, internal error");
        }
      };
      int index = 0;
      for (auto atom_it = ob_mol->BeginAtoms(); atom_it < ob_mol->EndAtoms(); ++atom_it) {
        const auto atom = *atom_it;

        autodock_ff adt;
        if (atom->GetAtomicNum() == OpenBabel::OBElements::Hydrogen) {
          adt = autodock_ff::HD;
        } else if ((atom->GetAtomicNum() == OpenBabel::OBElements::Carbon) && (atom->IsAromatic())) {
          adt = autodock_ff::A;
        } else if (atom->GetAtomicNum() == OpenBabel::OBElements::Oxygen) {
          adt = autodock_ff::OA;
        } else if ((atom->GetAtomicNum() == OpenBabel::OBElements::Nitrogen) && (atom->IsHbondAcceptor())) {
          adt = autodock_ff::NA;
        } else if ((atom->GetAtomicNum() == OpenBabel::OBElements::Sulfur) && (atom->IsHbondAcceptor())) {
          adt = autodock_ff::SA;
        } else {
          adt = convert_obelem(atom->GetAtomicNum());
        }
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
    }
  } // namespace openbabel

  template<class molecule_type>
    requires is_molecule<molecule_type>
  void parse(molecule_type&& molecule,
             const std::filesystem::path input_path,
             const bool assing_autodock_types = false) {
    const auto ob_mol = parser(input_path);
    convert(molecule, ob_mol);

    if (assing_autodock_types)
      openbabel::apply_autodock_forcefield(molecule, ob_mol);
  }

  template<supported_format format, class molecule_type>
    requires is_molecule<molecule_type>
  void parse(molecule_type&& molecule,
             const std::string_view description,
             const bool assing_autodock_types = false) {
    const auto ob_mol = format_parser<format>(description);
    convert(molecule, ob_mol);

    if (assing_autodock_types)
      openbabel::apply_autodock_forcefield(molecule, ob_mol);
  }

} // namespace mudock
