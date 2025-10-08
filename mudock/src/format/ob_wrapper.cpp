#include "mudock/format/supported_format.hpp"
#include "mudock/molecule.hpp"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <memory>
#include <mudock/format/ob_wrapper.hpp>
#include <mudock/utils.hpp>
#include <openbabel/atom.h>
#include <openbabel/babelconfig.h>
#include <openbabel/chargemodel.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>
#include <openbabel/obutil.h>
#include <openbabel/plugin.h>
#include <stdexcept>
#include <string_view>

namespace mudock {

  void writer(const ob_mol_wrapper& mol, const std::filesystem::path out_path) {
    const auto format = parse_supported_format(out_path);

    constexpr_switch<0, get_num_supported_format(), 1>(
        [&](const auto format_index) {
          format_writer<static_cast<supported_format>(format_index())>(mol, out_path);
        },
        format);
  }

  template<supported_format format>
  ob_mol_wrapper ob_parser(const std::string_view description) {
    OpenBabel::obErrorLog.SetOutputLevel(OpenBabel::obMessageLevel::obError); // Silence everything
    OpenBabel::OBConversion conv;
    const std::string ext{parse_supported_format(format)};
    conv.SetInFormat(ext.c_str());

    std::istringstream desc{std::string(description)};
    auto mol = std::make_unique<OpenBabel::OBMol>();
    if (!conv.Read(mol.get(), &desc)) {
      mudock::error(std::format("Couldn't open {} file", ext));
      throw std::runtime_error(std::format("{} Parser failed, look to logs for details", ext));
    }
    return mol;
  }

  template<>
  ob_mol_wrapper format_parser<supported_format::PDBQT>(const std::string_view description) {
    auto mol = ob_parser<supported_format::PDBQT>(description);

    return mol;
  }

  template<>
  ob_mol_wrapper format_parser<supported_format::PDB>(const std::string_view description) {
    auto mol = ob_parser<supported_format::PDB>(description);

    return mol;
  }

  template<>
  ob_mol_wrapper format_parser<supported_format::MOL2>(const std::string_view description) {
    auto mol = ob_parser<supported_format::MOL2>(description);
    return mol;
  }
  template<>
  ob_mol_wrapper format_parser<supported_format::MOL2X>(const std::string_view description) {
    static_molecule mol;
    mol2x::parse(mol, description);

    std::ostringstream oss;
    mol2::print(mol, oss);

    std::string s = oss.str();
    std::string_view mol2_description{s};
    return format_parser<supported_format::MOL2>(mol2_description);
  }
  template<>
  void format_parser<supported_format::MOL2X, static_molecule>(static_molecule& mol,
                                                               const std::string_view description) {
    mol2x::parse(mol, description);
  }
  template<>
  void format_parser<supported_format::MOL2X, dynamic_molecule>(dynamic_molecule& mol,
                                                                const std::string_view description) {
    mol2x::parse(mol, description);
  }

  void format_writer(const ob_mol_wrapper& mol, const std::filesystem::path out_path) {
    const auto format = parse_supported_format(out_path);

    constexpr_switch<0, get_num_supported_format(), 1>(
        [&](const auto format_index) {
          format_writer<static_cast<supported_format>(format_index())>(mol, out_path);
        },
        format);
  }

  template<supported_format format>
  void ob_writer(const ob_mol_wrapper& mol, std::ofstream& ofs) {
    OpenBabel::OBConversion conv;
    const std::string ext{parse_supported_format(format)};
    conv.SetOutFormat(ext.c_str());

    if (!ofs) {
      throw std::runtime_error("Error: Cannot open file for writing!");
    }

    if (!conv.Write(mol.get(), &ofs)) {
      throw std::runtime_error("Error: Failed to write molecule!");
    }
  }

  template<>
  void format_writer<supported_format::PDBQT>(const ob_mol_wrapper& mol, std::ofstream& ofs) {
    // Add charges using the Gasteiger method
    OpenBabel::OBChargeModel* chargeModel = OpenBabel::OBChargeModel::FindType("gasteiger");
    if (!chargeModel) {
      mudock::error("Error: Unable to find charge model.");
      throw std::runtime_error("Error in OpenBabel transformations");
    }
    chargeModel->ComputeCharges(*mol.get());

    mol.get()->ConnectTheDots();
    mol.get()->PerceiveBondOrders();

    ob_writer<supported_format::PDBQT>(mol, ofs);
  }

  template<>
  void format_writer<supported_format::PDBQT>(const ob_mol_wrapper& mol,
                                              const std::filesystem::path out_path) {
    std::ofstream ofs(out_path, std::ios::out);
    format_writer<supported_format::PDBQT>(mol, ofs);
  }

  template<>
  void format_writer<supported_format::MOL2>(const ob_mol_wrapper& mol,
                                             const std::filesystem::path out_path) {
    std::ofstream ofs(out_path, std::ios::out);
    ob_writer<supported_format::MOL2>(mol, ofs);
  }
  template<>
  void format_writer<supported_format::MOL2>(const ob_mol_wrapper& mol, std::ofstream& ofs) {
    ob_writer<supported_format::MOL2>(mol, ofs);
  }

  template<>
  void format_writer<supported_format::MOL2X>(const ob_mol_wrapper& mol, std::ofstream& ofs) {
    static_molecule s_mol;
    convert<rotate_check>(s_mol, mol);
    mudock::apply_autodock_forcefield(s_mol);
    mol2x::print(s_mol, ofs);
  }

  template<>
  void format_writer<supported_format::MOL2X>(const ob_mol_wrapper& mol,
                                              const std::filesystem::path out_path) {
    std::ofstream ofs(out_path, std::ios::out | std::ios::app);
    format_writer<supported_format::MOL2X>(mol, ofs);
  }

  template<>
  void format_writer<supported_format::PDB>(const ob_mol_wrapper& mol, std::ofstream& ofs) {
    ob_writer<supported_format::PDB>(mol, ofs);
  }

  template<>
  void format_writer<supported_format::PDB>(const ob_mol_wrapper& mol, const std::filesystem::path out_path) {
    std::ofstream ofs(out_path, std::ios::out);
    ob_writer<supported_format::PDB>(mol, ofs);
  }

  bond_type parse_ob_bond_type(const OpenBabel::OBBond& bond) {
    if (bond.IsAromatic())
      return bond_type::AROMATIC;
    else if (const_cast<OpenBabel::OBBond&>(bond).IsAmide())
      return bond_type::AMIDE;
    else
      switch (bond.GetBondOrder()) {
        case 1: return bond_type::SINGLE;
        case 2: return bond_type::DOUBLE;
        case 3: return bond_type::TRIPLE;
        default: throw std::runtime_error("Unsopported OpenBabel bond type");
      }
  }

  bool rotate_check(OpenBabel::OBBond& bond) { return bond.IsRotor(); }
  bool pdbqt_rotate_check(OpenBabel::OBBond& bond) {
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
    if ((bond.GetBondOrder() != 1 || bond.IsAromatic() || bond.IsAmide() || IsAmidine(&bond) ||
         bond.IsInRing()) ||
        (((bond.GetBeginAtom())->GetExplicitDegree() == 1) ||
         ((bond.GetEndAtom())->GetExplicitDegree() == 1))) {
      return false;
    }
    return true;
  }
  [[nodiscard]] autodock_ff convert_ob_elements(OpenBabel::OBAtom& ob_atom) {
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
    if (ob_atom.GetAtomicNum() == OpenBabel::OBElements::Hydrogen) {
      return autodock_ff::HD;
    } else if ((ob_atom.GetAtomicNum() == OpenBabel::OBElements::Carbon) && (ob_atom.IsAromatic())) {
      return autodock_ff::A;
    } else if (ob_atom.GetAtomicNum() == OpenBabel::OBElements::Oxygen) {
      return autodock_ff::OA;
    } else if ((ob_atom.GetAtomicNum() == OpenBabel::OBElements::Nitrogen) && (ob_atom.IsHbondAcceptor())) {
      return autodock_ff::NA;
    } else if ((ob_atom.GetAtomicNum() == OpenBabel::OBElements::Sulfur) && (ob_atom.IsHbondAcceptor())) {
      return autodock_ff::SA;
    } else {
      return convert_obelem(ob_atom.GetAtomicNum());
    }
  }

} // namespace mudock
