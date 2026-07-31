#include <cassert>
#include <mudock/format/ob_wrapper.hpp>
#include <mudock/format/supported_format.hpp>
#include <mudock/molecule.hpp>
#include <mudock/utils.hpp>
#include <openbabel/atom.h>
#include <openbabel/babelconfig.h>
#include <openbabel/chargemodel.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>
#include <openbabel/obutil.h>
#include <openbabel/plugin.h>
#include <openbabel/data.h>
#include <openbabel/parsmart.h>
#include <stdexcept>

namespace mudock {
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

  bool ob_rotate_check(OpenBabel::OBBond& bond) { return bond.IsRotor(); }
  bool pdbqt_rotate_check(OpenBabel::OBBond& b) {
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
    if ((b.GetBondOrder() != 1 || b.IsAromatic() || b.IsAmide() || IsAmidine(&b) || b.IsInRing()) ||
        (((b.GetBeginAtom())->GetExplicitDegree() == 1) || ((b.GetEndAtom())->GetExplicitDegree() == 1))) {
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

} // namespace mudock
