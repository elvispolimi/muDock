
#include "mudock/chem/drugscore/drugscore_layer.hpp"
#include "mudock/chem/sybyl_atom_types.hpp"

namespace mudock {

  template<class ds_layer>
    requires derived_from_drugscore_layer<ds_layer>
  void assign_drugscore_types_impl(ds_layer& mol) {
    const auto syb_el2drugscore = [&](const sybyl_atom_type& el) {
      switch (el) {
        case (sybyl_atom_type::Li): return ds_type::METAL;
        case (sybyl_atom_type::Be): return ds_type::METAL;
        case (sybyl_atom_type::C_2): return ds_type::C_2;
        case (sybyl_atom_type::C_3): return ds_type::C_3;
        case (sybyl_atom_type::C_ar): return ds_type::C_ar;
        case (sybyl_atom_type::C_1): return ds_type::C_1;
        case (sybyl_atom_type::C_cat): return ds_type::C_cat;
        case (sybyl_atom_type::N_3): return ds_type::N_3;
        case (sybyl_atom_type::N_pl3): return ds_type::N_pl3;
        case (sybyl_atom_type::N_2): return ds_type::N_2;
        case (sybyl_atom_type::N_ar): return ds_type::N_ar;
        case (sybyl_atom_type::N_am): return ds_type::N_am;
        case (sybyl_atom_type::N_1): return ds_type::N_1;
        case (sybyl_atom_type::N_4): return ds_type::N_4;
        case (sybyl_atom_type::O_3): return ds_type::O_3;
        case (sybyl_atom_type::O_2): return ds_type::O_2;
        case (sybyl_atom_type::O_co2): return ds_type::O_co2;
        case (sybyl_atom_type::F): return ds_type::F;
        case (sybyl_atom_type::Na): return ds_type::METAL;
        case (sybyl_atom_type::Mg): return ds_type::METAL;
        case (sybyl_atom_type::Al): return ds_type::METAL;
        case (sybyl_atom_type::P_3): return ds_type::P_3;
        case (sybyl_atom_type::S_3): return ds_type::S_3;
        case (sybyl_atom_type::S_2): return ds_type::S_2;
        case (sybyl_atom_type::S_o2): return ds_type::S_o2;
        case (sybyl_atom_type::S_o): return ds_type::S_o;
        case (sybyl_atom_type::Cl): return ds_type::Cl;
        case (sybyl_atom_type::K): return ds_type::METAL;
        case (sybyl_atom_type::Ca): return ds_type::METAL;
        case (sybyl_atom_type::Sc): return ds_type::METAL;
        case (sybyl_atom_type::Ti): return ds_type::METAL;
        case (sybyl_atom_type::V): return ds_type::METAL;
        case (sybyl_atom_type::Cr_oh): return ds_type::METAL;
        case (sybyl_atom_type::Cr_th): return ds_type::METAL;
        case (sybyl_atom_type::Mn): return ds_type::METAL;
        case (sybyl_atom_type::Fe): return ds_type::METAL;
        case (sybyl_atom_type::Co_oh): return ds_type::METAL;
        case (sybyl_atom_type::Ni): return ds_type::METAL;
        case (sybyl_atom_type::Cu): return ds_type::METAL;
        case (sybyl_atom_type::Zn): return ds_type::METAL;
        case (sybyl_atom_type::Ga): return ds_type::METAL;
        case (sybyl_atom_type::Ge): return ds_type::METAL;
        case (sybyl_atom_type::Br): return ds_type::Br;
        case (sybyl_atom_type::Rb): return ds_type::METAL;
        case (sybyl_atom_type::Sr): return ds_type::METAL;
        case (sybyl_atom_type::Y): return ds_type::METAL;
        case (sybyl_atom_type::Zr): return ds_type::METAL;
        case (sybyl_atom_type::Nb): return ds_type::METAL;
        case (sybyl_atom_type::Mo): return ds_type::METAL;
        case (sybyl_atom_type::Tc): return ds_type::METAL;
        case (sybyl_atom_type::Ru): return ds_type::METAL;
        case (sybyl_atom_type::Rh): return ds_type::METAL;
        case (sybyl_atom_type::Pd): return ds_type::METAL;
        case (sybyl_atom_type::Ag): return ds_type::METAL;
        case (sybyl_atom_type::Cd): return ds_type::METAL;
        case (sybyl_atom_type::In): return ds_type::METAL;
        case (sybyl_atom_type::Sn): return ds_type::METAL;
        case (sybyl_atom_type::Sb): return ds_type::METAL;
        case (sybyl_atom_type::I): return ds_type::I;
        case (sybyl_atom_type::Cs): return ds_type::METAL;
        case (sybyl_atom_type::Ba): return ds_type::METAL;
        case (sybyl_atom_type::La): return ds_type::METAL;
        case (sybyl_atom_type::Ce): return ds_type::METAL;
        case (sybyl_atom_type::Pr): return ds_type::METAL;
        case (sybyl_atom_type::Nd): return ds_type::METAL;
        case (sybyl_atom_type::Pm): return ds_type::METAL;
        case (sybyl_atom_type::Sm): return ds_type::METAL;
        case (sybyl_atom_type::Eu): return ds_type::METAL;
        case (sybyl_atom_type::Gd): return ds_type::METAL;
        case (sybyl_atom_type::Tb): return ds_type::METAL;
        case (sybyl_atom_type::Dy): return ds_type::METAL;
        case (sybyl_atom_type::Ho): return ds_type::METAL;
        case (sybyl_atom_type::Er): return ds_type::METAL;
        case (sybyl_atom_type::Tm): return ds_type::METAL;
        case (sybyl_atom_type::Yb): return ds_type::METAL;
        case (sybyl_atom_type::Lu): return ds_type::METAL;
        case (sybyl_atom_type::Hf): return ds_type::METAL;
        case (sybyl_atom_type::Ta): return ds_type::METAL;
        case (sybyl_atom_type::W): return ds_type::METAL;
        case (sybyl_atom_type::Re): return ds_type::METAL;
        case (sybyl_atom_type::Os): return ds_type::METAL;
        case (sybyl_atom_type::Ir): return ds_type::METAL;
        case (sybyl_atom_type::Pt): return ds_type::METAL;
        case (sybyl_atom_type::Au): return ds_type::METAL;
        case (sybyl_atom_type::Hg): return ds_type::METAL;
        case (sybyl_atom_type::Tl): return ds_type::METAL;
        case (sybyl_atom_type::Pb): return ds_type::METAL;
        case (sybyl_atom_type::Bi): return ds_type::METAL;
        case (sybyl_atom_type::Po): return ds_type::METAL;
        case (sybyl_atom_type::Fr): return ds_type::METAL;
        case (sybyl_atom_type::Ra): return ds_type::METAL;
        case (sybyl_atom_type::Ac): return ds_type::METAL;
        case (sybyl_atom_type::Th): return ds_type::METAL;
        case (sybyl_atom_type::Pa): return ds_type::METAL;
        case (sybyl_atom_type::U): return ds_type::METAL;
        case (sybyl_atom_type::Np): return ds_type::METAL;
        case (sybyl_atom_type::Pu): return ds_type::METAL;
        case (sybyl_atom_type::Am): return ds_type::METAL;
        case (sybyl_atom_type::Cm): return ds_type::METAL;
        case (sybyl_atom_type::Bk): return ds_type::METAL;
        case (sybyl_atom_type::Cf): return ds_type::METAL;
        case (sybyl_atom_type::Es): return ds_type::METAL;
        case (sybyl_atom_type::Fm): return ds_type::METAL;
        case (sybyl_atom_type::Md): return ds_type::METAL;
        case (sybyl_atom_type::No): return ds_type::METAL;
        case (sybyl_atom_type::Lr): return ds_type::METAL;
        case (sybyl_atom_type::Rf): return ds_type::METAL;
        case (sybyl_atom_type::Db): return ds_type::METAL;
        case (sybyl_atom_type::Sg): return ds_type::METAL;
        case (sybyl_atom_type::Bh): return ds_type::METAL;
        case (sybyl_atom_type::Hs): return ds_type::METAL;
        case (sybyl_atom_type::Mt): return ds_type::METAL;
        case (sybyl_atom_type::Ds): return ds_type::METAL;
        case (sybyl_atom_type::Rg): return ds_type::METAL;
        case (sybyl_atom_type::Cn): return ds_type::METAL;
        default: return ds_type::ND;
      };
    };

    const auto elements = mol.get_base_molecule().get_sybyl_type();
    auto types          = mol.drugscore_type();

    assert(elements.size() == mol.get_drugscore_type().size());

    for (size_t i = 0; i < elements.size(); ++i) { types[i] = syb_el2drugscore(elements[i]); }
  };

  template<>
  void drugscore_static_layer::assign_drugscore_types(std::function<void(drugscore_static_layer&)> f) {
    if (f)
      f(*this);
    else
      assign_drugscore_types_impl(*this);
  };

  template<>
  void drugscore_dynamic_layer::assign_drugscore_types(std::function<void(drugscore_dynamic_layer&)> f) {
    if (f)
      f(*this);
    else
      assign_drugscore_types_impl(*this);
  };
} // namespace mudock