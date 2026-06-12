#include <algorithm>
#include <cctype>
#include <mudock/chem/sybyl_atom_types.hpp>
#include <string>

namespace mudock {

  namespace {

    constexpr element unknown_element() {
      // muDock currently has no UNKNOWN element. Du/Xx/Lp are mapped to H only as
      // a safe in-range placeholder. Callers should not use get_element() for
      // sybyl_atom_type::UNKNOWN/Du/Xx/Lp to infer real chemistry.
      return element::H;
    }

    std::string lower_copy(std::string_view input) {
      std::string output(input);
      std::transform(output.begin(), output.end(), output.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
      });
      return output;
    }

  } // namespace

  const std::array<sybyl_atom_type_description, num_sybyl_atom_types()> SYBYL_ATOM_TYPE_DICTIONARY = {{
      {sybyl_atom_type::UNKNOWN, "", unknown_element(), false},

      {sybyl_atom_type::Du, "Du", unknown_element(), false},
      {sybyl_atom_type::Xx, "Xx", unknown_element(), false},
      {sybyl_atom_type::Lp, "Lp", unknown_element(), false},
      {sybyl_atom_type::X, "X", unknown_element(), false},
      {sybyl_atom_type::Any, "Any", unknown_element(), false},
      {sybyl_atom_type::Hal, "Hal", unknown_element(), false},
      {sybyl_atom_type::Het, "Het", unknown_element(), false},
      {sybyl_atom_type::Hev, "Hev", unknown_element(), false},

      {sybyl_atom_type::H, "H", element::H, false},
      {sybyl_atom_type::H_spc, "H.spc", element::H, false},
      {sybyl_atom_type::H_t3p, "H.t3p", element::H, false},
      {sybyl_atom_type::D, "D", element::H, false},
      {sybyl_atom_type::T, "T", element::H, false},

      {sybyl_atom_type::B, "B", element::B, false},
      {sybyl_atom_type::B_2, "B.2", element::B, false},
      {sybyl_atom_type::B_3, "B.3", element::B, false},

      {sybyl_atom_type::C, "C", element::C, false},
      {sybyl_atom_type::C_0, "C.0", element::C, false},
      {sybyl_atom_type::C_1, "C.1", element::C, false},
      {sybyl_atom_type::C_2, "C.2", element::C, false},
      {sybyl_atom_type::C_3, "C.3", element::C, false},
      {sybyl_atom_type::C_ar, "C.ar", element::C, true},
      {sybyl_atom_type::C_cat, "C.cat", element::C, false},
      {sybyl_atom_type::C_catC, "C.catC", element::C, false},
      {sybyl_atom_type::C_dot, "C.", element::C, false},
      {sybyl_atom_type::R_dot, "R.", element::C, false},
      {sybyl_atom_type::R_ar, "R.ar", element::C, true},
      {sybyl_atom_type::Nr_dot, "Nr.", element::C, false},

      {sybyl_atom_type::N, "N", element::N, false},
      {sybyl_atom_type::N_0, "N.0", element::N, false},
      {sybyl_atom_type::N_1, "N.1", element::N, false},
      {sybyl_atom_type::N_2, "N.2", element::N, false},
      {sybyl_atom_type::N_3, "N.3", element::N, false},
      {sybyl_atom_type::N_4, "N.4", element::N, false},
      {sybyl_atom_type::N_ar, "N.ar", element::N, true},
      {sybyl_atom_type::N_am, "N.am", element::N, false},
      {sybyl_atom_type::N_pl3, "N.pl3", element::N, false},
      {sybyl_atom_type::N_pl3N, "N.pl3N", element::N, false},
      {sybyl_atom_type::N_dot, "N.", element::N, false},

      {sybyl_atom_type::O, "O", element::O, false},
      {sybyl_atom_type::O_0, "O.0", element::O, false},
      {sybyl_atom_type::O_1, "O.1", element::O, false},
      {sybyl_atom_type::O_2, "O.2", element::O, false},
      {sybyl_atom_type::O_3, "O.3", element::O, false},
      {sybyl_atom_type::O_co2, "O.co2", element::O, false},
      {sybyl_atom_type::O_co2O, "O.co2O", element::O, false},
      {sybyl_atom_type::O_spc, "O.spc", element::O, false},
      {sybyl_atom_type::O_t3p, "O.t3p", element::O, false},
      {sybyl_atom_type::O_R, "O.R", element::O, true},
      {sybyl_atom_type::O_dot, "O.", element::O, false},

      {sybyl_atom_type::P, "P", element::P, false},
      {sybyl_atom_type::P_3, "P.3", element::P, false},

      {sybyl_atom_type::S, "S", element::S, false},
      {sybyl_atom_type::S_0, "S.0", element::S, false},
      {sybyl_atom_type::S_2, "S.2", element::S, false},
      {sybyl_atom_type::S_3, "S.3", element::S, false},
      {sybyl_atom_type::S_o, "S.o", element::S, false},
      {sybyl_atom_type::S_o2, "S.o2", element::S, false},
      {sybyl_atom_type::S_O, "S.O", element::S, false},
      {sybyl_atom_type::S_O2, "S.O2", element::S, false},
      {sybyl_atom_type::S_dot, "S.", element::S, false},

      {sybyl_atom_type::F, "F", element::F, false},
      {sybyl_atom_type::Cl, "Cl", element::Cl, false},
      {sybyl_atom_type::Br, "Br", element::Br, false},
      {sybyl_atom_type::I, "I", element::I, false},

      {sybyl_atom_type::He, "He", element::He, false},
      {sybyl_atom_type::Li, "Li", element::Li, false},
      {sybyl_atom_type::Be, "Be", element::Be, false},
      {sybyl_atom_type::Ne, "Ne", element::Ne, false},
      {sybyl_atom_type::Na, "Na", element::Na, false},
      {sybyl_atom_type::Mg, "Mg", element::Mg, false},
      {sybyl_atom_type::Al, "Al", element::Al, false},
      {sybyl_atom_type::Si, "Si", element::Si, false},
      {sybyl_atom_type::Ar, "Ar", element::Ar, false},
      {sybyl_atom_type::K, "K", element::K, false},
      {sybyl_atom_type::Ca, "Ca", element::Ca, false},
      {sybyl_atom_type::Sc, "Sc", element::Sc, false},
      {sybyl_atom_type::Ti, "Ti", element::Ti, false},
      {sybyl_atom_type::V, "V", element::V, false},
      {sybyl_atom_type::Cr, "Cr", element::Cr, false},
      {sybyl_atom_type::Cr_oh, "Cr.oh", element::Cr, false},
      {sybyl_atom_type::Cr_th, "Cr.th", element::Cr, false},
      {sybyl_atom_type::Mn, "Mn", element::Mn, false},
      {sybyl_atom_type::Fe, "Fe", element::Fe, false},
      {sybyl_atom_type::Co, "Co", element::Co, false},
      {sybyl_atom_type::Co_oh, "Co.oh", element::Co, false},
      {sybyl_atom_type::Ni, "Ni", element::Ni, false},
      {sybyl_atom_type::Cu, "Cu", element::Cu, false},
      {sybyl_atom_type::Zn, "Zn", element::Zn, false},
      {sybyl_atom_type::Ga, "Ga", element::Ga, false},
      {sybyl_atom_type::Ge, "Ge", element::Ge, false},
      {sybyl_atom_type::As, "As", element::As, false},
      {sybyl_atom_type::Se, "Se", element::Se, false},
      {sybyl_atom_type::Kr, "Kr", element::Kr, false},
      {sybyl_atom_type::Rb, "Rb", element::Rb, false},
      {sybyl_atom_type::Sr, "Sr", element::Sr, false},
      {sybyl_atom_type::Y, "Y", element::Y, false},
      {sybyl_atom_type::Zr, "Zr", element::Zr, false},
      {sybyl_atom_type::Nb, "Nb", element::Nb, false},
      {sybyl_atom_type::Mo, "Mo", element::Mo, false},
      {sybyl_atom_type::Tc, "Tc", element::Tc, false},
      {sybyl_atom_type::Ru, "Ru", element::Ru, false},
      {sybyl_atom_type::Rh, "Rh", element::Rh, false},
      {sybyl_atom_type::Pd, "Pd", element::Pd, false},
      {sybyl_atom_type::Ag, "Ag", element::Ag, false},
      {sybyl_atom_type::Cd, "Cd", element::Cd, false},
      {sybyl_atom_type::In, "In", element::In, false},
      {sybyl_atom_type::Sn, "Sn", element::Sn, false},
      {sybyl_atom_type::Sb, "Sb", element::Sb, false},
      {sybyl_atom_type::Te, "Te", element::Te, false},
      {sybyl_atom_type::Xe, "Xe", element::Xe, false},
      {sybyl_atom_type::Cs, "Cs", element::Cs, false},
      {sybyl_atom_type::Ba, "Ba", element::Ba, false},
      {sybyl_atom_type::La, "La", element::La, false},
      {sybyl_atom_type::Ce, "Ce", element::Ce, false},
      {sybyl_atom_type::Pr, "Pr", element::Pr, false},
      {sybyl_atom_type::Nd, "Nd", element::Nd, false},
      {sybyl_atom_type::Pm, "Pm", element::Pm, false},
      {sybyl_atom_type::Sm, "Sm", element::Sm, false},
      {sybyl_atom_type::Eu, "Eu", element::Eu, false},
      {sybyl_atom_type::Gd, "Gd", element::Gd, false},
      {sybyl_atom_type::Tb, "Tb", element::Tb, false},
      {sybyl_atom_type::Dy, "Dy", element::Dy, false},
      {sybyl_atom_type::Ho, "Ho", element::Ho, false},
      {sybyl_atom_type::Er, "Er", element::Er, false},
      {sybyl_atom_type::Tm, "Tm", element::Tm, false},
      {sybyl_atom_type::Yb, "Yb", element::Yb, false},
      {sybyl_atom_type::Lu, "Lu", element::Lu, false},
      {sybyl_atom_type::Hf, "Hf", element::Hf, false},
      {sybyl_atom_type::Ta, "Ta", element::Ta, false},
      {sybyl_atom_type::W, "W", element::W, false},
      {sybyl_atom_type::Re, "Re", element::Re, false},
      {sybyl_atom_type::Os, "Os", element::Os, false},
      {sybyl_atom_type::Ir, "Ir", element::Ir, false},
      {sybyl_atom_type::Pt, "Pt", element::Pt, false},
      {sybyl_atom_type::Au, "Au", element::Au, false},
      {sybyl_atom_type::Hg, "Hg", element::Hg, false},
      {sybyl_atom_type::Tl, "Tl", element::Tl, false},
      {sybyl_atom_type::Pb, "Pb", element::Pb, false},
      {sybyl_atom_type::Bi, "Bi", element::Bi, false},
      {sybyl_atom_type::Po, "Po", element::Po, false},
      {sybyl_atom_type::At, "At", element::At, false},
      {sybyl_atom_type::Rn, "Rn", element::Rn, false},
      {sybyl_atom_type::Fr, "Fr", element::Fr, false},
      {sybyl_atom_type::Ra, "Ra", element::Ra, false},
      {sybyl_atom_type::Ac, "Ac", element::Ac, false},
      {sybyl_atom_type::Th, "Th", element::Th, false},
      {sybyl_atom_type::Pa, "Pa", element::Pa, false},
      {sybyl_atom_type::U, "U", element::U, false},
      {sybyl_atom_type::Np, "Np", element::Np, false},
      {sybyl_atom_type::Pu, "Pu", element::Pu, false},
      {sybyl_atom_type::Am, "Am", element::Am, false},
      {sybyl_atom_type::Cm, "Cm", element::Cm, false},
      {sybyl_atom_type::Bk, "Bk", element::Bk, false},
      {sybyl_atom_type::Cf, "Cf", element::Cf, false},
      {sybyl_atom_type::Es, "Es", element::Es, false},
      {sybyl_atom_type::Fm, "Fm", element::Fm, false},
      {sybyl_atom_type::Md, "Md", element::Md, false},
      {sybyl_atom_type::No, "No", element::No, false},
      {sybyl_atom_type::Lr, "Lr", element::Lr, false},
      {sybyl_atom_type::Rf, "Rf", element::Rf, false},
      {sybyl_atom_type::Db, "Db", element::Db, false},
      {sybyl_atom_type::Sg, "Sg", element::Sg, false},
      {sybyl_atom_type::Bh, "Bh", element::Bh, false},
      {sybyl_atom_type::Hs, "Hs", element::Hs, false},
      {sybyl_atom_type::Mt, "Mt", element::Mt, false},
      {sybyl_atom_type::Ds, "Ds", element::Ds, false},
      {sybyl_atom_type::Rg, "Rg", element::Rg, false},
      {sybyl_atom_type::Cn, "Cn", element::Cn, false},
      {sybyl_atom_type::Nh, "Nh", element::Nh, false},
      {sybyl_atom_type::Fl, "Fl", element::Fl, false},
      {sybyl_atom_type::Mc, "Mc", element::Mc, false},
      {sybyl_atom_type::Lv, "Lv", element::Lv, false},
      {sybyl_atom_type::Ts, "Ts", element::Ts, false},
      {sybyl_atom_type::Og, "Og", element::Og, false},
      {sybyl_atom_type::Uue, "Uue", element::Uue, false},
  }};

  sybyl_atom_type parse_sybyl_atom_type(const std::string_view symbol) {
    for (const auto& entry: SYBYL_ATOM_TYPE_DICTIONARY) {
      if (entry.name == symbol) {
        return entry.value;
      }
    }

    // Common aliases and case normalizations.
    const auto lower = lower_copy(symbol);

    if (lower == "lp")
      return sybyl_atom_type::Lp;
    if (lower == "du")
      return sybyl_atom_type::Du;
    if (lower == "xx")
      return sybyl_atom_type::Xx;

    if (lower == "c.ar")
      return sybyl_atom_type::C_ar;
    if (lower == "c.cat")
      return sybyl_atom_type::C_cat;
    if (lower == "n.ar")
      return sybyl_atom_type::N_ar;
    if (lower == "n.am")
      return sybyl_atom_type::N_am;
    if (lower == "n.pl3" || lower == "n.p13")
      return sybyl_atom_type::N_pl3;
    if (lower == "o.co2")
      return sybyl_atom_type::O_co2;
    if (lower == "o.spc")
      return sybyl_atom_type::O_spc;
    if (lower == "o.t3p")
      return sybyl_atom_type::O_t3p;
    if (lower == "s.o")
      return sybyl_atom_type::S_o;
    if (lower == "s.o2")
      return sybyl_atom_type::S_o2;
    if (lower == "h.spc")
      return sybyl_atom_type::H_spc;
    if (lower == "h.t3p")
      return sybyl_atom_type::H_t3p;

    if (symbol == "CL" || symbol == "cl")
      return sybyl_atom_type::Cl;
    if (symbol == "BR" || symbol == "br")
      return sybyl_atom_type::Br;
    if (symbol == "NA" || symbol == "na")
      return sybyl_atom_type::Na;
    if (symbol == "CA" || symbol == "ca")
      return sybyl_atom_type::Ca;

    return sybyl_atom_type::UNKNOWN;
  }

} // namespace mudock
