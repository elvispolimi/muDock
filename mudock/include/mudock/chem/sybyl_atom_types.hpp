#pragma once

#include <array>
#include <cassert>
#include <cstddef>
#include <mudock/chem/elements.hpp>
#include <string_view>

namespace mudock {

  // MOL2/SYBYL atom types seen in the Tripos/SYBYL family and in Open Babel's
  // SYB column. This is intended as a compact classifier. For exact MOL2
  // round-tripping, keep the original token too, because MOL2 atom_type is not a
  // perfectly closed vocabulary in real-world files.
  enum class sybyl_atom_type : int {
    UNKNOWN = 0,

    // Pseudo / generic
    Du,
    Xx,
    Lp,
    X,
    Any,
    Hal,
    Het,
    Hev,

    // Hydrogen / isotopes / water models
    H,
    H_spc,
    H_t3p,
    D,
    T,

    // Boron
    B,
    B_2,
    B_3,

    // Carbon
    C,
    C_0,
    C_1,
    C_2,
    C_3,
    C_ar,
    C_cat,
    C_catC,
    C_dot,
    R_dot,
    R_ar,
    Nr_dot,

    // Nitrogen
    N,
    N_0,
    N_1,
    N_2,
    N_3,
    N_4,
    N_ar,
    N_am,
    N_pl3,
    N_pl3N,
    N_dot,

    // Oxygen
    O,
    O_0,
    O_1,
    O_2,
    O_3,
    O_co2,
    O_co2O,
    O_spc,
    O_t3p,
    O_R,
    O_dot,

    // Phosphorus
    P,
    P_3,

    // Sulfur
    S,
    S_0,
    S_2,
    S_3,
    S_o,
    S_o2,
    S_O,
    S_O2,
    S_dot,

    // Halogens
    F,
    Cl,
    Br,
    I,

    // Element-symbol fallbacks. These appear in MOL2 files, and Open Babel's
    // SYB column contains many of them directly for metals/inorganics.
    He,
    Li,
    Be,
    Ne,
    Na,
    Mg,
    Al,
    Si,
    Ar,
    K,
    Ca,
    Sc,
    Ti,
    V,
    Cr,
    Cr_oh,
    Cr_th,
    Mn,
    Fe,
    Co,
    Co_oh,
    Ni,
    Cu,
    Zn,
    Ga,
    Ge,
    As,
    Se,
    Kr,
    Rb,
    Sr,
    Y,
    Zr,
    Nb,
    Mo,
    Tc,
    Ru,
    Rh,
    Pd,
    Ag,
    Cd,
    In,
    Sn,
    Sb,
    Te,
    Xe,
    Cs,
    Ba,
    La,
    Ce,
    Pr,
    Nd,
    Pm,
    Sm,
    Eu,
    Gd,
    Tb,
    Dy,
    Ho,
    Er,
    Tm,
    Yb,
    Lu,
    Hf,
    Ta,
    W,
    Re,
    Os,
    Ir,
    Pt,
    Au,
    Hg,
    Tl,
    Pb,
    Bi,
    Po,
    At,
    Rn,
    Fr,
    Ra,
    Ac,
    Th,
    Pa,
    U,
    Np,
    Pu,
    Am,
    Cm,
    Bk,
    Cf,
    Es,
    Fm,
    Md,
    No,
    Lr,
    Rf,
    Db,
    Sg,
    Bh,
    Hs,
    Mt,
    Ds,
    Rg,
    Cn,
    Nh,
    Fl,
    Mc,
    Lv,
    Ts,
    Og,
    Uue,

    _COUNT
  };

  struct sybyl_atom_type_description {
    sybyl_atom_type value;
    std::string_view name;
    element element_value;
    bool aromatic = false;
  };

  static constexpr auto num_sybyl_atom_types() { return static_cast<std::size_t>(sybyl_atom_type::_COUNT); }

  extern const std::array<sybyl_atom_type_description, num_sybyl_atom_types()> SYBYL_ATOM_TYPE_DICTIONARY;

  inline const sybyl_atom_type_description& get_description(const sybyl_atom_type a) {
    assert(SYBYL_ATOM_TYPE_DICTIONARY[static_cast<std::size_t>(a)].value == a);
    return SYBYL_ATOM_TYPE_DICTIONARY[static_cast<std::size_t>(a)];
  }

  sybyl_atom_type parse_sybyl_atom_type(std::string_view symbol);

  inline std::string_view to_string(const sybyl_atom_type type) { return get_description(type).name; }

  inline element get_element(const sybyl_atom_type type) { return get_description(type).element_value; }

  inline bool is_aromatic(const sybyl_atom_type type) { return get_description(type).aromatic; }

} // namespace mudock
