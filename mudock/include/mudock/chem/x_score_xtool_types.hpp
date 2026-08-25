#pragma once

#include <array>
#include <cassert>
#include <mudock/chem/x_score_hb.hpp>
#include <mudock/type_alias.hpp>
#include <string_view>

//===------------------------------------------------------------------------------------------------------
// WARNING: This file has been automatically generated from chem/x_score_xtool_types.json and subsequently modified manually
//===------------------------------------------------------------------------------------------------------

namespace mudock {

  enum class sybyl_atom_type : int;

  enum class xtool_ff : int {

    C3 = 0,

    C3x = 1,

    C3un = 2,

    C2 = 3,

    C2x = 4,

    C2un = 5,

    Car = 6,

    Carx = 7,

    Carun = 8,

    C1 = 9,

    C1x = 10,

    C1un = 11,

    Ccat = 12,

    N3h = 13,

    N3 = 14,

    N3un = 15,

    Npl3h = 16,

    Npl3 = 17,

    Npl3un = 18,

    N2h = 19,

    N2 = 20,

    N2un = 21,

    Narh = 22,

    Nar = 23,

    Narun = 24,

    N1 = 25,

    N1un = 26,

    N4 = 27,

    O3h = 28,

    O3 = 29,

    O3un = 30,

    O2 = 31,

    O2un = 32,

    Oco2 = 33,

    S3h = 34,

    S3 = 35,

    S3un = 36,

    S2 = 37,

    S2un = 38,

    So = 39,

    P3 = 40,

    F = 41,

    Cl = 42,

    Br = 43,

    I = 44,

    H = 45,

    Hhb = 46,

    Si = 47,

    Ow = 48,

    Mplus = 49,

    Un = 50,

    // non xtool types added for type mismatch in residue (added manually)

    //xlogp
    Nam = 51,

    So2 = 52,

    //ions
    Li = 53,

    Na = 54,

    K = 55,

    Mg = 56,

    Ca = 57,

    Mn = 58,

    Fe = 59,

    Co = 60,

    Ni = 61,

    Cu = 62,

    Zn = 63,

    Cd = 64,

    Hg = 65,

    Al = 66,

    U = 67,

    Fminus = 68,

    Clminus = 69,

    Brminus = 70,

    Iminus = 71,

    };

    struct xtool_ff_description {
    xtool_ff value;
    std::string_view name;
    fp_type atomic_weight;
    fp_type vdw_radius;
    fp_type vdw_potential;
    fp_type par_charge;
    x_score_hb hbond;
    };
    static constexpr auto num_xtool_ff() { return 72; }
    extern const std::array<xtool_ff_description, num_xtool_ff()> XTOOL_FF_DICTIONARY;

  inline const xtool_ff_description& get_description(const xtool_ff a) {
    assert(XTOOL_FF_DICTIONARY[static_cast<int>(a)].value == a);
    return XTOOL_FF_DICTIONARY[static_cast<int>(a)];
  }

  xtool_ff parse_xtool_type(const std::string_view symbol);

  xtool_ff xtool_type_from_sybyl(sybyl_atom_type type);
} // namespace mudock