#pragma once

#include <array>
#include <cassert>
#include <cstdint>
#include <mudock/chem/elements.hpp>
#include <optional>
#include <string_view>

//===------------------------------------------------------------------------------------------------------
// WARNING: This file has been automatically generated from chem/autodock_babel_types.json
//===------------------------------------------------------------------------------------------------------

namespace mudock {

  // this is the list of all the known atoms
  enum class autodock_babel_ff : int {
    H       = 0,   // H
    HC      = 1,   // HC
    He      = 2,   // He
    Li      = 3,   // Li
    Be      = 4,   // Be
    B       = 5,   // B
    C       = 6,   // C
    C1      = 7,   // C1
    C2      = 8,   // C2
    C3      = 9,   // C3
    Cac     = 10,  // Cac
    C_plus  = 11,  // C+
    N       = 12,  // N
    Nox     = 13,  // Nox
    N1      = 14,  // N1
    N2      = 15,  // N2
    N3      = 16,  // N3
    N3_plus = 17,  // N3+
    Ntr     = 18,  // Ntr
    Npl     = 19,  // Npl
    Ng_plus = 20,  // Ng+
    Nam     = 21,  // Nam
    O       = 22,  // O
    O_minus = 23,  // O-
    O3      = 24,  // O3
    O2      = 25,  // O2
    F       = 26,  // F
    Ne      = 27,  // Ne
    Na      = 28,  // Na
    Mg      = 29,  // Mg
    Al      = 30,  // Al
    Si      = 31,  // Si
    P       = 32,  // P
    Pac     = 33,  // Pac
    Pox     = 34,  // Pox
    P3_plus = 35,  // P3+
    S       = 36,  // S
    S2      = 37,  // S2
    S3      = 38,  // S3
    S3_plus = 39,  // S3+
    Sac     = 40,  // Sac
    Sox     = 41,  // Sox
    Cl      = 42,  // Cl
    Ar      = 43,  // Ar
    K       = 44,  // K
    Ca      = 45,  // Ca
    Sc      = 46,  // Sc
    Ti      = 47,  // Ti
    V       = 48,  // V
    Cr      = 49,  // Cr
    Mn      = 50,  // Mn
    Fe      = 51,  // Fe
    Co      = 52,  // Co
    Ni      = 53,  // Ni
    Cu      = 54,  // Cu
    Zn      = 55,  // Zn
    Ga      = 56,  // Ga
    Ge      = 57,  // Ge
    As      = 58,  // As
    Se      = 59,  // Se
    Br      = 60,  // Br
    Kr      = 61,  // Kr
    Rb      = 62,  // Rb
    Sr      = 63,  // Sr
    Y       = 64,  // Y
    Zr      = 65,  // Zr
    Nb      = 66,  // Nb
    Mo      = 67,  // Mo
    Tc      = 68,  // Tc
    Ru      = 69,  // Ru
    Rh      = 70,  // Rh
    Pd      = 71,  // Pd
    Ag      = 72,  // Ag
    Cd      = 73,  // Cd
    In      = 74,  // In
    Sn      = 75,  // Sn
    Sb      = 76,  // Sb
    Te      = 77,  // Te
    I       = 78,  // I
    Xe      = 79,  // Xe
    Cs      = 80,  // Cs
    Ba      = 81,  // Ba
    La      = 82,  // La
    Ce      = 83,  // Ce
    Pr      = 84,  // Pr
    Nd      = 85,  // Nd
    Pm      = 86,  // Pm
    Sm      = 87,  // Sm
    Eu      = 88,  // Eu
    Gd      = 89,  // Gd
    Tb      = 90,  // Tb
    Dy      = 91,  // Dy
    Ho      = 92,  // Ho
    Er      = 93,  // Er
    Tm      = 94,  // Tm
    Yb      = 95,  // Yb
    Lu      = 96,  // Lu
    Hf      = 97,  // Hf
    Ta      = 98,  // Ta
    W       = 99,  // W
    Re      = 100, // Re
    Os      = 101, // Os
    Ir      = 102, // Ir
    Pt      = 103, // Pt
    Au      = 104, // Au
    Hg      = 105, // Hg
    Tl      = 106, // Tl
    Pb      = 107, // Pb
    Bi      = 108, // Bi
    Po      = 109, // Po
    At      = 110, // At
    Rn      = 111, // Rn
    Fr      = 112, // Fr
    Ra      = 113, // Ra
    Ac      = 114, // Ac
    Th      = 115, // Th
    Pa      = 116, // Pa
    U       = 117, // U
    Np      = 118, // Np
    Pu      = 119, // Pu
    Am      = 120, // Am
    Cm      = 121, // Cm
    Bk      = 122, // Bk
    Cf      = 123, // Cf
    Es      = 124, // Es
    Fm      = 125, // Fm
    Md      = 126, // Md
    No      = 127, // No
    Lr      = 128, // Lr
    Rf      = 129, // Rf
    Db      = 130, // Db
    Sg      = 131, // Sg
    Bh      = 132, // Bh
    Hs      = 133, // Hs
    Mt      = 134, // Mt
    Ds      = 135, // Ds
    Rg      = 136, // Rg
    Cn      = 137, // Cn
    Nh      = 138, // Nh
    Fl      = 139, // Fl
    Mc      = 140, // Mc
    Lv      = 141, // Lv
    Ts      = 142, // Ts
    Og      = 143, // Og
    Uue     = 144, // Uue
  };

  // this is knowledge that we have about all the elements (that we need at least)
  struct autodock_babel_ff_description {
    autodock_babel_ff value;
    std::string_view name;
    element base_element;
  };
  extern const std::array<autodock_babel_ff_description, 145> AUTODOCK_BABEL_FF_DICTIONARY;

  // utility functions to work with them
  inline const autodock_babel_ff_description& get_description(const autodock_babel_ff a) {
    assert(AUTODOCK_BABEL_FF_DICTIONARY[static_cast<int>(a)].value == a);
    return AUTODOCK_BABEL_FF_DICTIONARY[static_cast<int>(a)];
  }

} // namespace mudock
