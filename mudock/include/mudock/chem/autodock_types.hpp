#pragma once

#include <array>
#include <cassert>
#include <cstdint>
#include <mudock/chem/elements.hpp>
#include <mudock/type_alias.hpp>
#include <optional>
#include <string_view>

//===------------------------------------------------------------------------------------------------------
// WARNING: This file has been automatically generated from chem/autodock_babel_types.json
//===------------------------------------------------------------------------------------------------------

namespace mudock {

  // this is the list of all the known atoms
  enum class autodock_ff : std::size_t {
    H   = 0,   // H
    HD  = 1,   // HD
    HS  = 2,   // HS
    He  = 3,   // He
    Li  = 4,   // Li
    Be  = 5,   // Be
    B   = 6,   // B
    C   = 7,   // C
    A   = 8,   // A
    N   = 9,   // N
    NA  = 10,  // NA
    NS  = 11,  // NS
    OA  = 12,  // OA
    OS  = 13,  // OS
    F   = 14,  // F
    Ne  = 15,  // Ne
    Na  = 16,  // Na
    Mg  = 17,  // Mg
    Al  = 18,  // Al
    Si  = 19,  // Si
    P   = 20,  // P
    SA  = 21,  // SA
    S   = 22,  // S
    Cl  = 23,  // Cl
    Ar  = 24,  // Ar
    K   = 25,  // K
    Ca  = 26,  // Ca
    Sc  = 27,  // Sc
    Ti  = 28,  // Ti
    V   = 29,  // V
    Cr  = 30,  // Cr
    Mn  = 31,  // Mn
    Fe  = 32,  // Fe
    Co  = 33,  // Co
    Ni  = 34,  // Ni
    Cu  = 35,  // Cu
    Zn  = 36,  // Zn
    Ga  = 37,  // Ga
    Ge  = 38,  // Ge
    As  = 39,  // As
    Se  = 40,  // Se
    Br  = 41,  // Br
    Kr  = 42,  // Kr
    Rb  = 43,  // Rb
    Sr  = 44,  // Sr
    Y   = 45,  // Y
    Zr  = 46,  // Zr
    Nb  = 47,  // Nb
    Mo  = 48,  // Mo
    Tc  = 49,  // Tc
    Ru  = 50,  // Ru
    Rh  = 51,  // Rh
    Pd  = 52,  // Pd
    Ag  = 53,  // Ag
    Cd  = 54,  // Cd
    In  = 55,  // In
    Sn  = 56,  // Sn
    Sb  = 57,  // Sb
    Te  = 58,  // Te
    I   = 59,  // I
    Xe  = 60,  // Xe
    Cs  = 61,  // Cs
    Ba  = 62,  // Ba
    La  = 63,  // La
    Ce  = 64,  // Ce
    Pr  = 65,  // Pr
    Nd  = 66,  // Nd
    Pm  = 67,  // Pm
    Sm  = 68,  // Sm
    Eu  = 69,  // Eu
    Gd  = 70,  // Gd
    Tb  = 71,  // Tb
    Dy  = 72,  // Dy
    Ho  = 73,  // Ho
    Er  = 74,  // Er
    Tm  = 75,  // Tm
    Yb  = 76,  // Yb
    Lu  = 77,  // Lu
    Hf  = 78,  // Hf
    Ta  = 79,  // Ta
    W   = 80,  // W
    Re  = 81,  // Re
    Os  = 82,  // Os
    Ir  = 83,  // Ir
    Pt  = 84,  // Pt
    Au  = 85,  // Au
    Hg  = 86,  // Hg
    Tl  = 87,  // Tl
    Pb  = 88,  // Pb
    Bi  = 89,  // Bi
    Po  = 90,  // Po
    At  = 91,  // At
    Rn  = 92,  // Rn
    Fr  = 93,  // Fr
    Ra  = 94,  // Ra
    Ac  = 95,  // Ac
    Th  = 96,  // Th
    Pa  = 97,  // Pa
    U   = 98,  // U
    Np  = 99,  // Np
    Pu  = 100, // Pu
    Am  = 101, // Am
    Cm  = 102, // Cm
    Bk  = 103, // Bk
    Cf  = 104, // Cf
    Es  = 105, // Es
    Fm  = 106, // Fm
    Md  = 107, // Md
    No  = 108, // No
    Lr  = 109, // Lr
    Rf  = 110, // Rf
    Db  = 111, // Db
    Sg  = 112, // Sg
    Bh  = 113, // Bh
    Hs  = 114, // Hs
    Mt  = 115, // Mt
    Ds  = 116, // Ds
    Rg  = 117, // Rg
    Cn  = 118, // Cn
    Nh  = 119, // Nh
    Fl  = 120, // Fl
    Mc  = 121, // Mc
    Lv  = 122, // Lv
    Ts  = 123, // Ts
    Og  = 124, // Og
    Uue = 125, // Uue
    Z   = 126, // Z
    G   = 127, // G
    GA  = 128, // GA
    J   = 129, // J
    Q   = 130, // Q
  };

  // this is knowledge that we have about all the elements (that we need at least)
  struct autodock_ff_description {
    autodock_ff value;
    std::string_view name;
    coordinate_type Rii      = 0;
    coordinate_type epsii    = 0;
    coordinate_type vol      = 0;
    coordinate_type solpar   = 0;
    coordinate_type Rij_hb   = 0;
    coordinate_type epsij_hb = 0;
    std::size_t hbond        = 0;
  };
  extern const std::array<autodock_ff_description, 131> AUTODOCK_FF_DICTIONARY;

  // utility functions to work with them
  inline const autodock_ff_description& get_description(const autodock_ff a) {
    assert(AUTODOCK_FF_DICTIONARY[static_cast<std::size_t>(a)].value == a);
    return AUTODOCK_FF_DICTIONARY[static_cast<std::size_t>(a)];
  }
} // namespace mudock
