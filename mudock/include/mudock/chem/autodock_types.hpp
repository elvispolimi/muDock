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
  enum class autodock_ff : std::size_t {
    H = 0, // 
    HD = 1, // 
    He = 2, // 
    Li = 3, // 
    Be = 4, // 
    B = 5, // 
    C = 6, // 
    A = 7, // 
    N = 8, // 
    NA = 9, // 
    OA = 10, // 
    F = 11, // 
    Ne = 12, // 
    Na = 13, // 
    Mg = 14, // 
    Al = 15, // 
    Si = 16, // 
    P = 17, // 
    S = 18, // 
    SA = 19, // 
    Cl = 20, // 
    Ar = 21, // 
    K = 22, // 
    Ca = 23, // 
    Sc = 24, // 
    Ti = 25, // 
    V = 26, // 
    Cr = 27, // 
    Mn = 28, // 
    Fe = 29, // 
    Co = 30, // 
    Ni = 31, // 
    Cu = 32, // 
    Zn = 33, // 
    Ga = 34, // 
    Ge = 35, // 
    As = 36, // 
    Se = 37, // 
    Br = 38, // 
    Kr = 39, // 
    Rb = 40, // 
    Sr = 41, // 
    Y = 42, // 
    Zr = 43, // 
    Nb = 44, // 
    Mo = 45, // 
    Tc = 46, // 
    Ru = 47, // 
    Rh = 48, // 
    Pd = 49, // 
    Ag = 50, // 
    Cd = 51, // 
    In = 52, // 
    Sn = 53, // 
    Sb = 54, // 
    Te = 55, // 
    I = 56, // 
    Xe = 57, // 
    Cs = 58, // 
    Ba = 59, // 
    La = 60, // 
    Ce = 61, // 
    Pr = 62, // 
    Nd = 63, // 
    Pm = 64, // 
    Sm = 65, // 
    Eu = 66, // 
    Gd = 67, // 
    Tb = 68, // 
    Dy = 69, // 
    Ho = 70, // 
    Er = 71, // 
    Tm = 72, // 
    Yb = 73, // 
    Lu = 74, // 
    Hf = 75, // 
    Ta = 76, // 
    W = 77, // 
    Re = 78, // 
    Os = 79, // 
    Ir = 80, // 
    Pt = 81, // 
    Au = 82, // 
    Hg = 83, // 
    Tl = 84, // 
    Pb = 85, // 
    Bi = 86, // 
    Po = 87, // 
    At = 88, // 
    Rn = 89, // 
    Fr = 90, // 
    Ra = 91, // 
    Ac = 92, // 
    Th = 93, // 
    Pa = 94, // 
    U = 95, // 
    Np = 96, // 
    Pu = 97, // 
    Am = 98, // 
    Cm = 99, // 
    Bk = 100, // 
    Cf = 101, // 
    Es = 102, // 
    Fm = 103, // 
    Md = 104, // 
    No = 105, // 
    Lr = 106, // 
    Rf = 107, // 
    Db = 108, // 
    Sg = 109, // 
    Bh = 110, // 
    Hs = 111, // 
    Mt = 112, // 
    Ds = 113, // 
    Rg = 114, // 
    Cn = 115, // 
    Nh = 116, // 
    Fl = 117, // 
    Mc = 118, // 
    Lv = 119, // 
    Ts = 120, // 
    Og = 121, // 
    Uue = 122, // 
  };

  // this is knowledge that we have about all the elements (that we need at least)
  struct autodock_ff_description {
    autodock_ff value;
    std::string_view name;
  };
  extern const std::array<autodock_ff_description, 123> AUTODOCK_FF_DICTIONARY;
  
  // utility functions to work with them
  inline const autodock_ff_description& get_description(const autodock_ff a) {
    assert(AUTODOCK_FF_DICTIONARY[static_cast<std::size_t>(a)].value == a);
    return AUTODOCK_FF_DICTIONARY[static_cast<std::size_t>(a)];
  }
} // namespace mudock