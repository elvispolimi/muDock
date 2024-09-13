#pragma once

#include <array>
#include <cassert>
#include <cstdint>
#include <optional>
#include <string_view>

//===------------------------------------------------------------------------------------------------------
// WARNING: This file has been automatically generated from chem/periodic_table.json
//===------------------------------------------------------------------------------------------------------

namespace mudock {

  // this is the list of all the known atoms
  enum class element : int {
    H   = 0,   // Hydrogen
    He  = 1,   // Helium
    Li  = 2,   // Lithium
    Be  = 3,   // Beryllium
    B   = 4,   // Boron
    C   = 5,   // Carbon
    N   = 6,   // Nitrogen
    O   = 7,   // Oxygen
    F   = 8,   // Fluorine
    Ne  = 9,   // Neon
    Na  = 10,  // Sodium
    Mg  = 11,  // Magnesium
    Al  = 12,  // Aluminium
    Si  = 13,  // Silicon
    P   = 14,  // Phosphorus
    S   = 15,  // Sulfur
    Cl  = 16,  // Chlorine
    Ar  = 17,  // Argon
    K   = 18,  // Potassium
    Ca  = 19,  // Calcium
    Sc  = 20,  // Scandium
    Ti  = 21,  // Titanium
    V   = 22,  // Vanadium
    Cr  = 23,  // Chromium
    Mn  = 24,  // Manganese
    Fe  = 25,  // Iron
    Co  = 26,  // Cobalt
    Ni  = 27,  // Nickel
    Cu  = 28,  // Copper
    Zn  = 29,  // Zinc
    Ga  = 30,  // Gallium
    Ge  = 31,  // Germanium
    As  = 32,  // Arsenic
    Se  = 33,  // Selenium
    Br  = 34,  // Bromine
    Kr  = 35,  // Krypton
    Rb  = 36,  // Rubidium
    Sr  = 37,  // Strontium
    Y   = 38,  // Yttrium
    Zr  = 39,  // Zirconium
    Nb  = 40,  // Niobium
    Mo  = 41,  // Molybdenum
    Tc  = 42,  // Technetium
    Ru  = 43,  // Ruthenium
    Rh  = 44,  // Rhodium
    Pd  = 45,  // Palladium
    Ag  = 46,  // Silver
    Cd  = 47,  // Cadmium
    In  = 48,  // Indium
    Sn  = 49,  // Tin
    Sb  = 50,  // Antimony
    Te  = 51,  // Tellurium
    I   = 52,  // Iodine
    Xe  = 53,  // Xenon
    Cs  = 54,  // Cesium
    Ba  = 55,  // Barium
    La  = 56,  // Lanthanum
    Ce  = 57,  // Cerium
    Pr  = 58,  // Praseodymium
    Nd  = 59,  // Neodymium
    Pm  = 60,  // Promethium
    Sm  = 61,  // Samarium
    Eu  = 62,  // Europium
    Gd  = 63,  // Gadolinium
    Tb  = 64,  // Terbium
    Dy  = 65,  // Dysprosium
    Ho  = 66,  // Holmium
    Er  = 67,  // Erbium
    Tm  = 68,  // Thulium
    Yb  = 69,  // Ytterbium
    Lu  = 70,  // Lutetium
    Hf  = 71,  // Hafnium
    Ta  = 72,  // Tantalum
    W   = 73,  // Tungsten
    Re  = 74,  // Rhenium
    Os  = 75,  // Osmium
    Ir  = 76,  // Iridium
    Pt  = 77,  // Platinum
    Au  = 78,  // Gold
    Hg  = 79,  // Mercury
    Tl  = 80,  // Thallium
    Pb  = 81,  // Lead
    Bi  = 82,  // Bismuth
    Po  = 83,  // Polonium
    At  = 84,  // Astatine
    Rn  = 85,  // Radon
    Fr  = 86,  // Francium
    Ra  = 87,  // Radium
    Ac  = 88,  // Actinium
    Th  = 89,  // Thorium
    Pa  = 90,  // Protactinium
    U   = 91,  // Uranium
    Np  = 92,  // Neptunium
    Pu  = 93,  // Plutonium
    Am  = 94,  // Americium
    Cm  = 95,  // Curium
    Bk  = 96,  // Berkelium
    Cf  = 97,  // Californium
    Es  = 98,  // Einsteinium
    Fm  = 99,  // Fermium
    Md  = 100, // Mendelevium
    No  = 101, // Nobelium
    Lr  = 102, // Lawrencium
    Rf  = 103, // Rutherfordium
    Db  = 104, // Dubnium
    Sg  = 105, // Seaborgium
    Bh  = 106, // Bohrium
    Hs  = 107, // Hassium
    Mt  = 108, // Meitnerium
    Ds  = 109, // Darmstadtium
    Rg  = 110, // Roentgenium
    Cn  = 111, // Copernicium
    Nh  = 112, // Nihonium
    Fl  = 113, // Flerovium
    Mc  = 114, // Moscovium
    Lv  = 115, // Livermorium
    Ts  = 116, // Tennessine
    Og  = 117, // Oganesson
    Uue = 118, // Ununennium
  };

  // this is knowledge that we have about all the elements (that we need at least)
  struct element_description {
    element value;
    std::string_view symbol;
    std::string_view name;
    int number;
    int valence;
  };
  extern const std::array<element_description, 119> ELEMENT_DICTIONARY;

  // utility functions to work with them
  inline const element_description& get_description(const element e) {
    assert(ELEMENT_DICTIONARY[static_cast<int>(e)].value == e);
    return ELEMENT_DICTIONARY[static_cast<int>(e)];
  }
  std::optional<element> parse_element_symbol(const std::string_view symbol);

} // namespace mudock
