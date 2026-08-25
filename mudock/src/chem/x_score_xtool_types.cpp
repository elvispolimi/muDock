#include <algorithm>
#include <iostream>
#include <mudock/chem/sybyl_atom_types.hpp>
#include <mudock/chem/x_score_xtool_types.hpp>
#include <stdexcept>

//===---------------------------------------------------------------------------------------------------------------
// WARNING: This file has been automatically generated from chem/xtool_types.json and subsequently modified manually
//===---------------------------------------------------------------------------------------------------------------

namespace mudock {
  namespace {
    // SYBYL -> X-Tool mapping. Every entry defaults to xtool_ff::Un.
    // Any SYBYL type without an X-Tool equivalent is left untyped.
    constexpr auto make_sybyl_to_xtool_table() {
      std::array<xtool_ff, num_sybyl_atom_types()> table{};
      table.fill(xtool_ff::Un);

      const auto set = [&table](sybyl_atom_type s, xtool_ff x) { table[static_cast<std::size_t>(s)] = x; };

      set(sybyl_atom_type::H, xtool_ff::H);

      set(sybyl_atom_type::C_1, xtool_ff::C1);
      set(sybyl_atom_type::C_2, xtool_ff::C2);
      set(sybyl_atom_type::C_3, xtool_ff::C3);
      set(sybyl_atom_type::C_ar, xtool_ff::Car);
      set(sybyl_atom_type::C_cat, xtool_ff::Ccat);

      set(sybyl_atom_type::N_1, xtool_ff::N1);
      set(sybyl_atom_type::N_2, xtool_ff::N2);
      set(sybyl_atom_type::N_3, xtool_ff::N3);
      set(sybyl_atom_type::N_4, xtool_ff::N4);
      set(sybyl_atom_type::N_ar, xtool_ff::Nar);
      set(sybyl_atom_type::N_am, xtool_ff::Nam);
      set(sybyl_atom_type::N_pl3, xtool_ff::Npl3);

      set(sybyl_atom_type::O_2, xtool_ff::O2);
      set(sybyl_atom_type::O_3, xtool_ff::O3);
      set(sybyl_atom_type::O_co2, xtool_ff::Oco2);

      set(sybyl_atom_type::P_3, xtool_ff::P3);

      set(sybyl_atom_type::S_2, xtool_ff::S2);
      set(sybyl_atom_type::S_3, xtool_ff::S3);
      set(sybyl_atom_type::S_o, xtool_ff::So);
      set(sybyl_atom_type::S_o2, xtool_ff::So2);
      set(sybyl_atom_type::S_O, xtool_ff::So);   // non-standard upper-case descriptor "S.O"
      set(sybyl_atom_type::S_O2, xtool_ff::So2); // non-standard upper-case descriptor "S.O2"

      set(sybyl_atom_type::F, xtool_ff::F);
      set(sybyl_atom_type::Cl, xtool_ff::Cl);
      set(sybyl_atom_type::Br, xtool_ff::Br);
      set(sybyl_atom_type::I, xtool_ff::I);

      set(sybyl_atom_type::Si, xtool_ff::Si);

      // Metal / ion element-symbol fallbacks
      set(sybyl_atom_type::Li, xtool_ff::Li);
      set(sybyl_atom_type::Na, xtool_ff::Na);
      set(sybyl_atom_type::K, xtool_ff::K);
      set(sybyl_atom_type::Mg, xtool_ff::Mg);
      set(sybyl_atom_type::Ca, xtool_ff::Ca);
      set(sybyl_atom_type::Mn, xtool_ff::Mn);
      set(sybyl_atom_type::Fe, xtool_ff::Fe);
      set(sybyl_atom_type::Co, xtool_ff::Co);
      set(sybyl_atom_type::Ni, xtool_ff::Ni);
      set(sybyl_atom_type::Cu, xtool_ff::Cu);
      set(sybyl_atom_type::Zn, xtool_ff::Zn);
      set(sybyl_atom_type::Cd, xtool_ff::Cd);
      set(sybyl_atom_type::Hg, xtool_ff::Hg);
      set(sybyl_atom_type::Al, xtool_ff::Al);
      set(sybyl_atom_type::U, xtool_ff::U);

      return table;
    }

    constexpr auto SYBYL_TO_XTOOL_TABLE = make_sybyl_to_xtool_table();
  } // namespace

  xtool_ff xtool_type_from_sybyl(const sybyl_atom_type type) {
    return SYBYL_TO_XTOOL_TABLE[static_cast<std::size_t>(type)];
  }

  xtool_ff parse_xtool_type(const std::string_view symbol) {
    std::string_view search_symbol = symbol;
    
    if (symbol == "HO" || symbol == "HN") {
        search_symbol = "H";
    }

    const auto element_it = std::find_if(std::begin(XTOOL_FF_DICTIONARY),
                                         std::end(XTOOL_FF_DICTIONARY),
                                         [&search_symbol](const auto& e) { return e.name == search_symbol; });
    if (element_it != std::end(XTOOL_FF_DICTIONARY))
      return element_it->value;
    else {
      std::cout << "Missing xtool type: " << symbol << std::endl;
      throw std::runtime_error(std::string("Missing xtool type: ") + std::string(symbol));
    }
  }
  const std::array<xtool_ff_description, 72> XTOOL_FF_DICTIONARY = {{
    {
      xtool_ff::C3,
      "C3",
      12.01f,
      2.100f,
      0.000f,
      0.000f,
      x_score_hb::H
    },

    {
      xtool_ff::C3x,
      "C3x",
      12.01f,
      2.100f,
      0.000f,
      0.000f,
      x_score_hb::P
    },

    {
      xtool_ff::C3un,
      "C3un",
      12.01f,
      2.100f,
      0.000f,
      0.000f,
      x_score_hb::P
    },

    {
      xtool_ff::C2,
      "C2",
      12.01f,
      1.900f,
      0.000f,
      0.000f,
      x_score_hb::H
    },

    {
      xtool_ff::C2x,
      "C2x",
      12.01f,
      1.900f,
      0.000f,
      0.000f,
      x_score_hb::P
    },

    {
      xtool_ff::C2un,
      "C2un",
      12.01f,
      1.900f,
      0.000f,
      0.000f,
      x_score_hb::P
    },

    {
      xtool_ff::Car,
      "Car",
      12.01f,
      2.000f,
      0.000f,
      0.000f,
      x_score_hb::H
    },

    {
      xtool_ff::Carx,
      "Carx",
      12.01f,
      2.000f,
      0.000f,
      0.000f,
      x_score_hb::P
    },

    {
      xtool_ff::Carun,
      "Carun",
      12.01f,
      2.000f,
      0.000f,
      0.000f,
      x_score_hb::P
    },

    {
      xtool_ff::C1,
      "C1",
      12.01f,
      1.800f,
      0.000f,
      0.000f,
      x_score_hb::P
    },

    {
      xtool_ff::C1x,
      "C1x",
      12.01f,
      1.800f,
      0.000f,
      0.000f,
      x_score_hb::P
    },

    {
      xtool_ff::C1un,
      "C1un",
      12.01f,
      1.800f,
      0.000f,
      0.000f,
      x_score_hb::P
    },

    {
      xtool_ff::Ccat,
      "Ccat",
      12.01f,
      1.900f,
      0.000f,
      1.000f,
      x_score_hb::P
    },

    {
      xtool_ff::N3h,
      "N3h",
      14.01f,
      1.800f,
      0.000f,
      0.000f,
      x_score_hb::D
    },

    {
      xtool_ff::N3,
      "N3",
      14.01f,
      1.800f,
      0.000f,
      0.000f,
      x_score_hb::P
    },

    {
      xtool_ff::N3un,
      "N3un",
      14.01f,
      1.800f,
      0.000f,
      0.000f,
      x_score_hb::P
    },

    {
      xtool_ff::Npl3h,
      "Npl3h",
      14.01f,
      1.750f,
      0.000f,
      0.000f,
      x_score_hb::D
    },

    {
      xtool_ff::Npl3,
      "Npl3",
      14.01f,
      1.750f,
      0.000f,
      0.000f,
      x_score_hb::P
    },

    {
      xtool_ff::Npl3un,
      "Npl3un",
      14.01f,
      1.750f,
      0.000f,
      0.000f,
      x_score_hb::P
    },

    {
      xtool_ff::N2h,
      "N2h",
      14.01f,
      1.750f,
      0.000f,
      0.000f,
      x_score_hb::DA
    },

    {
      xtool_ff::N2,
      "N2",
      14.01f,
      1.750f,
      0.000f,
      0.000f,
      x_score_hb::A
    },

    {
      xtool_ff::N2un,
      "N2un",
      14.01f,
      1.750f,
      0.000f,
      0.000f,
      x_score_hb::P
    },

    {
      xtool_ff::Narh,
      "Narh",
      14.01f,
      1.750f,
      0.000f,
      0.000f,
      x_score_hb::D
    },

    {
      xtool_ff::Nar,
      "Nar",
      14.01f,
      1.750f,
      0.000f,
      0.000f,
      x_score_hb::A
    },

    {
      xtool_ff::Narun,
      "Narun",
      14.01f,
      1.750f,
      0.000f,
      0.000f,
      x_score_hb::P
    },

    {
      xtool_ff::N1,
      "N1",
      14.01f,
      1.750f,
      0.000f,
      0.000f,
      x_score_hb::A
    },

    {
      xtool_ff::N1un,
      "N1un",
      14.01f,
      1.750f,
      0.000f,
      0.000f,
      x_score_hb::P
    },

    {
      xtool_ff::N4,
      "N4",
      14.01f,
      1.800f,
      0.000f,
      1.000f,
      x_score_hb::D
    },

    {
      xtool_ff::O3h,
      "O3h",
      16.00f,
      1.650f,
      0.000f,
      0.000f,
      x_score_hb::DA
    },

    {
      xtool_ff::O3,
      "O3",
      16.00f,
      1.650f,
      0.000f,
      0.000f,
      x_score_hb::A
    },

    {
      xtool_ff::O3un,
      "O3un",
      16.00f,
      1.650f,
      0.000f,
      0.000f,
      x_score_hb::P
    },

    {
      xtool_ff::O2,
      "O2",
      16.00f,
      1.550f,
      0.000f,
      0.000f,
      x_score_hb::A
    },

    {
      xtool_ff::O2un,
      "O2un",
      16.00f,
      1.550f,
      0.000f,
      0.000f,
      x_score_hb::P
    },

    {
      xtool_ff::Oco2,
      "Oco2",
      16.00f,
      1.550f,
      0.000f,
      -0.500f,
      x_score_hb::DA
    },

    {
      xtool_ff::S3h,
      "S3h",
      32.07f,
      2.100f,
      0.000f,
      0.000f,
      x_score_hb::H
    },

    {
      xtool_ff::S3,
      "S3",
      32.07f,
      2.100f,
      0.000f,
      0.000f,
      x_score_hb::H
    },

    {
      xtool_ff::S3un,
      "S3un",
      32.07f,
      2.100f,
      0.000f,
      0.000f,
      x_score_hb::P
    },

    {
      xtool_ff::S2,
      "S2",
      32.07f,
      2.000f,
      0.000f,
      0.000f,
      x_score_hb::P
    },

    {
      xtool_ff::S2un,
      "S2un",
      32.07f,
      2.000f,
      0.000f,
      0.000f,
      x_score_hb::P
    },

    {
      xtool_ff::So,
      "So",
      32.07f,
      2.000f,
      0.000f,
      0.000f,
      x_score_hb::P
    },

    {
      xtool_ff::P3,
      "P3",
      30.97f,
      2.000f,
      0.000f,
      0.000f,
      x_score_hb::P
    },

    {
      xtool_ff::F,
      "F",
      19.00f,
      1.500f,
      0.000f,
      0.000f,
      x_score_hb::P
    },

    {
      xtool_ff::Cl,
      "Cl",
      35.45f,
      1.750f,
      0.000f,
      0.000f,
      x_score_hb::H
    },

    {
      xtool_ff::Br,
      "Br",
      79.90f,
      1.900f,
      0.000f,
      0.000f,
      x_score_hb::H
    },

    {
      xtool_ff::I,
      "I",
      126.90f,
      2.050f,
      0.000f,
      0.000f,
      x_score_hb::H
    },

    {
      xtool_ff::H,
      "H",
      1.00f,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N
    },

    {
      xtool_ff::Hhb,
      "Hhb",
      1.00f,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH
    },

    {
      xtool_ff::Si,
      "Si",
      28.09f,
      2.000f,
      0.000f,
      0.000f,
      x_score_hb::N
    },

    {
      xtool_ff::Ow,
      "Ow",
      16.00f,
      1.750f,
      0.000f,
      0.000f,
      x_score_hb::DA
    },

    {
      xtool_ff::Mplus,
      "Mplus",
      0.00f,
      1.250f,
      0.000f,
      2.000f,
      x_score_hb::M
    },

    {
      xtool_ff::Un,
      "Un",
      0.00f,
      0.000f,
      0.000f,
      0.000f,
      x_score_hb::N
    },

    // non xtool types added for type mismatch in residue
    //xlogp
    {
      xtool_ff::Nam,
      "Nam",
      14.01f,
      1.750f,
      0.000f,
      0.000f,
      x_score_hb::D
    },

    //xlogp
    {
      xtool_ff::So2,
      "So2",
      32.07f,
      2.000f,
      0.000f,
      0.000f,
      x_score_hb::P
    },

    //ions
    {
      xtool_ff::Li,
      "Li",
      0.00f,
      1.250f,
      0.000f,
      1.000f,
      x_score_hb::M
    },

    {
      xtool_ff::Na,
      "Na",
      0.00f,
      1.250f,
      0.000f,
      1.000f,
      x_score_hb::M
    },

    {
      xtool_ff::K,
      "K",
      0.00f,
      1.250f,
      0.000f,
      1.000f,
      x_score_hb::M
    },

    {
      xtool_ff::Mg,
      "Mg",
      0.00f,
      1.250f,
      0.000f,
      2.000f,
      x_score_hb::M
    },

    {
      xtool_ff::Ca,
      "Ca",
      0.00f,
      1.250f,
      0.000f,
      2.000f,
      x_score_hb::M
    },

    {
      xtool_ff::Mn,
      "Mn",
      0.00f,
      1.250f,
      0.000f,
      2.000f,
      x_score_hb::M
    },

    {
      xtool_ff::Fe,
      "Fe",
      0.00f,
      1.250f,
      0.000f,
      2.000f,
      x_score_hb::M
    },

    {
      xtool_ff::Co,
      "Co",
      0.00f,
      1.250f,
      0.000f,
      2.000f,
      x_score_hb::M
    },

    {
      xtool_ff::Ni,
      "Ni",
      0.00f,
      1.250f,
      0.000f,
      2.000f,
      x_score_hb::M
    },

    {
      xtool_ff::Cu,
      "Cu",
      0.00f,
      1.250f,
      0.000f,
      2.000f,
      x_score_hb::M
    },

    {
      xtool_ff::Zn,
      "Zn",
      0.00f,
      1.250f,
      0.000f,
      2.000f,
      x_score_hb::M
    },

    {
      xtool_ff::Cd,
      "Cd",
      0.00f,
      1.250f,
      0.000f,
      2.000f,
      x_score_hb::M
    },

    {
      xtool_ff::Hg,
      "Hg",
      0.00f,
      1.250f,
      0.000f,
      2.000f,
      x_score_hb::M
    },

    {
      xtool_ff::Al,
      "Al",
      0.00f,
      1.250f,
      0.000f,
      3.000f,
      x_score_hb::M
    },

    {
      xtool_ff::U,
      "U",
      0.00f,
      1.250f,
      0.000f,
      3.000f,
      x_score_hb::M
    },

    {
      xtool_ff::Fminus,
      "Fminus",
      19.00f,
      1.500f,
      0.000f,
      -1.000f,
      x_score_hb::P
    },

    {
      xtool_ff::Clminus,
      "Clminus",
      35.45f,
      1.750f,
      0.000f,
      -1.000f,
      x_score_hb::P
    },

    {
      xtool_ff::Brminus,
      "Brminus",
      79.90f,
      1.900f,
      0.000f,
      -1.000f,
      x_score_hb::P
    },

    {
      xtool_ff::Iminus,
      "Iminus",
      126.90f,
      2.050f,
      0.000f,
      -1.000f,
      x_score_hb::P
    }

  }};

} // namespace mudock