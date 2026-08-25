#pragma once

#include <array>
#include <cassert>
#include <mudock/chem/x_score_hb.hpp>
#include <mudock/type_alias.hpp>
#include <string_view>

//===------------------------------------------------------------------------------------------------------
// WARNING: This file has been automatically generated from chem/x_score_xlogp_types.json and subsequently modified manually
//===------------------------------------------------------------------------------------------------------

namespace mudock {

  enum class xlogp_ff : int {

    C_3_h3_pi_eq_0 = 0,
  
    C_3_h3_pi_eq_1 = 1,

    C_3_h3_x = 2,

    C_3_h2_pi_eq_0 = 3,

    C_3_h2_pi_eq_1 = 4,

    C_3_h2_pi_eq_2 = 5,

    C_3_h2_x_pi_eq_0 = 6,

    C_3_h2_x_pi_eq_1 = 7,

    C_3_h2_x_pi_eq_2 = 8,

    C_3_h_pi_eq_0 = 9,

    C_3_h_pi_eq_1 = 10,

    C_3_h_pi_gt_1 = 11,

    C_3_h_x_pi_eq_0 = 12,

    C_3_h_x_pi_eq_1 = 13,

    C_3_h_x_pi_gt_1 = 14,

    C_3_pi_eq_0 = 15,

    C_3_pi_eq_1 = 16,

    C_3_pi_gt_1 = 17,

    C_3_x_pi_eq_0 = 18,

    C_3_x_pi_gt_0 = 19,

    C_3_unknown = 20,

    C_2_h2 = 21,

    C_2_h_pi_eq_0 = 22,

    C_2_h_pi_eq_1 = 23,

    C_2_h_x_pi_eq_0 = 24,

    C_2_h_x_pi_eq_1 = 25,

    C_2_pi_eq_0 = 26,

    C_2_pi_gt_0 = 27,

    C_2_x_pi_eq_0 = 28,

    C_2_x_pi_gt_0 = 29,

    C_2_x2_pi_eq_0 = 30,

    C_2_x2_pi_gt_0 = 31,

    C_2_unknown = 32,

    C_ar_h = 33,

    C_ar_h_X = 34,

    C_ar = 35,

    C_ar_x = 36,

    C_ar_X = 37,

    C_ar_X_x = 38,

    C_ar_unknown = 39,

    C_1_h = 40,

    C_1 = 41,

    C_1_eq_eq = 42,

    C_1_unknown = 43,

    C_cat = 44,

    N_3_h2_pi_eq_0 = 45,

    N_3_h2_pi_eq_1 = 46,

    N_3_h2_x = 47,

    N_3_h_pi_eq_0 = 48,

    N_3_h_pi_gt_0 = 49,

    N_3_h_ring = 50,

    N_3_h_x = 51,

    N_3_h_x_ring = 52,

    N_3_pi_eq_0 = 53,

    N_3_pi_gt_0 = 54,

    N_3_ring = 55,

    N_3_x = 56,

    N_3_x_ring = 57,

    N_3_unknown = 58,

    N_am_h2 = 59,

    N_am_h = 60,

    N_am_h_x = 61,

    N_am = 62,

    N_am_x = 63,

    N_am_unknown = 64,

    N_2_eq_C_pi_eq_0 = 65,

    N_2_eq_C_pi_eq_1 = 66,

    N_2_eq_C_x_pi_eq_0 = 67,

    N_2_eq_C_x_pi_eq_1 = 68,

    N_2_eq_N = 69,

    N_2_eq_N_x = 70,

    N_2_o = 71,

    N_2_o2 = 72,

    N_2_unknown = 73,

    N_ar = 74,

    N_1 = 75,

    N_4 = 76,

    O_3_h_pi_eq_0 = 77,

    O_3_h_pi_eq_1 = 78,

    O_3_h_x = 79,

    O_3_pi_eq_0 = 80,

    O_3_pi_gt_0 = 81,

    O_3_x = 82,

    O_3_unknown = 83,

    O_2 = 84,

    O_co2 = 85,

    S_3_h = 86,

    S_3 = 87,

    S_3_unknown = 88,

    S_2 = 89,

    S_o = 90,

    S_o2 = 91,

    P_3_eq_O = 92,

    P_3_eq_S = 93,

    P_3_unknown = 94,

    F_pi_eq_0 = 95,

    F_pi_eq_1 = 96,

    F_unknown = 97,

    Cl_pi_eq_0 = 98,

    Cl_pi_eq_1 = 99,

    Cl_unknown = 100,

    Br_pi_eq_0 = 101,

    Br_pi_eq_1 = 102,

    Br_unknown = 103,

    I_pi_eq_0 = 104,

    I_pi_eq_1 = 105,

    I_unknown = 106,

    H = 107,

    H_hb = 108,

    Si = 109,

    Un = 110,

    Du = 111,

  };

  struct xlogp_ff_description {
    xlogp_ff value;
    std::string_view name;
    x_score_hb hbond;
    fp_type hydrophobic_scale;
  };
  static constexpr auto num_xlogp_ff() { return 112; }
  extern const std::array<xlogp_ff_description, num_xlogp_ff()> XLOGP_FF_DICTIONARY;

  inline const xlogp_ff_description& get_description(const xlogp_ff a) {
    assert(XLOGP_FF_DICTIONARY[static_cast<int>(a)].value == a);
    return XLOGP_FF_DICTIONARY[static_cast<int>(a)];
  }

  xlogp_ff parse_xlogp_type(const std::string_view symbol);


} // namespace mudock