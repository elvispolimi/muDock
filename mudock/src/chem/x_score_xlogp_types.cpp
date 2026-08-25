#include <algorithm>
#include <mudock/chem/x_score_xlogp_types.hpp>
#include <stdexcept>

//===------------------------------------------------------------------------------------------------------
// WARNING: This file has been automatically generated from chem/xlogp_types.json and subsequently modified manually
//===------------------------------------------------------------------------------------------------------

namespace mudock {
  xlogp_ff parse_xlogp_type(const std::string_view symbol) {
    const auto element_it = std::find_if(std::begin(XLOGP_FF_DICTIONARY),
                                         std::end(XLOGP_FF_DICTIONARY),
                                         [&symbol](const auto& e) { return e.name == symbol; });
    if (element_it != std::end(XLOGP_FF_DICTIONARY))
      return element_it->value;
    else
      throw std::runtime_error("Missing xlogp type");
  }
  const std::array<xlogp_ff_description, 112> XLOGP_FF_DICTIONARY = {{

    {
      xlogp_ff::C_3_h3_pi_eq_0,
      "C_3_h3_pi_eq_0",
      x_score_hb::H,
      0.528f
    },

    {
      xlogp_ff::C_3_h3_pi_eq_1,
      "C_3_h3_pi_eq_1",
      x_score_hb::H,
      0.267f
    },

    {
      xlogp_ff::C_3_h3_x,
      "C_3_h3_x",
      x_score_hb::P,
      -0.032f
    },

    {
      xlogp_ff::C_3_h2_pi_eq_0,
      "C_3_h2_pi_eq_0",
      x_score_hb::H,
      0.358f
    },

    {
      xlogp_ff::C_3_h2_pi_eq_1,
      "C_3_h2_pi_eq_1",
      x_score_hb::H,
      -0.008f
    },

    {
      xlogp_ff::C_3_h2_pi_eq_2,
      "C_3_h2_pi_eq_2",
      x_score_hb::P,
      -0.185f
    },

    {
      xlogp_ff::C_3_h2_x_pi_eq_0,
      "C_3_h2_x_pi_eq_0",
      x_score_hb::P,
      -0.137f
    },

    {
      xlogp_ff::C_3_h2_x_pi_eq_1,
      "C_3_h2_x_pi_eq_1",
      x_score_hb::P,
      -0.303f
    },

    {
      xlogp_ff::C_3_h2_x_pi_eq_2,
      "C_3_h2_x_pi_eq_2",
      x_score_hb::P,
      -0.815f
    },

    {
      xlogp_ff::C_3_h_pi_eq_0,
      "C_3_h_pi_eq_0",
      x_score_hb::H,
      0.127f
    },

    {
      xlogp_ff::C_3_h_pi_eq_1,
      "C_3_h_pi_eq_1",
      x_score_hb::H,
      -0.243f
    },

    {
      xlogp_ff::C_3_h_pi_gt_1,
      "C_3_h_pi_gt_1",
      x_score_hb::P,
      -0.499f
    },

    {
      xlogp_ff::C_3_h_x_pi_eq_0,
      "C_3_h_x_pi_eq_0",
      x_score_hb::P,
      -0.205f
    },

    {
      xlogp_ff::C_3_h_x_pi_eq_1,
      "C_3_h_x_pi_eq_1",
      x_score_hb::P,
      -0.305f
    },

    {
      xlogp_ff::C_3_h_x_pi_gt_1,
      "C_3_h_x_pi_gt_1",
      x_score_hb::P,
      -0.709f
    },

    {
      xlogp_ff::C_3_pi_eq_0,
      "C_3_pi_eq_0",
      x_score_hb::H,
      -0.006f
    },

    {
      xlogp_ff::C_3_pi_eq_1,
      "C_3_pi_eq_1",
      x_score_hb::H,
      -0.570f
    },

    {
      xlogp_ff::C_3_pi_gt_1,
      "C_3_pi_gt_1",
      x_score_hb::P,
      -0.317f
    },

    {
      xlogp_ff::C_3_x_pi_eq_0,
      "C_3_x_pi_eq_0",
      x_score_hb::P,
      -0.316f
    },

    {
      xlogp_ff::C_3_x_pi_gt_0,
      "C_3_x_pi_gt_0",
      x_score_hb::P,
      -0.723f
    },

    {
      xlogp_ff::C_3_unknown,
      "C_3_unknown",
      x_score_hb::P,
      0.528f
    },

    {
      xlogp_ff::C_2_h2,
      "C_2_h2",
      x_score_hb::H,
      0.420f
    },

    {
      xlogp_ff::C_2_h_pi_eq_0,
      "C_2_h_pi_eq_0",
      x_score_hb::H,
      0.466f
    },

    {
      xlogp_ff::C_2_h_pi_eq_1,
      "C_2_h_pi_eq_1",
      x_score_hb::H,
      0.136f
    },

    {
      xlogp_ff::C_2_h_x_pi_eq_0,
      "C_2_h_x_pi_eq_0",
      x_score_hb::P,
      0.001f
    },

    {
      xlogp_ff::C_2_h_x_pi_eq_1,
      "C_2_h_x_pi_eq_1",
      x_score_hb::P,
      -0.310f
    },

    {
      xlogp_ff::C_2_pi_eq_0,
      "C_2_pi_eq_0",
      x_score_hb::H,
      0.050f
    },

    {
      xlogp_ff::C_2_pi_gt_0,
      "C_2_pi_gt_0",
      x_score_hb::H,
      0.013f
    },

    {
      xlogp_ff::C_2_x_pi_eq_0,
      "C_2_x_pi_eq_0",
      x_score_hb::P,
      -0.030f
    },

    {
      xlogp_ff::C_2_x_pi_gt_0,
      "C_2_x_pi_gt_0",
      x_score_hb::P,
      -0.027f
    },

    {
      xlogp_ff::C_2_x2_pi_eq_0,
      "C_2_x2_pi_eq_0",
      x_score_hb::P,
      0.005f
    },

    {
      xlogp_ff::C_2_x2_pi_gt_0,
      "C_2_x2_pi_gt_0",
      x_score_hb::P,
      -0.315f
    },

    {
      xlogp_ff::C_2_unknown,
      "C_2_unknown",
      x_score_hb::P,
      0.050f
    },

    {
      xlogp_ff::C_ar_h,
      "C_ar_h",
      x_score_hb::H,
      0.337f
    },

    {
      xlogp_ff::C_ar_h_X,
      "C_ar_h_X",
      x_score_hb::H,
      0.126f
    },

    {
      xlogp_ff::C_ar,
      "C_ar",
      x_score_hb::H,
      0.296f
    },

    {
      xlogp_ff::C_ar_x,
      "C_ar_x",
      x_score_hb::P,
      -0.151f
    },

    {
      xlogp_ff::C_ar_X,
      "C_ar_X",
      x_score_hb::H,
      0.174f
    },

    {
      xlogp_ff::C_ar_X_x,
      "C_ar_X_x",
      x_score_hb::H,
      0.366f
    },

    {
      xlogp_ff::C_ar_unknown,
      "C_ar_unknown",
      x_score_hb::P,
      0.296f
    },

    {
      xlogp_ff::C_1_h,
      "C_1_h",
      x_score_hb::H,
      0.209f
    },

    {
      xlogp_ff::C_1,
      "C_1",
      x_score_hb::H,
      0.330f
    },

    {
      xlogp_ff::C_1_eq_eq,
      "C_1_eq_eq",
      x_score_hb::H,
      2.073f
    },

    {
      xlogp_ff::C_1_unknown,
      "C_1_unknown",
      x_score_hb::H,
      0.330f
    },

    {
      xlogp_ff::C_cat,
      "C_cat",
      x_score_hb::P,
      -0.315f
    },

    {
      xlogp_ff::N_3_h2_pi_eq_0,
      "N_3_h2_pi_eq_0",
      x_score_hb::D,
      -0.534f
    },

    {
      xlogp_ff::N_3_h2_pi_eq_1,
      "N_3_h2_pi_eq_1",
      x_score_hb::D,
      -0.329f
    },

    {
      xlogp_ff::N_3_h2_x,
      "N_3_h2_x",
      x_score_hb::D,
      -1.082f
    },

    {
      xlogp_ff::N_3_h_pi_eq_0,
      "N_3_h_pi_eq_0",
      x_score_hb::D,
      -0.112f
    },

    {
      xlogp_ff::N_3_h_pi_gt_0,
      "N_3_h_pi_gt_0",
      x_score_hb::D,
      0.166f
    },

    {
      xlogp_ff::N_3_h_ring,
      "N_3_h_ring",
      x_score_hb::D,
      0.545f
    },

    {
      xlogp_ff::N_3_h_x,
      "N_3_h_x",
      x_score_hb::D,
      0.324f
    },

    {
      xlogp_ff::N_3_h_x_ring,
      "N_3_h_x_ring",
      x_score_hb::D,
      0.153f
    },

    {
      xlogp_ff::N_3_pi_eq_0,
      "N_3_pi_eq_0",
      x_score_hb::P,
      0.159f
    },

    {
      xlogp_ff::N_3_pi_gt_0,
      "N_3_pi_gt_0",
      x_score_hb::P,
      0.761f
    },

    {
      xlogp_ff::N_3_ring,
      "N_3_ring",
      x_score_hb::P,
      0.881f
    },

    {
      xlogp_ff::N_3_x,
      "N_3_x",
      x_score_hb::P,
      -0.239f
    },

    {
      xlogp_ff::N_3_x_ring,
      "N_3_x_ring",
      x_score_hb::P,
      -0.010f
    },

    {
      xlogp_ff::N_3_unknown,
      "N_3_unknown",
      x_score_hb::P,
      0.159f
    },

    {
      xlogp_ff::N_am_h2,
      "N_am_h2",
      x_score_hb::D,
      -0.646f
    },

    {
      xlogp_ff::N_am_h,
      "N_am_h",
      x_score_hb::D,
      -0.096f
    },

    {
      xlogp_ff::N_am_h_x,
      "N_am_h_x",
      x_score_hb::D,
      -0.044f
    },

    {
      xlogp_ff::N_am,
      "N_am",
      x_score_hb::P,
      0.078f
    },

    {
      xlogp_ff::N_am_x,
      "N_am_x",
      x_score_hb::P,
      -0.118f
    },

    {
      xlogp_ff::N_am_unknown,
      "N_am_unknown",
      x_score_hb::P,
      0.078f
    },

    {
      xlogp_ff::N_2_eq_C_pi_eq_0,
      "N_2_eq_C_pi_eq_0",
      x_score_hb::A,
      0.007f
    },

    {
      xlogp_ff::N_2_eq_C_pi_eq_1,
      "N_2_eq_C_pi_eq_1",
      x_score_hb::A,
      -0.275f
    },

    {
      xlogp_ff::N_2_eq_C_x_pi_eq_0,
      "N_2_eq_C_x_pi_eq_0",
      x_score_hb::A,
      0.366f
    },

    {
      xlogp_ff::N_2_eq_C_x_pi_eq_1,
      "N_2_eq_C_x_pi_eq_1",
      x_score_hb::A,
      0.251f
    },

    {
      xlogp_ff::N_2_eq_N,
      "N_2_eq_N",
      x_score_hb::A,
      0.536f
    },

    {
      xlogp_ff::N_2_eq_N_x,
      "N_2_eq_N_x",
      x_score_hb::A,
      -0.597f
    },

    {
      xlogp_ff::N_2_o,
      "N_2_o",
      x_score_hb::P,
      0.427f
    },

    {
      xlogp_ff::N_2_o2,
      "N_2_o2",
      x_score_hb::P,
      1.178f
    },

    {
      xlogp_ff::N_2_unknown,
      "N_2_unknown",
      x_score_hb::P,
      0.007f
    },

    {
      xlogp_ff::N_ar,
      "N_ar",
      x_score_hb::A,
      -0.493f
    },

    {
      xlogp_ff::N_1,
      "N_1",
      x_score_hb::A,
      -0.566f
    },

    {
      xlogp_ff::N_4,
      "N_4",
      x_score_hb::D,
      -0.534f
    },

    {
      xlogp_ff::O_3_h_pi_eq_0,
      "O_3_h_pi_eq_0",
      x_score_hb::DA,
      -0.467f
    },

    {
      xlogp_ff::O_3_h_pi_eq_1,
      "O_3_h_pi_eq_1",
      x_score_hb::DA,
      0.082f
    },

    {
      xlogp_ff::O_3_h_x,
      "O_3_h_x",
      x_score_hb::DA,
      -0.522f
    },

    {
      xlogp_ff::O_3_pi_eq_0,
      "O_3_pi_eq_0",
      x_score_hb::A,
      0.084f
    },

    {
      xlogp_ff::O_3_pi_gt_0,
      "O_3_pi_gt_0",
      x_score_hb::A,
      0.435f
    },

    {
      xlogp_ff::O_3_x,
      "O_3_x",
      x_score_hb::A,
      0.105f
    },

    {
      xlogp_ff::O_3_unknown,
      "O_3_unknown",
      x_score_hb::P,
      0.084f
    },

    {
      xlogp_ff::O_2,
      "O_2",
      x_score_hb::A,
      -0.399f
    },

    {
      xlogp_ff::O_co2,
      "O_co2",
      x_score_hb::DA,
      -0.399f
    },

    {
      xlogp_ff::S_3_h,
      "S_3_h",
      x_score_hb::H,
      0.419f
    },

    {
      xlogp_ff::S_3,
      "S_3",
      x_score_hb::H,
      0.255f
    },

    {
      xlogp_ff::S_3_unknown,
      "S_3_unknown",
      x_score_hb::H,
      0.255f
    },

    {
      xlogp_ff::S_2,
      "S_2",
      x_score_hb::P,
      -0.148f
    },

    {
      xlogp_ff::S_o,
      "S_o",
      x_score_hb::P,
      -1.375f
    },

    {
      xlogp_ff::S_o2,
      "S_o2",
      x_score_hb::P,
      -0.168f
    },

    {
      xlogp_ff::P_3_eq_O,
      "P_3_eq_O",
      x_score_hb::P,
      -0.447f
    },

    {
      xlogp_ff::P_3_eq_S,
      "P_3_eq_S",
      x_score_hb::P,
      1.253f
    },

    {
      xlogp_ff::P_3_unknown,
      "P_3_unknown",
      x_score_hb::P,
      -0.447f
    },

    {
      xlogp_ff::F_pi_eq_0,
      "F_pi_eq_0",
      x_score_hb::H,
      0.375f
    },

    {
      xlogp_ff::F_pi_eq_1,
      "F_pi_eq_1",
      x_score_hb::H,
      0.202f
    },

    {
      xlogp_ff::F_unknown,
      "F_unknown",
      x_score_hb::H,
      0.375f
    },

    {
      xlogp_ff::Cl_pi_eq_0,
      "Cl_pi_eq_0",
      x_score_hb::H,
      0.512f
    },

    {
      xlogp_ff::Cl_pi_eq_1,
      "Cl_pi_eq_1",
      x_score_hb::H,
      0.663f
    },

    {
      xlogp_ff::Cl_unknown,
      "Cl_unknown",
      x_score_hb::H,
      0.512f
    },

    {
      xlogp_ff::Br_pi_eq_0,
      "Br_pi_eq_0",
      x_score_hb::H,
      0.850f
    },

    {
      xlogp_ff::Br_pi_eq_1,
      "Br_pi_eq_1",
      x_score_hb::H,
      0.839f
    },

    {
      xlogp_ff::Br_unknown,
      "Br_unknown",
      x_score_hb::H,
      0.850f
    },

    {
      xlogp_ff::I_pi_eq_0,
      "I_pi_eq_0",
      x_score_hb::H,
      1.050f
    },

    {
      xlogp_ff::I_pi_eq_1,
      "I_pi_eq_1",
      x_score_hb::H,
      1.109f
    },

    {
      xlogp_ff::I_unknown,
      "I_unknown",
      x_score_hb::H,
      1.050f
    },

    {
      xlogp_ff::H,
      "H",
      x_score_hb::N,
      0.000f
    },

    {
      xlogp_ff::H_hb,
      "H_hb",
      x_score_hb::DH,
      0.000f
    },

    {
      xlogp_ff::Si,
      "Si",
      x_score_hb::N,
      0.000f
    },

    {
      xlogp_ff::Un,
      "Un",
      x_score_hb::N,
      0.000f
    },

    {
      xlogp_ff::Du,
      "Du",
      x_score_hb::N,
      0.000f
    }

  }};

} // namespace mudock