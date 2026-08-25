#include <mudock/chem/x_score_residue_xtool_types.hpp>
#include <mudock/chem/residue.hpp>


namespace mudock {

  static const xtool_residue_atom_description ACE_atoms[] = {

    {
      "C",
      xtool_ff::C2,
      xtool_ff::C2x,
      1.900f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.030f,
      0.023f,
      0,
      "CP"
    },

    {
      "O",
      xtool_ff::O2,
      xtool_ff::O2,
      1.550f,
      0.000f,
      0.000f,
      x_score_hb::A,
      -0.399f,
      -0.035f,
      0,
      "OA"
    },

    {
      "CA",
      xtool_ff::C3,
      xtool_ff::C3,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.267f,
      0.021f,
      0,
      "CF"
    },

    {
      "HA1",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HA2",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HA3",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    }

  };

  static const xtool_residue_bond_description ACE_bonds[] = {

    { "C", "O", 2 },

    { "C", "CA", 1 },

    { "CA", "HA1", 1 },

    { "CA", "HA2", 1 },

    { "CA", "HA3", 1 }

  };

  static const xtool_residue_atom_description ALA_atoms[] = {

    {
      "N",
      xtool_ff::Nam,
      xtool_ff::Npl3h,
      1.750f,
      0.000f,
      0.000f,
      x_score_hb::D,
      -0.096f,
      -0.052f,
      0,
      "ND"
    },

    {
      "H",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      -0.115f,
      0,
      "HH"
    },

    {
      "C",
      xtool_ff::C2,
      xtool_ff::C2x,
      1.900f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.030f,
      0.023f,
      0,
      "CP"
    },

    {
      "O",
      xtool_ff::O2,
      xtool_ff::O2,
      1.550f,
      0.000f,
      0.000f,
      x_score_hb::A,
      -0.399f,
      -0.035f,
      0,
      "OA"
    },

    {
      "CA",
      xtool_ff::C3,
      xtool_ff::C3x,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.305f,
      0.005f,
      0,
      "CP"
    },

    {
      "HA",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CB",
      xtool_ff::C3,
      xtool_ff::C3,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.528f,
      0.021f,
      0,
      "CF"
    },

    {
      "HB1",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HB2",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HB3",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    }

  };

  static const xtool_residue_bond_description ALA_bonds[] = {

    { "N", "H", 1 },

    { "N", "CA", 1 },

    { "C", "O", 2 },

    { "C", "CA", 1 },

    { "CA", "HA", 1 },

    { "CA", "CB", 1 },

    { "CB", "HB1", 1 },

    { "CB", "HB2", 1 },

    { "CB", "HB3", 1 }

  };

  static const xtool_residue_atom_description ARG_atoms[] = {

    {
      "N",
      xtool_ff::Nam,
      xtool_ff::Npl3h,
      1.750f,
      0.000f,
      0.000f,
      x_score_hb::D,
      -0.096f,
      -0.052f,
      0,
      "ND"
    },

    {
      "H",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      -0.115f,
      0,
      "HH"
    },

    {
      "C",
      xtool_ff::C2,
      xtool_ff::C2x,
      1.900f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.030f,
      0.023f,
      0,
      "CP"
    },

    {
      "O",
      xtool_ff::O2,
      xtool_ff::O2,
      1.550f,
      0.000f,
      0.000f,
      x_score_hb::A,
      -0.399f,
      -0.035f,
      0,
      "OA"
    },

    {
      "CA",
      xtool_ff::C3,
      xtool_ff::C3x,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.305f,
      0.005f,
      0,
      "CP"
    },

    {
      "HA",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CB",
      xtool_ff::C3,
      xtool_ff::C3,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.358f,
      0.021f,
      0,
      "CF"
    },

    {
      "HB1",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HB2",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CG",
      xtool_ff::C3,
      xtool_ff::C3,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.358f,
      0.021f,
      0,
      "CF"
    },

    {
      "HG1",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HG2",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CD",
      xtool_ff::C3,
      xtool_ff::C3x,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.137f,
      0.005f,
      0,
      "CN"
    },

    {
      "HD1",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HD2",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "NE",
      xtool_ff::Npl3,
      xtool_ff::Npl3h,
      1.750f,
      0.000f,
      0.000f,
      x_score_hb::D,
      -0.096f,
      0.041f,
      0,
      "NC"
    },

    {
      "HE",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      -0.115f,
      0,
      "HH"
    },

    {
      "CZ",
      xtool_ff::Ccat,
      xtool_ff::Ccat,
      1.800f,
      0.000f,
      0.000f,
      x_score_hb::P,
      0.005f,
      0.023f,
      0,
      "CN"
    },

    {
      "NH1",
      xtool_ff::Npl3,
      xtool_ff::Npl3h,
      1.750f,
      0.000f,
      0.500f,
      x_score_hb::D,
      -0.646f,
      0.041f,
      0,
      "NC"
    },

    {
      "HH11",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HH12",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "NH2",
      xtool_ff::Npl3,
      xtool_ff::Npl3h,
      1.750f,
      0.000f,
      0.500f,
      x_score_hb::D,
      -0.646f,
      0.041f,
      0,
      "NC"
    },

    {
      "HH21",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HH22",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      0.000f,
      0,
      "HH"
    }

  };

  static const xtool_residue_bond_description ARG_bonds[] = {

    { "N", "H", 1 },

    { "N", "CA", 1 },

    { "C", "O", 2 },

    { "C", "CA", 1 },

    { "CA", "HA", 1 },

    { "CA", "CB", 1 },

    { "CB", "HB1", 1 },

    { "CB", "HB2", 1 },

    { "CB", "CG", 1 },

    { "CG", "HG1", 1 },

    { "CG", "HG2", 1 },

    { "CG", "CD", 1 },

    { "CD", "HD1", 1 },

    { "CD", "HD2", 1 },

    { "CD", "NE", 1 },

    { "NE", "HE", 1 },

    { "NE", "CZ", 1 },

    { "CZ", "NH1", 1 },

    { "CZ", "NH2", 1 },

    { "NH1", "HH11", 1 },

    { "NH1", "HH12", 1 },

    { "NH2", "HH21", 1 },

    { "NH2", "HH22", 1 }

  };

  static const xtool_residue_atom_description ASN_atoms[] = {

    {
      "N",
      xtool_ff::Nam,
      xtool_ff::Npl3h,
      1.750f,
      0.000f,
      0.000f,
      x_score_hb::D,
      -0.096f,
      -0.052f,
      0,
      "ND"
    },

    {
      "H",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      -0.115f,
      0,
      "HH"
    },

    {
      "C",
      xtool_ff::C2,
      xtool_ff::C2x,
      1.900f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.030f,
      0.023f,
      0,
      "CP"
    },

    {
      "O",
      xtool_ff::O2,
      xtool_ff::O2,
      1.550f,
      0.000f,
      0.000f,
      x_score_hb::A,
      -0.399f,
      -0.035f,
      0,
      "OA"
    },

    {
      "CA",
      xtool_ff::C3,
      xtool_ff::C3x,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.305f,
      0.005f,
      0,
      "CP"
    },

    {
      "HA",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CB",
      xtool_ff::C3,
      xtool_ff::C3,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::H,
      -0.008f,
      0.021f,
      0,
      "CF"
    },

    {
      "HB1",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HB2",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CG",
      xtool_ff::C2,
      xtool_ff::C2x,
      1.900f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.030f,
      0.023f,
      0,
      "CP"
    },

    {
      "OD1",
      xtool_ff::O2,
      xtool_ff::O2,
      1.550f,
      0.000f,
      0.000f,
      x_score_hb::A,
      -0.399f,
      -0.035f,
      0,
      "OA"
    },

    {
      "ND2",
      xtool_ff::Nam,
      xtool_ff::Npl3h,
      1.750f,
      0.000f,
      0.000f,
      x_score_hb::D,
      -0.646f,
      -0.052f,
      0,
      "ND"
    },

    {
      "HD21",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      -0.115f,
      0,
      "HH"
    },

    {
      "HD22",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      -0.115f,
      0,
      "HH"
    }

  };

  static const xtool_residue_bond_description ASN_bonds[] = {

    { "N", "H", 1 },

    { "N", "CA", 1 },

    { "C", "O", 2 },

    { "C", "CA", 1 },

    { "CA", "HA", 1 },

    { "CA", "CB", 1 },

    { "CB", "HB1", 1 },

    { "CB", "HB2", 1 },

    { "CB", "CG", 1 },

    { "CG", "OD1", 2 },

    { "CG", "ND2", 3 },

    { "ND2", "HD21", 1 },

    { "ND2", "HD22", 1 }

  };

  static const xtool_residue_atom_description ASP_atoms[] = {

    {
      "N",
      xtool_ff::Nam,
      xtool_ff::Npl3h,
      1.750f,
      0.000f,
      0.000f,
      x_score_hb::D,
      -0.096f,
      -0.052f,
      0,
      "ND"
    },

    {
      "H",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      -0.115f,
      0,
      "HH"
    },

    {
      "C",
      xtool_ff::C2,
      xtool_ff::C2x,
      1.900f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.030f,
      0.023f,
      0,
      "CP"
    },

    {
      "O",
      xtool_ff::O2,
      xtool_ff::O2,
      1.550f,
      0.000f,
      0.000f,
      x_score_hb::A,
      -0.399f,
      -0.035f,
      0,
      "OA"
    },

    {
      "CA",
      xtool_ff::C3,
      xtool_ff::C3x,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.305f,
      0.005f,
      0,
      "CP"
    },

    {
      "HA",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CB",
      xtool_ff::C3,
      xtool_ff::C3,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::H,
      -0.008f,
      0.021f,
      0,
      "CF"
    },

    {
      "HB1",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HB2",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CG",
      xtool_ff::C2,
      xtool_ff::C2x,
      1.900f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.030f,
      0.023f,
      0,
      "CO"
    },

    {
      "OD1",
      xtool_ff::Oco2,
      xtool_ff::Oco2,
      1.550f,
      0.000f,
      -0.500f,
      x_score_hb::A,
      -0.880f,
      0.000f,
      0,
      "OC"
    },

    {
      "OD2",
      xtool_ff::Oco2,
      xtool_ff::Oco2,
      1.550f,
      0.000f,
      -0.500f,
      x_score_hb::A,
      -0.880f,
      0.000f,
      0,
      "OC"
    }

  };

  static const xtool_residue_bond_description ASP_bonds[] = {

    { "N", "H", 1 },

    { "N", "CA", 1 },

    { "C", "O", 2 },

    { "C", "CA", 1 },

    { "CA", "HA", 1 },

    { "CA", "CB", 1 },

    { "CB", "HB1", 1 },

    { "CB", "HB2", 1 },

    { "CB", "CG", 1 },

    { "CG", "OD1", 2 },

    { "CG", "OD2", 2 }

  };

  static const xtool_residue_atom_description CYS_atoms[] = {

    {
      "N",
      xtool_ff::Nam,
      xtool_ff::Npl3h,
      1.750f,
      0.000f,
      0.000f,
      x_score_hb::D,
      -0.096f,
      -0.052f,
      0,
      "ND"
    },

    {
      "H",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      -0.115f,
      0,
      "HH"
    },

    {
      "C",
      xtool_ff::C2,
      xtool_ff::C2x,
      1.900f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.030f,
      0.023f,
      0,
      "CP"
    },

    {
      "O",
      xtool_ff::O2,
      xtool_ff::O2,
      1.550f,
      0.000f,
      0.000f,
      x_score_hb::A,
      -0.399f,
      -0.035f,
      0,
      "OA"
    },

    {
      "CA",
      xtool_ff::C3,
      xtool_ff::C3x,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.305f,
      0.005f,
      0,
      "CP"
    },

    {
      "HA",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CB",
      xtool_ff::C3,
      xtool_ff::C3,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.358f,
      0.005f,
      0,
      "CP"
    },

    {
      "HB1",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HB2",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "SG",
      xtool_ff::S3,
      xtool_ff::S3,
      2.100f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.419f,
      0.037f,
      0,
      "SA"
    },

    {
      "HG",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      -0.106f,
      0,
      "HH"
    }

  };

  static const xtool_residue_bond_description CYS_bonds[] = {

    { "N", "H", 1 },

    { "N", "CA", 1 },

    { "C", "O", 2 },

    { "C", "CA", 1 },

    { "CA", "HA", 1 },

    { "CA", "CB", 1 },

    { "CB", "HB1", 1 },

    { "CB", "HB2", 1 },

    { "CB", "SG", 1 },

    { "SG", "HG", 1 }

  };

  static const xtool_residue_atom_description GLN_atoms[] = {

    {
      "N",
      xtool_ff::Nam,
      xtool_ff::Npl3h,
      1.750f,
      0.000f,
      0.000f,
      x_score_hb::D,
      -0.096f,
      -0.052f,
      0,
      "ND"
    },

    {
      "H",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      -0.115f,
      0,
      "HH"
    },

    {
      "C",
      xtool_ff::C2,
      xtool_ff::C2x,
      1.900f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.030f,
      0.023f,
      0,
      "CP"
    },

    {
      "O",
      xtool_ff::O2,
      xtool_ff::O2,
      1.550f,
      0.000f,
      0.000f,
      x_score_hb::A,
      -0.399f,
      -0.035f,
      0,
      "OA"
    },

    {
      "CA",
      xtool_ff::C3,
      xtool_ff::C3x,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.305f,
      0.005f,
      0,
      "CP"
    },

    {
      "HA",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CB",
      xtool_ff::C3,
      xtool_ff::C3,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.358f,
      0.021f,
      0,
      "CF"
    },

    {
      "HB1",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HB2",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CG",
      xtool_ff::C3,
      xtool_ff::C3,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::H,
      -0.008f,
      0.021f,
      0,
      "CF"
    },

    {
      "HG1",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HG2",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CD",
      xtool_ff::C2,
      xtool_ff::C2x,
      1.900f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.030f,
      0.023f,
      0,
      "CP"
    },

    {
      "OE1",
      xtool_ff::O2,
      xtool_ff::O2,
      1.550f,
      0.000f,
      0.000f,
      x_score_hb::A,
      -0.399f,
      -0.035f,
      0,
      "OA"
    },

    {
      "NE2",
      xtool_ff::Nam,
      xtool_ff::Npl3h,
      1.750f,
      0.000f,
      0.000f,
      x_score_hb::D,
      -0.646f,
      -0.052f,
      0,
      "ND"
    },

    {
      "HE21",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      -0.115f,
      0,
      "HH"
    },

    {
      "HE22",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      -0.115f,
      0,
      "HH"
    }

  };

  static const xtool_residue_bond_description GLN_bonds[] = {

    { "N", "H", 1 },

    { "N", "CA", 1 },

    { "C", "O", 2 },

    { "C", "CA", 1 },

    { "CA", "HA", 1 },

    { "CA", "CB", 1 },

    { "CB", "HB1", 1 },

    { "CB", "HB2", 1 },

    { "CB", "CG", 1 },

    { "CG", "HG1", 1 },

    { "CG", "HG2", 1 },

    { "CG", "CD", 1 },

    { "CD", "OE1", 2 },

    { "CD", "NE2", 3 },

    { "NE2", "HE21", 1 },

    { "NE2", "HE22", 1 }

  };

  static const xtool_residue_atom_description GLU_atoms[] = {

    {
      "N",
      xtool_ff::Nam,
      xtool_ff::Npl3h,
      1.750f,
      0.000f,
      0.000f,
      x_score_hb::D,
      -0.096f,
      -0.052f,
      0,
      "ND"
    },

    {
      "H",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      -0.115f,
      0,
      "HH"
    },

    {
      "C",
      xtool_ff::C2,
      xtool_ff::C2x,
      1.900f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.030f,
      0.023f,
      0,
      "CP"
    },

    {
      "O",
      xtool_ff::O2,
      xtool_ff::O2,
      1.550f,
      0.000f,
      0.000f,
      x_score_hb::A,
      -0.399f,
      -0.035f,
      0,
      "OA"
    },

    {
      "CA",
      xtool_ff::C3,
      xtool_ff::C3x,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.305f,
      0.005f,
      0,
      "CP"
    },

    {
      "HA",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CB",
      xtool_ff::C3,
      xtool_ff::C3,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.358f,
      0.021f,
      0,
      "CF"
    },

    {
      "HB1",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HB2",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CG",
      xtool_ff::C3,
      xtool_ff::C3,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::H,
      -0.008f,
      0.021f,
      0,
      "CF"
    },

    {
      "HG1",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HG2",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CD",
      xtool_ff::C2,
      xtool_ff::C2x,
      1.900f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.030f,
      0.023f,
      0,
      "CO"
    },

    {
      "OE1",
      xtool_ff::Oco2,
      xtool_ff::Oco2,
      1.550f,
      0.000f,
      -0.500f,
      x_score_hb::A,
      -0.880f,
      0.000f,
      0,
      "OC"
    },

    {
      "OE2",
      xtool_ff::Oco2,
      xtool_ff::Oco2,
      1.550f,
      0.000f,
      -0.500f,
      x_score_hb::A,
      -0.880f,
      0.000f,
      0,
      "OC"
    }

  };

  static const xtool_residue_bond_description GLU_bonds[] = {

    { "N", "H", 1 },

    { "N", "CA", 1 },

    { "C", "O", 2 },

    { "C", "CA", 1 },

    { "CA", "HA", 1 },

    { "CA", "CB", 1 },

    { "CB", "HB1", 1 },

    { "CB", "HB2", 1 },

    { "CB", "CG", 1 },

    { "CG", "HG1", 1 },

    { "CG", "HG2", 1 },

    { "CG", "CD", 1 },

    { "CD", "OE1", 2 },

    { "CD", "OE2", 2 }

  };

  static const xtool_residue_atom_description GLY_atoms[] = {

    {
      "N",
      xtool_ff::Nam,
      xtool_ff::Npl3h,
      1.750f,
      0.000f,
      0.000f,
      x_score_hb::D,
      -0.096f,
      -0.052f,
      0,
      "ND"
    },

    {
      "H",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      -0.115f,
      0,
      "HH"
    },

    {
      "C",
      xtool_ff::C2,
      xtool_ff::C2x,
      1.900f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.030f,
      0.023f,
      0,
      "CP"
    },

    {
      "O",
      xtool_ff::O2,
      xtool_ff::O2,
      1.550f,
      0.000f,
      0.000f,
      x_score_hb::A,
      -0.399f,
      -0.035f,
      0,
      "OA"
    },

    {
      "CA",
      xtool_ff::C3,
      xtool_ff::C3x,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.303f,
      0.005f,
      0,
      "CP"
    },

    {
      "HA1",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HA2",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    }

  };

  static const xtool_residue_bond_description GLY_bonds[] = {

    { "N", "H", 1 },

    { "N", "CA", 1 },

    { "C", "O", 2 },

    { "C", "CA", 1 },

    { "CA", "HA1", 1 },

    { "CA", "HA2", 1 }

  };

  static const xtool_residue_atom_description HIS_atoms[] = {

    {
      "N",
      xtool_ff::Nam,
      xtool_ff::Npl3h,
      1.750f,
      0.000f,
      0.000f,
      x_score_hb::D,
      -0.096f,
      -0.052f,
      0,
      "ND"
    },

    {
      "H",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      -0.115f,
      0,
      "HH"
    },

    {
      "C",
      xtool_ff::C2,
      xtool_ff::C2x,
      1.900f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.030f,
      0.023f,
      0,
      "CP"
    },

    {
      "O",
      xtool_ff::O2,
      xtool_ff::O2,
      1.550f,
      0.000f,
      0.000f,
      x_score_hb::A,
      -0.399f,
      -0.035f,
      0,
      "OA"
    },

    {
      "CA",
      xtool_ff::C3,
      xtool_ff::C3x,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.305f,
      0.005f,
      0,
      "CP"
    },

    {
      "HA",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CB",
      xtool_ff::C3,
      xtool_ff::C3,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::H,
      -0.008f,
      0.021f,
      0,
      "CF"
    },

    {
      "HB1",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HB2",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CG",
      xtool_ff::C2,
      xtool_ff::Carx,
      2.000f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.027f,
      0.017f,
      2,
      "cP"
    },

    {
      "ND1",
      xtool_ff::N2,
      xtool_ff::Narh,
      1.750f,
      0.000f,
      0.500f,
      x_score_hb::DA,
      0.135f,
      -0.078f,
      2,
      "NR"
    },

    {
      "HD1",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      -0.115f,
      0,
      "HH"
    },

    {
      "CD2",
      xtool_ff::C2,
      xtool_ff::Carx,
      2.000f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.310f,
      0.017f,
      2,
      "cP"
    },

    {
      "HD2",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CE1",
      xtool_ff::C2,
      xtool_ff::Carx,
      2.000f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.310f,
      0.017f,
      2,
      "cP"
    },

    {
      "HE1",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "NE2",
      xtool_ff::N2,
      xtool_ff::Narh,
      1.750f,
      0.000f,
      0.500f,
      x_score_hb::DA,
      0.135f,
      -0.078f,
      2,
      "NR"
    },

    {
      "HE2",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      -0.115f,
      0,
      "HH"
    }

  };

  static const xtool_residue_bond_description HIS_bonds[] = {

    { "N", "H", 1 },

    { "N", "CA", 1 },

    { "C", "O", 2 },

    { "C", "CA", 1 },

    { "CA", "HA", 1 },

    { "CA", "CB", 1 },

    { "CB", "HB1", 1 },

    { "CB", "HB2", 1 },

    { "CB", "CG", 1 },

    { "CG", "ND1", 4 },

    { "ND1", "HD1", 1 },

    { "CG", "CD2", 4 },

    { "CD2", "HD2", 1 },

    { "ND1", "CE1", 4 },

    { "CE1", "HE1", 1 },

    { "CD2", "NE2", 4 },

    { "NE2", "HE2", 1 },

    { "CE1", "NE2", 4 }

  };

  static const xtool_residue_atom_description ILE_atoms[] = {

    {
      "N",
      xtool_ff::Nam,
      xtool_ff::Npl3h,
      1.750f,
      0.000f,
      0.000f,
      x_score_hb::D,
      -0.096f,
      -0.052f,
      0,
      "ND"
    },

    {
      "H",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      -0.115f,
      0,
      "HH"
    },

    {
      "C",
      xtool_ff::C2,
      xtool_ff::C2x,
      1.900f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.030f,
      0.023f,
      0,
      "CP"
    },

    {
      "O",
      xtool_ff::O2,
      xtool_ff::O2,
      1.550f,
      0.000f,
      0.000f,
      x_score_hb::A,
      -0.399f,
      -0.035f,
      0,
      "OA"
    },

    {
      "CA",
      xtool_ff::C3,
      xtool_ff::C3x,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.305f,
      0.005f,
      0,
      "CP"
    },

    {
      "HA",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CB",
      xtool_ff::C3,
      xtool_ff::C3,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.127f,
      0.021f,
      0,
      "CF"
    },

    {
      "HB",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CG2",
      xtool_ff::C3,
      xtool_ff::C3,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.528f,
      0.031f,
      0,
      "CF"
    },

    {
      "HG21",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HG22",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HG23",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CG1",
      xtool_ff::C3,
      xtool_ff::C3,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.358f,
      0.031f,
      0,
      "CF"
    },

    {
      "HG11",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HG12",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CD1",
      xtool_ff::C3,
      xtool_ff::C3,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.528f,
      0.031f,
      0,
      "CF"
    },

    {
      "HD11",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HD12",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HD13",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    }

  };

  static const xtool_residue_bond_description ILE_bonds[] = {

    { "N", "H", 1 },

    { "N", "CA", 1 },

    { "C", "O", 2 },

    { "C", "CA", 1 },

    { "CA", "HA", 1 },

    { "CA", "CB", 1 },

    { "CB", "HB", 1 },

    { "CB", "CG2", 1 },

    { "CG2", "HG21", 1 },

    { "CG2", "HG22", 1 },

    { "CG2", "HG23", 1 },

    { "CB", "CG1", 1 },

    { "CG1", "HG11", 1 },

    { "CG1", "HG12", 1 },

    { "CG1", "CD1", 1 },

    { "CD1", "HD11", 1 },

    { "CD1", "HD12", 1 },

    { "CD1", "HD13", 1 }

  };

  static const xtool_residue_atom_description LEU_atoms[] = {

    {
      "N",
      xtool_ff::Nam,
      xtool_ff::Npl3h,
      1.750f,
      0.000f,
      0.000f,
      x_score_hb::D,
      -0.096f,
      -0.052f,
      0,
      "ND"
    },

    {
      "H",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      -0.115f,
      0,
      "HH"
    },

    {
      "C",
      xtool_ff::C2,
      xtool_ff::C2x,
      1.900f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.030f,
      0.023f,
      0,
      "CP"
    },

    {
      "O",
      xtool_ff::O2,
      xtool_ff::O2,
      1.550f,
      0.000f,
      0.000f,
      x_score_hb::A,
      -0.399f,
      -0.035f,
      0,
      "OA"
    },

    {
      "CA",
      xtool_ff::C3,
      xtool_ff::C3x,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.305f,
      0.005f,
      0,
      "CP"
    },

    {
      "HA",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CB",
      xtool_ff::C3,
      xtool_ff::C3,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.358f,
      0.021f,
      0,
      "CF"
    },

    {
      "HB1",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HB2",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CG",
      xtool_ff::C3,
      xtool_ff::C3,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.127f,
      0.031f,
      0,
      "CF"
    },

    {
      "HG",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CD1",
      xtool_ff::C3,
      xtool_ff::C3,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.528f,
      0.031f,
      0,
      "CF"
    },

    {
      "HD11",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HD12",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HD13",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CD2",
      xtool_ff::C3,
      xtool_ff::C3,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.528f,
      0.031f,
      0,
      "CF"
    },

    {
      "HD21",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HD22",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HD23",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    }

  };

  static const xtool_residue_bond_description LEU_bonds[] = {

    { "N", "H", 1 },

    { "N", "CA", 1 },

    { "C", "O", 2 },

    { "C", "CA", 1 },

    { "CA", "HA", 1 },

    { "CA", "CB", 1 },

    { "CB", "HB1", 1 },

    { "CB", "HB2", 1 },

    { "CB", "CG", 1 },

    { "CG", "HG", 1 },

    { "CG", "CD1", 1 },

    { "CD1", "HD11", 1 },

    { "CD1", "HD12", 1 },

    { "CD1", "HD13", 1 },

    { "CG", "CD2", 1 },

    { "CD2", "HD21", 1 },

    { "CD2", "HD22", 1 },

    { "CD2", "HD23", 1 }

  };

  static const xtool_residue_atom_description LYS_atoms[] = {

    {
      "N",
      xtool_ff::Nam,
      xtool_ff::Npl3h,
      1.750f,
      0.000f,
      0.000f,
      x_score_hb::D,
      -0.096f,
      -0.052f,
      0,
      "ND"
    },

    {
      "H",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      -0.115f,
      0,
      "HH"
    },

    {
      "C",
      xtool_ff::C2,
      xtool_ff::C2x,
      1.900f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.030f,
      0.023f,
      0,
      "CP"
    },

    {
      "O",
      xtool_ff::O2,
      xtool_ff::O2,
      1.550f,
      0.000f,
      0.000f,
      x_score_hb::A,
      -0.399f,
      -0.035f,
      0,
      "OA"
    },

    {
      "CA",
      xtool_ff::C3,
      xtool_ff::C3x,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.305f,
      0.005f,
      0,
      "CP"
    },

    {
      "HA",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CB",
      xtool_ff::C3,
      xtool_ff::C3,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.358f,
      0.021f,
      0,
      "CF"
    },

    {
      "HB1",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HB2",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CG",
      xtool_ff::C3,
      xtool_ff::C3,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.358f,
      0.021f,
      0,
      "CF"
    },

    {
      "HG1",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HG2",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CD",
      xtool_ff::C3,
      xtool_ff::C3,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.358f,
      0.021f,
      0,
      "CF"
    },

    {
      "HD1",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HD2",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CE",
      xtool_ff::C3,
      xtool_ff::C3x,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.137f,
      0.005f,
      0,
      "CN"
    },

    {
      "HE1",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HE2",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "NZ",
      xtool_ff::N4,
      xtool_ff::N4,
      1.800f,
      0.000f,
      1.000f,
      x_score_hb::D,
      -1.200f,
      0.041f,
      0,
      "NC"
    },

    {
      "HZ1",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HZ2",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HZ3",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      0.000f,
      0,
      "HH"
    }

  };

  static const xtool_residue_bond_description LYS_bonds[] = {

    { "N", "H", 1 },

    { "N", "CA", 1 },

    { "C", "O", 2 },

    { "C", "CA", 1 },

    { "CA", "HA", 1 },

    { "CA", "CB", 1 },

    { "CB", "HB1", 1 },

    { "CB", "HB2", 1 },

    { "CB", "CG", 1 },

    { "CG", "HG1", 1 },

    { "CG", "HG2", 1 },

    { "CG", "CD", 1 },

    { "CD", "HD1", 1 },

    { "CD", "HD2", 1 },

    { "CD", "CE", 1 },

    { "CE", "HE1", 1 },

    { "CE", "HE2", 1 },

    { "CE", "NZ", 1 },

    { "NZ", "HZ1", 1 },

    { "NZ", "HZ2", 1 },

    { "NZ", "HZ3", 1 }

  };

  static const xtool_residue_atom_description MET_atoms[] = {

    {
      "N",
      xtool_ff::Nam,
      xtool_ff::Npl3h,
      1.750f,
      0.000f,
      0.000f,
      x_score_hb::D,
      -0.096f,
      -0.052f,
      0,
      "ND"
    },

    {
      "H",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      -0.115f,
      0,
      "HH"
    },

    {
      "C",
      xtool_ff::C2,
      xtool_ff::C2x,
      1.900f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.030f,
      0.023f,
      0,
      "CP"
    },

    {
      "O",
      xtool_ff::O2,
      xtool_ff::O2,
      1.550f,
      0.000f,
      0.000f,
      x_score_hb::A,
      -0.399f,
      -0.035f,
      0,
      "OA"
    },

    {
      "CA",
      xtool_ff::C3,
      xtool_ff::C3x,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.305f,
      0.005f,
      0,
      "CP"
    },

    {
      "HA",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CB",
      xtool_ff::C3,
      xtool_ff::C3,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.358f,
      0.021f,
      0,
      "CF"
    },

    {
      "HB1",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HB2",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CG",
      xtool_ff::C3,
      xtool_ff::C3,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.358f,
      0.005f,
      0,
      "CP"
    },

    {
      "HG1",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HG2",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "SD",
      xtool_ff::S3,
      xtool_ff::S3,
      2.100f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.255f,
      0.037f,
      0,
      "SA"
    },

    {
      "CE",
      xtool_ff::C3,
      xtool_ff::C3,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.528f,
      0.005f,
      0,
      "CP"
    },

    {
      "HE1",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HE2",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HE3",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    }

  };

  static const xtool_residue_bond_description MET_bonds[] = {

    { "N", "H", 1 },

    { "N", "CA", 1 },

    { "C", "O", 2 },

    { "C", "CA", 1 },

    { "CA", "HA", 1 },

    { "CA", "CB", 1 },

    { "CB", "HB1", 1 },

    { "CB", "HB2", 1 },

    { "CB", "CG", 1 },

    { "CG", "HG1", 1 },

    { "CG", "HG2", 1 },

    { "CG", "SD", 1 },

    { "SD", "CE", 1 },

    { "CE", "HE1", 1 },

    { "CE", "HE2", 1 },

    { "CE", "HE3", 1 }

  };

  static const xtool_residue_atom_description NME_atoms[] = {

    {
      "N",
      xtool_ff::Nam,
      xtool_ff::Npl3h,
      1.750f,
      0.000f,
      0.000f,
      x_score_hb::D,
      -0.096f,
      -0.052f,
      0,
      "ND"
    },

    {
      "H",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      -0.115f,
      0,
      "HH"
    },

    {
      "CA",
      xtool_ff::C3,
      xtool_ff::C3x,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.032f,
      0.005f,
      0,
      "CP"
    },

    {
      "HA1",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HA2",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HA3",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    }

  };

  static const xtool_residue_bond_description NME_bonds[] = {

    { "N", "H", 1 },

    { "N", "CA", 1 },

    { "CA", "HA1", 1 },

    { "CA", "HA2", 1 },

    { "CA", "HA3", 1 }

  };

  static const xtool_residue_atom_description PHE_atoms[] = {

    {
      "N",
      xtool_ff::Nam,
      xtool_ff::Npl3h,
      1.750f,
      0.000f,
      0.000f,
      x_score_hb::D,
      -0.096f,
      -0.052f,
      0,
      "ND"
    },

    {
      "H",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      -0.115f,
      0,
      "HH"
    },

    {
      "C",
      xtool_ff::C2,
      xtool_ff::C2x,
      1.900f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.030f,
      0.023f,
      0,
      "CP"
    },

    {
      "O",
      xtool_ff::O2,
      xtool_ff::O2,
      1.550f,
      0.000f,
      0.000f,
      x_score_hb::A,
      -0.399f,
      -0.035f,
      0,
      "OA"
    },

    {
      "CA",
      xtool_ff::C3,
      xtool_ff::C3x,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.305f,
      0.005f,
      0,
      "CP"
    },

    {
      "HA",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CB",
      xtool_ff::C3,
      xtool_ff::C3,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::H,
      -0.008f,
      0.021f,
      0,
      "CF"
    },

    {
      "HB1",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HB2",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CG",
      xtool_ff::Car,
      xtool_ff::Car,
      2.000f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.296f,
      0.017f,
      2,
      "cF"
    },

    {
      "CD1",
      xtool_ff::Car,
      xtool_ff::Car,
      2.000f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.337f,
      0.017f,
      2,
      "cF"
    },

    {
      "HD1",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CD2",
      xtool_ff::Car,
      xtool_ff::Car,
      2.000f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.337f,
      0.017f,
      2,
      "cF"
    },

    {
      "HD2",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CE1",
      xtool_ff::Car,
      xtool_ff::Car,
      2.000f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.337f,
      0.017f,
      2,
      "cF"
    },

    {
      "HE1",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CE2",
      xtool_ff::Car,
      xtool_ff::Car,
      2.000f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.337f,
      0.017f,
      2,
      "cF"
    },

    {
      "HE2",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CZ",
      xtool_ff::Car,
      xtool_ff::Car,
      2.000f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.337f,
      0.017f,
      2,
      "cF"
    },

    {
      "HZ",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    }

  };

  static const xtool_residue_bond_description PHE_bonds[] = {

    { "N", "H", 1 },

    { "N", "CA", 1 },

    { "C", "O", 2 },

    { "C", "CA", 1 },

    { "CA", "HA", 1 },

    { "CA", "CB", 1 },

    { "CB", "HB1", 1 },

    { "CB", "HB2", 1 },

    { "CB", "CG", 1 },

    { "CG", "CD1", 4 },

    { "CD1", "HD1", 1 },

    { "CG", "CD2", 4 },

    { "CD2", "HD2", 1 },

    { "CD1", "CE1", 4 },

    { "CE1", "HE1", 1 },

    { "CD2", "CE2", 4 },

    { "CE2", "HE2", 1 },

    { "CE1", "CZ", 4 },

    { "CE2", "CZ", 4 },

    { "CZ", "HZ", 1 }

  };

  static const xtool_residue_atom_description PRO_atoms[] = {

    {
      "N",
      xtool_ff::Nam,
      xtool_ff::Npl3,
      1.750f,
      0.000f,
      0.000f,
      x_score_hb::P,
      0.078f,
      -0.052f,
      1,
      "ND"
    },

    {
      "C",
      xtool_ff::C2,
      xtool_ff::C2x,
      1.900f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.030f,
      0.023f,
      0,
      "CP"
    },

    {
      "O",
      xtool_ff::O2,
      xtool_ff::O2,
      1.550f,
      0.000f,
      0.000f,
      x_score_hb::A,
      -0.399f,
      -0.035f,
      0,
      "OA"
    },

    {
      "CD",
      xtool_ff::C3,
      xtool_ff::C3x,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.137f,
      0.005f,
      1,
      "CP"
    },

    {
      "HD1",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HD2",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CA",
      xtool_ff::C3,
      xtool_ff::C3x,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.305f,
      0.005f,
      1,
      "CP"
    },

    {
      "HA",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CB",
      xtool_ff::C3,
      xtool_ff::C3,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.358f,
      0.021f,
      1,
      "CF"
    },

    {
      "HB1",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HB2",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CG",
      xtool_ff::C3,
      xtool_ff::C3,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.358f,
      0.021f,
      1,
      "CF"
    },

    {
      "HG1",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HG2",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    }

  };

  static const xtool_residue_bond_description PRO_bonds[] = {

    { "N", "H", 1 },

    { "N", "CA", 1 },

    { "C", "O", 2 },

    { "C", "CA", 1 },

    { "CA", "HA", 1 },

    { "CA", "CB", 1 },

    { "CB", "HB1", 1 },

    { "CB", "HB2", 1 },

    { "CB", "CG", 1 },

    { "CG", "HG1", 1 },

    { "CG", "HG2", 1 },

    { "CG", "CD", 1 },

    { "CD", "HD1", 1 },

    { "CD", "HD2", 1 },

    { "N", "CD", 1 }

  };

  static const xtool_residue_atom_description SER_atoms[] = {

    {
      "N",
      xtool_ff::Nam,
      xtool_ff::Npl3h,
      1.750f,
      0.000f,
      0.000f,
      x_score_hb::D,
      -0.096f,
      -0.052f,
      0,
      "ND"
    },

    {
      "H",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      -0.115f,
      0,
      "HH"
    },

    {
      "C",
      xtool_ff::C2,
      xtool_ff::C2x,
      1.900f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.030f,
      0.023f,
      0,
      "CP"
    },

    {
      "O",
      xtool_ff::O2,
      xtool_ff::O2,
      1.550f,
      0.000f,
      0.000f,
      x_score_hb::A,
      -0.399f,
      -0.035f,
      0,
      "OA"
    },

    {
      "CA",
      xtool_ff::C3,
      xtool_ff::C3x,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.305f,
      0.005f,
      0,
      "CP"
    },

    {
      "HA",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CB",
      xtool_ff::C3,
      xtool_ff::C3x,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.137f,
      0.005f,
      0,
      "CP"
    },

    {
      "HB1",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HB2",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "OG",
      xtool_ff::O3,
      xtool_ff::O3h,
      1.650f,
      0.000f,
      0.000f,
      x_score_hb::DA,
      -0.467f,
      -0.014f,
      0,
      "OD"
    },

    {
      "HG",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      -0.115f,
      0,
      "HH"
    }

  };

  static const xtool_residue_bond_description SER_bonds[] = {

    { "N", "H", 1 },

    { "N", "CA", 1 },

    { "C", "O", 2 },

    { "C", "CA", 1 },

    { "CA", "HA", 1 },

    { "CA", "CB", 1 },

    { "CB", "HB1", 1 },

    { "CB", "HB2", 1 },

    { "CB", "OG", 1 },

    { "OG", "HG", 1 }

  };

  static const xtool_residue_atom_description THR_atoms[] = {

    {
      "N",
      xtool_ff::Nam,
      xtool_ff::Npl3h,
      1.750f,
      0.000f,
      0.000f,
      x_score_hb::D,
      -0.096f,
      -0.052f,
      0,
      "ND"
    },

    {
      "H",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      -0.115f,
      0,
      "HH"
    },

    {
      "C",
      xtool_ff::C2,
      xtool_ff::C2x,
      1.900f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.030f,
      0.023f,
      0,
      "CP"
    },

    {
      "O",
      xtool_ff::O2,
      xtool_ff::O2,
      1.550f,
      0.000f,
      0.000f,
      x_score_hb::A,
      -0.399f,
      -0.035f,
      0,
      "OA"
    },

    {
      "CA",
      xtool_ff::C3,
      xtool_ff::C3x,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.305f,
      0.005f,
      0,
      "CP"
    },

    {
      "HA",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CB",
      xtool_ff::C3,
      xtool_ff::C3x,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.205f,
      0.005f,
      0,
      "CP"
    },

    {
      "HB",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "OG1",
      xtool_ff::O3,
      xtool_ff::O3h,
      1.650f,
      0.000f,
      0.000f,
      x_score_hb::DA,
      -0.467f,
      -0.014f,
      0,
      "OD"
    },

    {
      "HG1",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      -0.115f,
      0,
      "HH"
    },

    {
      "CG2",
      xtool_ff::C3,
      xtool_ff::C3,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.528f,
      0.021f,
      0,
      "CF"
    },

    {
      "HG21",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HG22",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HG23",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    }

  };

  static const xtool_residue_bond_description THR_bonds[] = {

    { "N", "H", 1 },

    { "N", "CA", 1 },

    { "C", "O", 2 },

    { "C", "CA", 1 },

    { "CA", "HA", 1 },

    { "CA", "CB", 1 },

    { "CB", "HB", 1 },

    { "CB", "OG1", 1 },

    { "OG1", "HG1", 1 },

    { "CB", "CG2", 1 },

    { "CG2", "HG21", 1 },

    { "CG2", "HG22", 1 },

    { "CG2", "HG23", 1 }

  };

  static const xtool_residue_atom_description TRP_atoms[] = {

    {
      "N",
      xtool_ff::Nam,
      xtool_ff::Npl3h,
      1.750f,
      0.000f,
      0.000f,
      x_score_hb::D,
      -0.096f,
      -0.052f,
      0,
      "ND"
    },

    {
      "H",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      -0.115f,
      0,
      "HH"
    },

    {
      "C",
      xtool_ff::C2,
      xtool_ff::C2x,
      1.900f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.030f,
      0.023f,
      0,
      "CP"
    },

    {
      "O",
      xtool_ff::O2,
      xtool_ff::O2,
      1.550f,
      0.000f,
      0.000f,
      x_score_hb::A,
      -0.399f,
      -0.035f,
      0,
      "OA"
    },

    {
      "CA",
      xtool_ff::C3,
      xtool_ff::C3x,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.305f,
      0.005f,
      0,
      "CP"
    },

    {
      "HA",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CB",
      xtool_ff::C3,
      xtool_ff::C3,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::H,
      -0.008f,
      0.021f,
      0,
      "CF"
    },

    {
      "HB1",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HB2",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CG",
      xtool_ff::C2,
      xtool_ff::Car,
      2.000f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.013f,
      0.017f,
      2,
      "cF"
    },

    {
      "CD2",
      xtool_ff::Car,
      xtool_ff::Car,
      2.000f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.296f,
      0.017f,
      2,
      "cF"
    },

    {
      "CE2",
      xtool_ff::Car,
      xtool_ff::Carx,
      2.000f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.151f,
      0.017f,
      2,
      "cP"
    },

    {
      "CE3",
      xtool_ff::Car,
      xtool_ff::Car,
      2.000f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.337f,
      0.017f,
      2,
      "cF"
    },

    {
      "HE3",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CD1",
      xtool_ff::C2,
      xtool_ff::Carx,
      2.000f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.310f,
      0.017f,
      2,
      "cP"
    },

    {
      "HD1",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "NE1",
      xtool_ff::Npl3,
      xtool_ff::Narh,
      1.750f,
      0.000f,
      0.000f,
      x_score_hb::D,
      0.545f,
      -0.078f,
      2,
      "ND"
    },

    {
      "HE1",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      -0.115f,
      0,
      "HH"
    },

    {
      "CZ2",
      xtool_ff::Car,
      xtool_ff::Car,
      2.000f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.337f,
      0.017f,
      2,
      "cF"
    },

    {
      "HZ2",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CZ3",
      xtool_ff::Car,
      xtool_ff::Car,
      2.000f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.337f,
      0.017f,
      2,
      "cF"
    },

    {
      "HZ3",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CH2",
      xtool_ff::Car,
      xtool_ff::Car,
      2.000f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.337f,
      0.017f,
      2,
      "cF"
    },

    {
      "HH2",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    }

  };

  static const xtool_residue_bond_description TRP_bonds[] = {

    { "N", "H", 1 },

    { "N", "CA", 1 },

    { "C", "O", 2 },

    { "C", "CA", 1 },

    { "CA", "HA", 1 },

    { "CA", "CB", 1 },

    { "CB", "HB1", 1 },

    { "CB", "HB2", 1 },

    { "CB", "CG", 1 },

    { "CG", "CD1", 4 },

    { "CD1", "HD1", 1 },

    { "CG", "CD2", 4 },

    { "CD1", "NE1", 4 },

    { "NE1", "HE1", 1 },

    { "NE1", "CE2", 4 },

    { "CD2", "CE2", 4 },

    { "CD2", "CE3", 4 },

    { "CE3", "HE3", 1 },

    { "CE2", "CZ2", 4 },

    { "CZ2", "HZ2", 1 },

    { "CE3", "CZ3", 4 },

    { "CZ3", "HZ3", 1 },

    { "CZ2", "CH2", 4 },

    { "CZ3", "CH2", 4 },

    { "CH2", "HH2", 1 }

  };

  static const xtool_residue_atom_description TYR_atoms[] = {

    {
      "N",
      xtool_ff::Nam,
      xtool_ff::Npl3h,
      1.750f,
      0.000f,
      0.000f,
      x_score_hb::D,
      -0.096f,
      -0.052f,
      0,
      "ND"
    },

    {
      "H",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      -0.115f,
      0,
      "HH"
    },

    {
      "C",
      xtool_ff::C2,
      xtool_ff::C2x,
      1.900f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.030f,
      0.023f,
      0,
      "CP"
    },

    {
      "O",
      xtool_ff::O2,
      xtool_ff::O2,
      1.550f,
      0.000f,
      0.000f,
      x_score_hb::A,
      -0.399f,
      -0.035f,
      0,
      "OA"
    },

    {
      "CA",
      xtool_ff::C3,
      xtool_ff::C3x,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.305f,
      0.005f,
      0,
      "CP"
    },

    {
      "HA",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CB",
      xtool_ff::C3,
      xtool_ff::C3,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::H,
      -0.008f,
      0.021f,
      0,
      "CF"
    },

    {
      "HB1",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HB2",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CG",
      xtool_ff::Car,
      xtool_ff::Car,
      2.000f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.296f,
      0.017f,
      2,
      "cF"
    },

    {
      "CD1",
      xtool_ff::Car,
      xtool_ff::Car,
      2.000f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.337f,
      0.017f,
      2,
      "cF"
    },

    {
      "HD1",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CD2",
      xtool_ff::Car,
      xtool_ff::Car,
      2.000f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.337f,
      0.017f,
      2,
      "cF"
    },

    {
      "HD2",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CE1",
      xtool_ff::Car,
      xtool_ff::Car,
      2.000f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.337f,
      0.017f,
      2,
      "cF"
    },

    {
      "HE1",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CE2",
      xtool_ff::Car,
      xtool_ff::Car,
      2.000f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.337f,
      0.017f,
      2,
      "cF"
    },

    {
      "HE2",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CZ",
      xtool_ff::Car,
      xtool_ff::Carx,
      2.000f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.151f,
      0.017f,
      2,
      "cP"
    },

    {
      "OH",
      xtool_ff::O3,
      xtool_ff::O3h,
      1.650f,
      0.000f,
      0.000f,
      x_score_hb::DA,
      0.082f,
      -0.014f,
      0,
      "OD"
    },

    {
      "HH",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      -0.115f,
      0,
      "HH"
    }

  };

  static const xtool_residue_bond_description TYR_bonds[] = {

    { "N", "H", 1 },

    { "N", "CA", 1 },

    { "C", "O", 2 },

    { "C", "CA", 1 },

    { "CA", "HA", 1 },

    { "CA", "CB", 1 },

    { "CB", "HB1", 1 },

    { "CB", "HB2", 1 },

    { "CB", "CG", 1 },

    { "CG", "CD1", 4 },

    { "CD1", "HD1", 1 },

    { "CG", "CD2", 4 },

    { "CD2", "HD2", 1 },

    { "CD1", "CE1", 4 },

    { "CE1", "HE1", 1 },

    { "CD2", "CE2", 4 },

    { "CE2", "HE2", 1 },

    { "CE1", "CZ", 4 },

    { "CE2", "CZ", 4 },

    { "CZ", "OH", 1 },

    { "OH", "HH", 1 }

  };

  static const xtool_residue_atom_description VAL_atoms[] = {

    {
      "N",
      xtool_ff::Nam,
      xtool_ff::Npl3h,
      1.750f,
      0.000f,
      0.000f,
      x_score_hb::D,
      -0.096f,
      -0.052f,
      0,
      "ND"
    },

    {
      "H",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      -0.115f,
      0,
      "HH"
    },

    {
      "C",
      xtool_ff::C2,
      xtool_ff::C2x,
      1.900f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.030f,
      0.023f,
      0,
      "CP"
    },

    {
      "O",
      xtool_ff::O2,
      xtool_ff::O2,
      1.550f,
      0.000f,
      0.000f,
      x_score_hb::A,
      -0.399f,
      -0.035f,
      0,
      "OA"
    },

    {
      "CA",
      xtool_ff::C3,
      xtool_ff::C3x,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::P,
      -0.305f,
      0.005f,
      0,
      "CP"
    },

    {
      "HA",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CB",
      xtool_ff::C3,
      xtool_ff::C3,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.127f,
      0.021f,
      0,
      "CF"
    },

    {
      "HB",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CG1",
      xtool_ff::C3,
      xtool_ff::C3,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.528f,
      0.031f,
      0,
      "CF"
    },

    {
      "HG11",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HG12",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HG13",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "CG2",
      xtool_ff::C3,
      xtool_ff::C3,
      2.200f,
      0.000f,
      0.000f,
      x_score_hb::H,
      0.528f,
      0.031f,
      0,
      "CF"
    },

    {
      "HG21",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HG22",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    },

    {
      "HG23",
      xtool_ff::H,
      xtool_ff::H,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::N,
      0.000f,
      0.000f,
      0,
      "HH"
    }

  };

  static const xtool_residue_bond_description VAL_bonds[] = {

    { "N", "H", 1 },

    { "N", "CA", 1 },

    { "C", "O", 2 },

    { "C", "CA", 1 },

    { "CA", "HA", 1 },

    { "CA", "CB", 1 },

    { "CB", "HB", 1 },

    { "CB", "CG1", 1 },

    { "CG1", "HG11", 1 },

    { "CG1", "HG12", 1 },

    { "CG1", "HG13", 1 },

    { "CB", "CG2", 1 },

    { "CG2", "HG21", 1 },

    { "CG2", "HG22", 1 },

    { "CG2", "HG23", 1 }

  };

  static const xtool_residue_atom_description TER_atoms[] = {

    {
      "OXT",
      xtool_ff::Oco2,
      xtool_ff::Oco2,
      1.550f,
      0.000f,
      -0.500f,
      x_score_hb::A,
      -0.399f,
      0.000f,
      0,
      "OC"
    },

    {
      "HOCA",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      -0.115f,
      0,
      "HH"
    },

    {
      "HN1",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      -0.115f,
      0,
      "HH"
    },

    {
      "HN2",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      -0.115f,
      0,
      "HH"
    },

    {
      "HN3",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      -0.115f,
      0,
      "HH"
    }

  };

  // declared as span because no bonds are present in the .dat file
  static const std::span<const xtool_residue_bond_description> TER_bonds;

  static const xtool_residue_atom_description HOH_atoms[] = {

    {
      "O",
      xtool_ff::Ow,
      xtool_ff::Ow,
      1.750f,
      0.000f,
      0.000f,
      x_score_hb::DA,
      -1.000f,
      -0.014f,
      0,
      "OW"
    },

    {
      "H1",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      -0.115f,
      0,
      "HH"
    },

    {
      "H2",
      xtool_ff::H,
      xtool_ff::Hhb,
      1.000f,
      0.000f,
      0.000f,
      x_score_hb::DH,
      0.000f,
      -0.115f,
      0,
      "HH"
    }

  };

  static const xtool_residue_bond_description HOH_bonds[] = {

    { "O", "H1", 1 },

    { "O", "H2", 1 }

  };

  static const xtool_residue_atom_description SO4_atoms[] = {

    {
      "S",
      xtool_ff::So2,
      xtool_ff::So,
      2.000f,
      0.000f,
      0.000f,
      x_score_hb::P,
      0.000f,
      0.000f,
      0,
      "UN"
    },

    {
      "O1",
      xtool_ff::Oco2,
      xtool_ff::Oco2,
      1.550f,
      0.000f,
      -0.500f,
      x_score_hb::A,
      -0.399f,
      0.000f,
      0,
      "OC"
    },

    {
      "O2",
      xtool_ff::Oco2,
      xtool_ff::Oco2,
      1.550f,
      0.000f,
      -0.500f,
      x_score_hb::A,
      -0.399f,
      0.000f,
      0,
      "OC"
    },

    {
      "O3",
      xtool_ff::Oco2,
      xtool_ff::Oco2,
      1.550f,
      0.000f,
      -0.500f,
      x_score_hb::A,
      -0.399f,
      0.000f,
      0,
      "OC"
    },

    {
      "O4",
      xtool_ff::Oco2,
      xtool_ff::Oco2,
      1.550f,
      0.000f,
      -0.500f,
      x_score_hb::A,
      -0.399f,
      0.000f,
      0,
      "OC"
    }

  };

  static const xtool_residue_bond_description SO4_bonds[] = {

    { "S", "O1", 2 },

    { "S", "O2", 2 },

    { "S", "O3", 2 },

    { "S", "O4", 2 }

  };

  static const xtool_residue_atom_description PO4_atoms[] = {

    {
      "P",
      xtool_ff::P3,
      xtool_ff::P3,
      2.000f,
      0.000f,
      0.000f,
      x_score_hb::P,
      0.000f,
      0.000f,
      0,
      "UN"
    },

    {
      "O1",
      xtool_ff::Oco2,
      xtool_ff::Oco2,
      1.550f,
      0.000f,
      -0.500f,
      x_score_hb::A,
      -0.399f,
      0.000f,
      0,
      "OC"
    },

    {
      "O2",
      xtool_ff::Oco2,
      xtool_ff::Oco2,
      1.550f,
      0.000f,
      -0.500f,
      x_score_hb::A,
      -0.399f,
      0.000f,
      0,
      "OC"
    },

    {
      "O3",
      xtool_ff::Oco2,
      xtool_ff::Oco2,
      1.550f,
      0.000f,
      -0.500f,
      x_score_hb::A,
      -0.399f,
      0.000f,
      0,
      "OC"
    },

    {
      "O4",
      xtool_ff::Oco2,
      xtool_ff::Oco2,
      1.550f,
      0.000f,
      -0.500f,
      x_score_hb::A,
      -0.399f,
      0.000f,
      0,
      "OC"
    }

  };

  static const xtool_residue_bond_description PO4_bonds[] = {

    { "P", "O1", 2 },

    { "P", "O2", 2 },

    { "P", "O3", 2 },

    { "P", "O4", 2 }

  };

  static const xtool_residue_atom_description HET_atoms[] = {

    {
      "LI",
      xtool_ff::Li,
      xtool_ff::Mplus,
      1.250f,
      0.000f,
      1.000f,
      x_score_hb::M,
      -2.000f,
      0.000f,
      0,
      "UN"
    },

    {
      "NA",
      xtool_ff::Na,
      xtool_ff::Mplus,
      1.250f,
      0.000f,
      1.000f,
      x_score_hb::M,
      -2.000f,
      0.000f,
      0,
      "UN"
    },

    {
      "K",
      xtool_ff::K,
      xtool_ff::Mplus,
      1.250f,
      0.000f,
      1.000f,
      x_score_hb::M,
      -2.000f,
      0.000f,
      0,
      "UN"
    },

    {
      "CA",
      xtool_ff::Ca,
      xtool_ff::Mplus,
      1.250f,
      0.000f,
      2.000f,
      x_score_hb::M,
      -2.000f,
      0.000f,
      0,
      "UN"
    },

    {
      "MG",
      xtool_ff::Mg,
      xtool_ff::Mplus,
      1.250f,
      0.000f,
      2.000f,
      x_score_hb::M,
      -2.000f,
      0.000f,
      0,
      "UN"
    },

    {
      "AL",
      xtool_ff::Al,
      xtool_ff::Mplus,
      1.250f,
      0.000f,
      3.000f,
      x_score_hb::M,
      -2.000f,
      0.000f,
      0,
      "UN"
    },

    {
      "MN",
      xtool_ff::Mn,
      xtool_ff::Mplus,
      1.250f,
      0.000f,
      2.000f,
      x_score_hb::M,
      -2.000f,
      0.000f,
      0,
      "UN"
    },

    {
      "FE",
      xtool_ff::Fe,
      xtool_ff::Mplus,
      1.250f,
      0.000f,
      2.000f,
      x_score_hb::M,
      -2.000f,
      0.000f,
      0,
      "UN"
    },

    {
      "NI",
      xtool_ff::Ni,
      xtool_ff::Mplus,
      1.250f,
      0.000f,
      2.000f,
      x_score_hb::M,
      -2.000f,
      0.000f,
      0,
      "UN"
    },

    {
      "CD",
      xtool_ff::Cd,
      xtool_ff::Mplus,
      1.250f,
      0.000f,
      2.000f,
      x_score_hb::M,
      -2.000f,
      0.000f,
      0,
      "UN"
    },

    {
      "CO",
      xtool_ff::Co,
      xtool_ff::Mplus,
      1.250f,
      0.000f,
      2.000f,
      x_score_hb::M,
      -2.000f,
      0.000f,
      0,
      "UN"
    },

    {
      "CU",
      xtool_ff::Cu,
      xtool_ff::Mplus,
      1.250f,
      0.000f,
      2.000f,
      x_score_hb::M,
      -2.000f,
      0.000f,
      0,
      "UN"
    },

    {
      "ZN",
      xtool_ff::Zn,
      xtool_ff::Mplus,
      1.250f,
      0.000f,
      2.000f,
      x_score_hb::M,
      -2.000f,
      0.000f,
      0,
      "UN"
    },

    {
      "HG",
      xtool_ff::Hg,
      xtool_ff::Mplus,
      1.250f,
      0.000f,
      2.000f,
      x_score_hb::M,
      -2.000f,
      0.000f,
      0,
      "UN"
    },

    {
      "U",
      xtool_ff::U,
      xtool_ff::Mplus,
      1.250f,
      0.000f,
      3.000f,
      x_score_hb::M,
      -2.000f,
      0.000f,
      0,
      "UN"
    },

    {
      "F",
      xtool_ff::F,
      xtool_ff::Fminus,
      1.500f,
      0.000f,
      -1.000f,
      x_score_hb::P,
      -1.000f,
      0.000f,
      0,
      "UN"
    },

    {
      "CL",
      xtool_ff::Cl,
      xtool_ff::Clminus,
      1.750f,
      0.000f,
      -1.000f,
      x_score_hb::P,
      -1.000f,
      0.000f,
      0,
      "UN"
    },

    {
      "BR",
      xtool_ff::Br,
      xtool_ff::Brminus,
      1.900f,
      0.000f,
      -1.000f,
      x_score_hb::P,
      -1.000f,
      0.000f,
      0,
      "UN"
    },

    {
      "I",
      xtool_ff::I,
      xtool_ff::Iminus,
      2.050f,
      0.000f,
      -1.000f,
      x_score_hb::P,
      -1.000f,
      0.000f,
      0,
      "UN"
    }

  };

  // declared as span because no bonds are present in the .dat file
  static const std::span<const xtool_residue_bond_description> HET_bonds;


  // --- Main Dictionary ---
  const std::array<xtool_residue_description, num_residues()> XTOOL_RESIDUE_DICTIONARY = {{

    {
      residue::ACE,
      "ACE",
      0.00f,
      "Acetyl",
      ACE_atoms,
      ACE_bonds
    },

    {
      residue::ALA,
      "ALA",
      0.00f,
      "Alanine",
      ALA_atoms,
      ALA_bonds
    },

    {
      residue::ARG,
      "ARG",
      1.00f,
      "Arginine",
      ARG_atoms,
      ARG_bonds
    },

    {
      residue::ASN,
      "ASN",
      0.00f,
      "Asparagine",
      ASN_atoms,
      ASN_bonds
    },

    {
      residue::ASP,
      "ASP",
      -1.00f,
      "Aspartic acid",
      ASP_atoms,
      ASP_bonds
    },

    {
      residue::CYS,
      "CYS",
      0.00f,
      "Cysteine",
      CYS_atoms,
      CYS_bonds
    },

    {
      residue::GLN,
      "GLN",
      0.00f,
      "Glutamine",
      GLN_atoms,
      GLN_bonds
    },

    {
      residue::GLU,
      "GLU",
      -1.00f,
      "Glutamic acid",
      GLU_atoms,
      GLU_bonds
    },

    {
      residue::GLY,
      "GLY",
      0.00f,
      "Glycine",
      GLY_atoms,
      GLY_bonds
    },

    {
      residue::HIS,
      "HIS",
      0.00f,
      "Histidine (all variations)",
      HIS_atoms,
      HIS_bonds
    },

    {
      residue::ILE,
      "ILE",
      0.00f,
      "Isoleucine",
      ILE_atoms,
      ILE_bonds
    },

    {
      residue::LEU,
      "LEU",
      0.00f,
      "Leucine",
      LEU_atoms,
      LEU_bonds
    },

    {
      residue::LYS,
      "LYS",
      1.00f,
      "Lysine",
      LYS_atoms,
      LYS_bonds
    },

    {
      residue::MET,
      "MET",
      0.00f,
      "Methionine",
      MET_atoms,
      MET_bonds
    },

    {
      residue::NME,
      "NME",
      0.00f,
      "N-Methyl",
      NME_atoms,
      NME_bonds
    },

    {
      residue::PHE,
      "PHE",
      0.00f,
      "Phenylalanine",
      PHE_atoms,
      PHE_bonds
    },

    {
      residue::PRO,
      "PRO",
      0.00f,
      "Proline",
      PRO_atoms,
      PRO_bonds
    },

    {
      residue::SER,
      "SER",
      0.00f,
      "Serine",
      SER_atoms,
      SER_bonds
    },

    {
      residue::THR,
      "THR",
      0.00f,
      "Threonine",
      THR_atoms,
      THR_bonds
    },

    {
      residue::TRP,
      "TRP",
      0.00f,
      "Tryptophan",
      TRP_atoms,
      TRP_bonds
    },

    {
      residue::TYR,
      "TYR",
      0.00f,
      "Tyrosine",
      TYR_atoms,
      TYR_bonds
    },

    {
      residue::VAL,
      "VAL",
      0.00f,
      "Valine",
      VAL_atoms,
      VAL_bonds
    },

    {
      residue::TER,
      "TER",
      0.00f,
      "N- and C-terminal atoms",
      TER_atoms,
      TER_bonds
    },

    {
      residue::HOH,
      "HOH",
      0.00f,
      "Water",
      HOH_atoms,
      HOH_bonds
    },

    {
      residue::SO4,
      "SO4",
      -2.00f,
      "SO4--",
      SO4_atoms,
      SO4_bonds
    },

    {
      residue::PO4,
      "PO4",
      -2.00f,
      "PO4--",
      PO4_atoms,
      PO4_bonds
    },

    {
      residue::HET,
      "HET",
      0.00f,
      "ions",
      HET_atoms,
      HET_bonds
    }

  }};

  // --- Lookup Map ---
  const std::unordered_map<std::string_view, residue> XTOOL_RESIDUE_LOOKUP = {

    { "ACE", residue::ACE },

    { "ALA", residue::ALA },

    { "ARG", residue::ARG },

    { "ASN", residue::ASN },

    { "ASP", residue::ASP },

    { "CYS", residue::CYS },

    { "GLN", residue::GLN },

    { "GLU", residue::GLU },

    { "GLY", residue::GLY },

    { "HIS", residue::HIS },

    { "ILE", residue::ILE },

    { "LEU", residue::LEU },

    { "LYS", residue::LYS },

    { "MET", residue::MET },

    { "NME", residue::NME },

    { "PHE", residue::PHE },

    { "PRO", residue::PRO },

    { "SER", residue::SER },

    { "THR", residue::THR },

    { "TRP", residue::TRP },

    { "TYR", residue::TYR },

    { "VAL", residue::VAL },

    { "TER", residue::TER },

    { "HOH", residue::HOH },

    { "SO4", residue::SO4 },

    { "PO4", residue::PO4 },

    { "HET", residue::HET }

  };

} // namespace mudock