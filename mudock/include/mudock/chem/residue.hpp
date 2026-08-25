#pragma once

#include <cassert>
#include <mudock/type_alias.hpp>
#include <string_view>

namespace mudock {

  // List of all known residues 
  enum class residue : int {

    ACE = 0, // Acetyl

    ALA = 1, // Alanine

    ARG = 2, // Arginine

    ASN = 3, // Asparagine

    ASP = 4, // Aspartic acid

    CYS = 5, // Cysteine

    GLN = 6, // Glutamine

    GLU = 7, // Glutamic acid

    GLY = 8, // Glycine

    HIS = 9, // Histidine (all variations)

    ILE = 10, // Isoleucine

    LEU = 11, // Leucine

    LYS = 12, // Lysine

    MET = 13, // Methionine

    NME = 14, // N-Methyl

    PHE = 15, // Phenylalanine

    PRO = 16, // Proline

    SER = 17, // Serine

    THR = 18, // Threonine

    TRP = 19, // Tryptophan

    TYR = 20, // Tyrosine

    VAL = 21, // Valine

    TER = 22, // N- and C-terminal atoms

    HOH = 23, // Water

    SO4 = 24, // SO4--

    PO4 = 25, // PO4--

    HET = 26, // ions

    UNKNOWN = 27

  };

  static constexpr auto num_residues() { return 27; }

  //string-to-enum translation
  constexpr residue parse_residue_name(std::string_view name) {
    if (name == "ACE") return residue::ACE;
    if (name == "ALA") return residue::ALA;
    if (name == "ARG") return residue::ARG;
    if (name == "ASN") return residue::ASN;
    if (name == "ASP") return residue::ASP;
    if (name == "CYS") return residue::CYS;
    if (name == "GLN") return residue::GLN;
    if (name == "GLU") return residue::GLU;
    if (name == "GLY") return residue::GLY;
    if (name == "HIS") return residue::HIS;
    if (name == "ILE") return residue::ILE;
    if (name == "LEU") return residue::LEU;
    if (name == "LYS") return residue::LYS;
    if (name == "MET") return residue::MET;
    if (name == "NME") return residue::NME;
    if (name == "PHE") return residue::PHE;
    if (name == "PRO") return residue::PRO;
    if (name == "SER") return residue::SER;
    if (name == "THR") return residue::THR;
    if (name == "TRP") return residue::TRP;
    if (name == "TYR") return residue::TYR;
    if (name == "VAL") return residue::VAL;
    if (name == "TER") return residue::TER;
    if (name == "HOH") return residue::HOH;
    if (name == "SO4") return residue::SO4;
    if (name == "PO4") return residue::PO4;
    if (name == "HET") return residue::HET;
    return residue::UNKNOWN;
  }

} // namespace mudock