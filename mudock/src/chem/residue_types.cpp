#include <algorithm>
#include <cctype>
#include <mudock/chem/residue_types.hpp>
#include <stdexcept>
#include <string>

namespace mudock {

  namespace {

    std::string upper_copy(std::string_view input) {
      std::string output(input);
      std::transform(output.begin(), output.end(), output.begin(), [](unsigned char c) {
        return static_cast<char>(std::toupper(c));
      });
      return output;
    }

    std::string trim_copy(std::string_view input) {
      const auto begin = input.find_first_not_of(" \t\r\n");
      if (begin == std::string_view::npos) {
        return {};
      }
      const auto end = input.find_last_not_of(" \t\r\n");
      return std::string(input.substr(begin, end - begin + 1));
    }

  } // namespace

  const std::array<residue_type_description, num_residue_types()> RESIDUE_TYPE_DICTIONARY = {{
      {residue_type::UNKNOWN, "", false, false, false, false, false},

      {residue_type::ALA, "ALA", true, false, false, false, false},
      {residue_type::ARG, "ARG", true, false, false, false, false},
      {residue_type::ASN, "ASN", true, false, false, false, false},
      {residue_type::ASP, "ASP", true, false, false, false, false},
      {residue_type::CYS, "CYS", true, false, false, false, false},
      {residue_type::GLN, "GLN", true, false, false, false, false},
      {residue_type::GLU, "GLU", true, false, false, false, false},
      {residue_type::GLY, "GLY", true, false, false, false, false},
      {residue_type::HIS, "HIS", true, false, false, false, false},
      {residue_type::HSD, "HSD", true, false, false, false, false},
      {residue_type::HSE, "HSE", true, false, false, false, false},
      {residue_type::HSP, "HSP", true, false, false, false, false},
      {residue_type::HID, "HID", true, false, false, false, false},
      {residue_type::HIE, "HIE", true, false, false, false, false},
      {residue_type::HIP, "HIP", true, false, false, false, false},
      {residue_type::ILE, "ILE", true, false, false, false, false},
      {residue_type::LEU, "LEU", true, false, false, false, false},
      {residue_type::LYS, "LYS", true, false, false, false, false},
      {residue_type::MET, "MET", true, false, false, false, false},
      {residue_type::PHE, "PHE", true, false, false, false, false},
      {residue_type::PRO, "PRO", true, false, false, false, false},
      {residue_type::SER, "SER", true, false, false, false, false},
      {residue_type::THR, "THR", true, false, false, false, false},
      {residue_type::TRP, "TRP", true, false, false, false, false},
      {residue_type::TYR, "TYR", true, false, false, false, false},
      {residue_type::VAL, "VAL", true, false, false, false, false},

      {residue_type::ASH, "ASH", true, false, false, false, false},
      {residue_type::CYM, "CYM", true, false, false, false, false},
      {residue_type::CYX, "CYX", true, false, false, false, false},
      {residue_type::GLH, "GLH", true, false, false, false, false},
      {residue_type::LYN, "LYN", true, false, false, false, false},
      {residue_type::MSE, "MSE", true, false, false, false, false},
      {residue_type::NME, "NME", false, false, false, false, true},
      {residue_type::ACE, "ACE", false, false, false, false, true},

      {residue_type::A, "A", false, true, false, false, false},
      {residue_type::C, "C", false, true, false, false, false},
      {residue_type::G, "G", false, true, false, false, false},
      {residue_type::T, "T", false, true, false, false, false},
      {residue_type::U, "U", false, true, false, false, false},
      {residue_type::DA, "DA", false, true, false, false, false},
      {residue_type::DC, "DC", false, true, false, false, false},
      {residue_type::DG, "DG", false, true, false, false, false},
      {residue_type::DT, "DT", false, true, false, false, false},
      {residue_type::DU, "DU", false, true, false, false, false},
      {residue_type::RA, "RA", false, true, false, false, false},
      {residue_type::RC, "RC", false, true, false, false, false},
      {residue_type::RG, "RG", false, true, false, false, false},
      {residue_type::RT, "RT", false, true, false, false, false},
      {residue_type::RU, "RU", false, true, false, false, false},
      {residue_type::ADE, "ADE", false, true, false, false, false},
      {residue_type::CYT, "CYT", false, true, false, false, false},
      {residue_type::GUA, "GUA", false, true, false, false, false},
      {residue_type::THY, "THY", false, true, false, false, false},
      {residue_type::URA, "URA", false, true, false, false, false},

      {residue_type::HOH, "HOH", false, false, true, false, false},
      {residue_type::WAT, "WAT", false, false, true, false, false},
      {residue_type::SOL, "SOL", false, false, true, false, false},
      {residue_type::TIP3, "TIP3", false, false, true, false, false},
      {residue_type::SPC, "SPC", false, false, true, false, false},
      {residue_type::NA, "NA", false, false, false, true, false},
      {residue_type::K, "K", false, false, false, true, false},
      {residue_type::CL, "CL", false, false, false, true, false},
      {residue_type::CA, "CA", false, false, false, true, false},
      {residue_type::MG, "MG", false, false, false, true, false},
      {residue_type::MN, "MN", false, false, false, true, false},
      {residue_type::FE, "FE", false, false, false, true, false},
      {residue_type::CO, "CO", false, false, false, true, false},
      {residue_type::NI, "NI", false, false, false, true, false},
      {residue_type::CU, "CU", false, false, false, true, false},
      {residue_type::ZN, "ZN", false, false, false, true, false},
      {residue_type::CD, "CD", false, false, false, true, false},
      {residue_type::HG, "HG", false, false, false, true, false},

      {residue_type::LIG, "LIG", false, false, false, false, true},
      {residue_type::UNL, "UNL", false, false, false, false, true},
      {residue_type::UNK, "UNK", false, false, false, false, true},
      {residue_type::MOL, "MOL", false, false, false, false, true},
      {residue_type::DRG, "DRG", false, false, false, false, true},
      {residue_type::INH, "INH", false, false, false, false, true},

      {residue_type::ATP, "ATP", false, false, false, false, true},
      {residue_type::ADP, "ADP", false, false, false, false, true},
      {residue_type::AMP, "AMP", false, false, false, false, true},
      {residue_type::GTP, "GTP", false, false, false, false, true},
      {residue_type::GDP, "GDP", false, false, false, false, true},
      {residue_type::GMP, "GMP", false, false, false, false, true},
      {residue_type::NAD, "NAD", false, false, false, false, true},
      {residue_type::NAP, "NAP", false, false, false, false, true},
      {residue_type::NDP, "NDP", false, false, false, false, true},
      {residue_type::FAD, "FAD", false, false, false, false, true},
      {residue_type::FMN, "FMN", false, false, false, false, true},
      {residue_type::HEM, "HEM", false, false, false, false, true},
      {residue_type::HEC, "HEC", false, false, false, false, true},
      {residue_type::HEB, "HEB", false, false, false, false, true},
      {residue_type::NAG, "NAG", false, false, false, false, true},
      {residue_type::NDG, "NDG", false, false, false, false, true},
      {residue_type::BMA, "BMA", false, false, false, false, true},
      {residue_type::MAN, "MAN", false, false, false, false, true},
      {residue_type::GAL, "GAL", false, false, false, false, true},
      {residue_type::GLC, "GLC", false, false, false, false, true},
      {residue_type::SO4, "SO4", false, false, false, false, true},
      {residue_type::PO4, "PO4", false, false, false, false, true},
      {residue_type::ACT, "ACT", false, false, false, false, true},
      {residue_type::EDO, "EDO", false, false, false, false, true},
      {residue_type::GOL, "GOL", false, false, false, false, true},
      {residue_type::DMS, "DMS", false, false, false, false, true},
      {residue_type::DMSO, "DMSO", false, false, false, false, true},
  }};

  residue_type parse_residue_type(std::string_view token) {
    const auto normalized = upper_copy(trim_copy(token));
    if (normalized.empty()) {
      return residue_type::UNKNOWN;
    }

    for (const auto& entry: RESIDUE_TYPE_DICTIONARY) {
      if (entry.name == normalized) {
        return entry.value;
      }
    }

    return residue_type::UNKNOWN;
  }

  residue_type parse_residue_type_or_throw(std::string_view token) {
    const auto type = parse_residue_type(token);
    if (type == residue_type::UNKNOWN) {
      throw std::runtime_error("Unsupported residue/substructure type: " + std::string(token));
    }
    return type;
  }

} // namespace mudock
