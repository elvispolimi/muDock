#pragma once

#include "mudock/chem/autodock_types.hpp"

#include <mudock/format/ob_wrapper.hpp>
#include <mudock/molecule.hpp>
#include <string>
#include <string_view>

namespace mudock {
  //
  // class pdbqt {
  // public:
  //   static constexpr auto PDBQT_ATOM_TOKEN = "ATOM";
  //
  //   template<class molecule_type>
  //     requires is_molecule<molecule_type>
  //   void parse(molecule_type&& molecule, std::string_view molecule_description) const;
  //
  //   template<class molecule_type>
  //     requires is_molecule<molecule_type>
  //   void adjust_aromaticity(molecule_type&& molecule, std::string_view molecule_description) const;
  // };
  //
  // //===------------------------------------------------------------------------------------------------------
  // // Out-of-class method definitions
  // //===------------------------------------------------------------------------------------------------------
  //
  // template<class molecule_type>
  //   requires is_molecule<molecule_type>
  // void pdbqt::parse(molecule_type&& molecule, std::string_view molecule_description) const {
  //   // NOTE: we use OpenBabel to parse everything and deal with all the chemistry complexity
  //   convert(molecule, format_parser<supported_format::PDBQT>(molecule_description));
  //   adjust_aromaticity(molecule, molecule_description);
  // }
  //
  // template<class molecule_type>
  //   requires is_molecule<molecule_type>
  // void pdbqt::adjust_aromaticity(molecule_type&& molecule, std::string_view molecule_description) const {
  //   std::stringstream desc{std::string{molecule_description}};
  //
  //   std::string line;
  //   // Read the header (first few lines) for the grid information
  //   while (std::getline(desc, line)) {
  //     std::istringstream stream{line};
  //     std::string tokens;
  //     int id;
  //     fp_type charge;
  //     std::string ad_type_s;
  //     if (line.find(PDBQT_ATOM_TOKEN) != std::string::npos) {
  //       stream >> tokens >> id >> tokens >> tokens >> tokens >> tokens >> tokens >> tokens >> tokens >>
  //           tokens >> tokens >> charge >> ad_type_s;
  //       // assert(charge == molecule.charge(id - 1));
  //       const auto pos = ad_type_s.find("A");
  //       if (pos != std::string::npos && (pos != 0 || ad_type_s == "A")) {
  //         // molecule.is_aromatic(id - 1) = true;
  //         // Not requires since forcefields are applied later on
  //         //molecule.autodock_type(id - 1) = parse_autodock_type(ad_type_s);
  //       }
  //     }
  //   }
  // }
} // namespace mudock
