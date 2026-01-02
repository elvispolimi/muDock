#pragma once

#include <cassert>
#include <cmath>
#include <mudock/format/ob_wrapper.hpp>
#include <mudock/format/supported_format.hpp>
#include <mudock/log.hpp>
#include <openbabel/bond.h>
#include <openbabel/data.h>
#include <openbabel/elements.h>
#include <openbabel/generic.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/oberror.h>
#include <sys/types.h>

namespace mudock {

  template<supported_format format, class molecule_type>
  void writer(const molecule_type& mol, std::ofstream& ofs);

  template<>
  void writer<supported_format::PDBQT>(const ob_mol_wrapper& mol, std::ofstream& ofs);
  template<>
  void writer<supported_format::MOL2>(const ob_mol_wrapper& mol, std::ofstream& ofs);
  template<>
  void writer<supported_format::ADTMOL2>(const ob_mol_wrapper& mol, std::ofstream& ofs);
  template<>
  void writer<supported_format::PDB>(const ob_mol_wrapper& mol, std::ofstream& ofs);
  template<>
  void writer<supported_format::PDBQT>(const static_molecule& mol, std::ofstream& ofs);
  template<>
  void writer<supported_format::MOL2>(const static_molecule& mol, std::ofstream& ofs);
  template<>
  void writer<supported_format::ADTMOL2>(const static_molecule& mol, std::ofstream& ofs);
  template<>
  void writer<supported_format::PDB>(const static_molecule& mol, std::ofstream& ofs);
  template<>
  void writer<supported_format::PDBQT>(const dynamic_molecule& mol, std::ofstream& ofs);
  template<>
  void writer<supported_format::MOL2>(const dynamic_molecule& mol, std::ofstream& ofs);
  template<>
  void writer<supported_format::ADTMOL2>(const dynamic_molecule& mol, std::ofstream& ofs);
  template<>
  void writer<supported_format::PDB>(const dynamic_molecule& mol, std::ofstream& ofs);

  template<class molecule_type>
    requires is_molecule<molecule_type>
  void writer(const molecule_type& mol, const std::filesystem::path out_path) {
    const auto out_format = parse_supported_format(out_path);

    std::ofstream ofs(out_path, std::ios::out);
    constexpr_switch<0, get_num_supported_format(), 1>(
        [&](const auto format_index) {
          const auto format = static_cast<supported_format>(format_index());
          writer<format, molecule_type>(mol, ofs);
          return;
        },
        out_format);

    mudock::error("The provided path " + out_path.string() + " extension is not yet supported");
  }

} // namespace mudock
