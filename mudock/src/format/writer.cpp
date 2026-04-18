#include "mudock/chem/assign_autodock_types.hpp"

#include <cassert>
#include <fstream>
#include <memory>
#include <mudock/format/adt_mol2.hpp>
#include <mudock/format/mol2.hpp>
#include <mudock/format/ob_wrapper.hpp>
#include <mudock/format/supported_format.hpp>
#include <mudock/format/writer.hpp>
#include <mudock/molecule.hpp>
#include <mudock/utils.hpp>
#include <openbabel/atom.h>
#include <openbabel/babelconfig.h>
#include <openbabel/chargemodel.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>
#include <openbabel/obutil.h>
#include <openbabel/plugin.h>
#include <stdexcept>

namespace mudock {
  template<supported_format format>
  void ob_writer(const ob_mol_wrapper& mol, std::ofstream& ofs) {
    OpenBabel::OBConversion conv;
    const std::string ext{parse_supported_format(format)};
    conv.SetOutFormat(ext.c_str());

    if (!ofs) {
      throw std::runtime_error("Error: Cannot open file for writing!");
    }

    if (!conv.Write(mol.get(), &ofs)) {
      throw std::runtime_error("Error: Failed to write molecule!");
    }
  }

  // PDBQT
  template<>
  void writer<supported_format::PDBQT>(const ob_mol_wrapper& mol, std::ofstream& ofs) {
    // Add charges using the Gasteiger method
    OpenBabel::OBChargeModel* chargeModel = OpenBabel::OBChargeModel::FindType("gasteiger");
    if (!chargeModel) {
      mudock::error("Error: Unable to find charge model.");
      throw std::runtime_error("Error in OpenBabel transformations");
    }
    chargeModel->ComputeCharges(*mol.get());

    mol.get()->ConnectTheDots();
    mol.get()->PerceiveBondOrders();

    ob_writer<supported_format::PDBQT>(mol, ofs);
  };
  template<class molecule_type>
    requires is_molecule<molecule_type>
  void writer_impl_pdbqt(const molecule_type& mol, std::ofstream& ofs) {
    ob_mol_wrapper ob_mol = std::make_unique<OpenBabel::OBMol>();
    convert(ob_mol, mol);
    writer<supported_format::PDBQT>(ob_mol, ofs);
  }
  template<>
  void writer<supported_format::PDBQT>(const static_molecule& mol, std::ofstream& ofs) {
    writer_impl_pdbqt(mol, ofs);
  };
  template<>
  void writer<supported_format::PDBQT>(const dynamic_molecule& mol, std::ofstream& ofs) {
    writer_impl_pdbqt(mol, ofs);
  }

  // Helper function to convert to ob mol from molecule
  template<supported_format format, class molecule_type>
    requires is_molecule<molecule_type>
  void writer_impl(const molecule_type& mol, std::ofstream& ofs) {
    ob_mol_wrapper ob_mol = std::make_unique<OpenBabel::OBMol>();
    convert(ob_mol, mol);
    ob_writer<format>(ob_mol, ofs);
  }

  // MOL2
  template<>
  void writer<supported_format::MOL2>(const ob_mol_wrapper& mol, std::ofstream& ofs) {
    dynamic_molecule s_mol;
    convert(s_mol, mol);
    writer<supported_format::MOL2>(s_mol, ofs);
  };
  template<>
  void writer<supported_format::MOL2>(const static_molecule& mol, std::ofstream& ofs) {
    mol2::print(mol, ofs);
  }
  template<>
  void writer<supported_format::MOL2>(const dynamic_molecule& mol, std::ofstream& ofs) {
    mol2::print(mol, ofs);
  }

  // PDB
  template<>
  void writer<supported_format::PDB>(const ob_mol_wrapper& mol, std::ofstream& ofs) {
    ob_writer<supported_format::PDB>(mol, ofs);
  }
  template<>
  void writer<supported_format::PDB>(const static_molecule& mol, std::ofstream& ofs) {
    writer_impl<supported_format::PDB>(mol, ofs);
  }
  template<>
  void writer<supported_format::PDB>(const dynamic_molecule& mol, std::ofstream& ofs) {
    writer_impl<supported_format::PDB>(mol, ofs);
  }

  template<>
  void writer<supported_format::ADTMOL2>(const ob_mol_wrapper& mol, std::ofstream& ofs) {
    dynamic_molecule s_mol;
    convert(s_mol, mol);
    writer<supported_format::ADTMOL2>(s_mol, ofs);
  };
  template<>
  void writer<supported_format::ADTMOL2>(const static_molecule& mol, std::ofstream& ofs) {
    auto temp = mol;
    assign_autodock_types(temp);
    adt_mol2::print(temp, ofs);
  }
  template<>
  void writer<supported_format::ADTMOL2>(const dynamic_molecule& mol, std::ofstream& ofs) {
    auto temp = mol;
    assign_autodock_types(temp);
    adt_mol2::print(temp, ofs);
  }
} // namespace mudock
