

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <memory>
#include <unordered_map>
#include <mudock/format/mol2.hpp>
#include <mudock/format/reader.hpp>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MonomerInfo.h>
#include <stdexcept>
#include <string_view>

namespace mudock {
  namespace {
    struct pdbqt_atom_record {
      int serial{};
      std::string atom_name;
      std::string residue_name;
      int residue_id{1};
      fp_type x{}, y{}, z{}, charge{};
      autodock_ff type{};
    };

    std::string pdbqt_element(const autodock_ff type) {
      const auto token = get_description(type).name;
      if (token == "A" || token == "C") return "C";
      if (token == "NA" || token == "NS" || token == "N") return "N";
      if (token == "OA" || token == "OS" || token == "O") return "O";
      if (token == "SA" || token == "S") return "S";
      if (token == "HD" || token == "HS" || token == "H") return "H";
      return std::string(token);
    }

    template<class molecule_type>
      requires is_molecule<molecule_type>
    molecule_type convert_rdkit_molecule(const RDKit::ROMol& source) {
      if (source.getNumAtoms() == 0) throw std::runtime_error("Molecular record contains no atoms");
      molecule_type molecule;
      molecule.resize(static_cast<int>(source.getNumAtoms()), static_cast<int>(source.getNumBonds()));
      for (unsigned int i = 0; i < source.getNumAtoms(); ++i) {
        const auto* atom = source.getAtomWithIdx(i);
        const auto index = static_cast<int>(i);
        molecule.elements(index) = parse_element_symbol(atom->getSymbol());
        molecule.autodock_type(index) = autodock_ff::C;
        molecule.num_hbond(index) = 0;
        molecule.charge(index) = static_cast<fp_type>(atom->getFormalCharge());
        const auto& position = source.getConformer().getAtomPos(i);
        molecule.x(index) = static_cast<fp_type>(position.x);
        molecule.y(index) = static_cast<fp_type>(position.y);
        molecule.z(index) = static_cast<fp_type>(position.z);
        molecule.is_aromatic(index) = atom->getIsAromatic();
        if (const auto* residue = dynamic_cast<const RDKit::AtomPDBResidueInfo*>(atom->getMonomerInfo())) {
          molecule.atom_name(index) = residue->getName();
          molecule.residue_name(index) = residue->getResidueName();
          molecule.residue_id(index) = residue->getResidueNumber();
        }
      }
      int index = 0;
      for (const auto* edge : source.bonds()) {
        auto type = bond_type::SINGLE;
        switch (edge->getBondType()) {
          case RDKit::Bond::DOUBLE: type = bond_type::DOUBLE; break;
          case RDKit::Bond::TRIPLE: type = bond_type::TRIPLE; break;
          case RDKit::Bond::AROMATIC: type = bond_type::AROMATIC; break;
          default: break;
        }
        molecule.bonds(index++) = {static_cast<int>(edge->getBeginAtomIdx()),
                                   static_cast<int>(edge->getEndAtomIdx()),
                                   type,
                                   false};
      }
      molecule.set_prepared_types(false);
      return molecule;
    }

    template<class molecule_type>
      requires is_molecule<molecule_type>
    molecule_type parse_pdbqt_native(const std::string_view description) {
      molecule_type mol;
      std::vector<pdbqt_atom_record> atoms;
      std::vector<std::pair<int, int>> branches;
      std::istringstream input{std::string(description)};
      std::string line;
      std::string molecule_name;
      int torsdof = -1;
      auto field = [](const std::string& text, const std::size_t offset, const std::size_t count) {
        if (offset >= text.size()) return std::string{};
        return text.substr(offset, std::min(count, text.size() - offset));
      };
      while (std::getline(input, line)) {
        if (line.rfind("ATOM  ", 0) == 0 || line.rfind("HETATM", 0) == 0) {
          try {
            if (line.size() < 70) throw std::runtime_error("short PDBQT atom record");
            const auto type_token = detail::trim_copy(field(line, 77, 2));
            pdbqt_atom_record atom;
            atom.serial = std::stoi(field(line, 6, 5));
            atom.atom_name = detail::trim_copy(field(line, 12, 4));
            atom.residue_name = detail::trim_copy(field(line, 17, 3));
            atom.residue_id = std::stoi(field(line, 22, 4));
            atom.x = static_cast<fp_type>(std::stod(field(line, 30, 8)));
            atom.y = static_cast<fp_type>(std::stod(field(line, 38, 8)));
            atom.z = static_cast<fp_type>(std::stod(field(line, 46, 8)));
            atom.charge = static_cast<fp_type>(std::stod(field(line, 70, 6)));
            atom.type = parse_autodock_type(type_token);
            atoms.push_back(std::move(atom));
          } catch (const std::exception& e) {
            throw std::runtime_error("Invalid PDBQT atom record: " + std::string(e.what()));
          }
        } else if (line.rfind("BRANCH", 0) == 0) {
          std::istringstream record(line.substr(6));
          int first{}, second{};
          if (!(record >> first >> second)) throw std::runtime_error("Invalid PDBQT BRANCH record");
          branches.emplace_back(first, second);
        } else if (line.rfind("TORSDOF", 0) == 0) {
          std::istringstream record(line.substr(7));
          if (!(record >> torsdof) || torsdof < 0) throw std::runtime_error("Invalid PDBQT TORSDOF record");
        } else if (line.rfind(pdbqt::molecule_name_token, 0) == 0 && molecule_name.empty()) {
          molecule_name = detail::trim_copy(line.substr(pdbqt::molecule_name_token.size()));
        }
      }
      if (atoms.empty()) throw std::runtime_error("PDBQT contains no atom records");

      // RDKit is used only to infer connectivity from the coordinates. All
      // PDBQT fields and atom ordering below come from the source records.
      std::ostringstream pdb_block;
      for (std::size_t i = 0; i < atoms.size(); ++i) {
        const auto& atom = atoms[i];
        const auto element = pdbqt_element(atom.type);
        char atom_name[5]{};
        std::snprintf(atom_name, sizeof(atom_name), "%-4.4s", atom.atom_name.c_str());
        pdb_block << "HETATM" << std::setw(5) << atom.serial << ' ' << atom_name << ' '
                 << std::setw(3) << std::left << atom.residue_name.substr(0, 3) << std::right << " A"
                 << std::setw(4) << atom.residue_id << "    " << std::fixed << std::setprecision(3)
                 << std::setw(8) << atom.x << std::setw(8) << atom.y << std::setw(8) << atom.z
                 << "  1.00  0.00          " << std::setw(2) << std::right << element << "\n";
      }
      std::unique_ptr<RDKit::RWMol> rdkit_mol(RDKit::PDBBlockToMol(pdb_block.str(), false, false, 0, true));
      if (!rdkit_mol || rdkit_mol->getNumAtoms() != atoms.size())
        throw std::runtime_error("RDKit could not infer PDBQT atom connectivity");
      for (std::size_t i = 0; i < atoms.size(); ++i) {
        const auto& atom = atoms[i];
        const auto& rd_atom = *rdkit_mol->getAtomWithIdx(static_cast<unsigned int>(i));
        const auto& position = rdkit_mol->getConformer().getAtomPos(static_cast<unsigned int>(i));
        if (parse_element_symbol(rd_atom.getSymbol()) != parse_element_symbol(pdbqt_element(atom.type)) ||
            std::abs(position.x - static_cast<double>(atom.x)) > 0.001 ||
            std::abs(position.y - static_cast<double>(atom.y)) > 0.001 ||
            std::abs(position.z - static_cast<double>(atom.z)) > 0.001)
          throw std::runtime_error("RDKit changed PDBQT atom order or coordinates");
      }

      std::unordered_map<int, int> serial_to_index;
      for (std::size_t i = 0; i < atoms.size(); ++i) {
        if (!serial_to_index.emplace(atoms[i].serial, static_cast<int>(i)).second)
          throw std::runtime_error("PDBQT contains duplicate atom serial numbers");
      }
      std::vector<bond> bonds;
      bonds.reserve(rdkit_mol->getNumBonds());
      for (const auto* rd_bond : rdkit_mol->bonds()) {
        const auto a = static_cast<int>(rd_bond->getBeginAtomIdx());
        const auto b = static_cast<int>(rd_bond->getEndAtomIdx());
        bonds.push_back({a, b, bond_type::SINGLE, false});
      }
      for (const auto& [first, second] : branches) {
        const auto a = serial_to_index.find(first), b = serial_to_index.find(second);
        if (a == serial_to_index.end() || b == serial_to_index.end())
          throw std::runtime_error("PDBQT BRANCH references an unknown atom serial");
        const auto it = std::find_if(bonds.begin(), bonds.end(), [&](const bond& edge) {
          return (edge.source == a->second && edge.dest == b->second) ||
                 (edge.source == b->second && edge.dest == a->second);
        });
        if (it == bonds.end())
          bonds.push_back({a->second, b->second, bond_type::SINGLE, true});
        else
          it->can_rotate = true;
      }
      mol.resize(static_cast<int>(atoms.size()), static_cast<int>(bonds.size()));
      for (std::size_t i = 0; i < atoms.size(); ++i) {
        const auto& atom = atoms[i];
        mol.elements(static_cast<int>(i)) = parse_element_symbol(pdbqt_element(atom.type));
        mol.x(static_cast<int>(i)) = atom.x;
        mol.y(static_cast<int>(i)) = atom.y;
        mol.z(static_cast<int>(i)) = atom.z;
        mol.charge(static_cast<int>(i)) = atom.charge;
        mol.autodock_type(static_cast<int>(i)) = atom.type;
        const auto& ff = get_description(atom.type);
        mol.num_hbond(static_cast<int>(i)) = ff.hbond;
        mol.atom_name(static_cast<int>(i)) = atom.atom_name;
        mol.residue_name(static_cast<int>(i)) = atom.residue_name;
        mol.residue_id(static_cast<int>(i)) = atom.residue_id;
      }
      std::copy(bonds.begin(), bonds.end(), mol.get_bonds().begin());
      mol.set_torsdof(torsdof >= 0 ? torsdof : static_cast<int>(branches.size()));
      mol.set_prepared_types(true);
      if (!molecule_name.empty()) mol.properties.assign(property_type::NAME, std::move(molecule_name));
      return mol;
    }
  } // namespace

  template<>
  static_molecule parser<supported_format::PDBQT, static_molecule>(const std::string_view description) {
    return parse_pdbqt_native<static_molecule>(description);
  }
  template<>
  dynamic_molecule parser<supported_format::PDBQT, dynamic_molecule>(const std::string_view description) {
    return parse_pdbqt_native<dynamic_molecule>(description);
  }
  template<>
  static_molecule parser<supported_format::MOL2, static_molecule>(const std::string_view description) {
    static_molecule mol;
    mol2::parse(mol, description);
    return mol;
  }
  template<>
  dynamic_molecule parser<supported_format::MOL2, dynamic_molecule>(const std::string_view description) {
    dynamic_molecule mol;
    mol2::parse(mol, description);
    return mol;
  }
  template<>
  static_molecule parser<supported_format::PDB, static_molecule>(const std::string_view description) {
    std::unique_ptr<RDKit::RWMol> parsed(RDKit::PDBBlockToMol(std::string(description), false, false, 0, true));
    if (!parsed) throw std::runtime_error("RDKit failed to parse PDB block");
    return convert_rdkit_molecule<static_molecule>(*parsed);
  }
  template<>
  dynamic_molecule parser<supported_format::PDB, dynamic_molecule>(const std::string_view description) {
    std::unique_ptr<RDKit::RWMol> parsed(RDKit::PDBBlockToMol(std::string(description), false, false, 0, true));
    if (!parsed) throw std::runtime_error("RDKit failed to parse PDB block");
    return convert_rdkit_molecule<dynamic_molecule>(*parsed);
  }
  template<>
  static_molecule parser<supported_format::SDF, static_molecule>(const std::string_view description) {
    std::unique_ptr<RDKit::RWMol> parsed(RDKit::MolBlockToMol(std::string(description), false, false, true));
    if (!parsed) throw std::runtime_error("RDKit failed to parse SDF record");
    return convert_rdkit_molecule<static_molecule>(*parsed);
  }
  template<>
  dynamic_molecule parser<supported_format::SDF, dynamic_molecule>(const std::string_view description) {
    std::unique_ptr<RDKit::RWMol> parsed(RDKit::MolBlockToMol(std::string(description), false, false, true));
    if (!parsed) throw std::runtime_error("RDKit failed to parse SDF record");
    return convert_rdkit_molecule<dynamic_molecule>(*parsed);
  }
  template<>
  static_molecule parser<supported_format::MOL, static_molecule>(const std::string_view description) {
    std::unique_ptr<RDKit::RWMol> parsed(RDKit::MolBlockToMol(std::string(description), false, false, true));
    if (!parsed) throw std::runtime_error("RDKit failed to parse MOL block");
    return convert_rdkit_molecule<static_molecule>(*parsed);
  }
  template<>
  dynamic_molecule parser<supported_format::MOL, dynamic_molecule>(const std::string_view description) {
    std::unique_ptr<RDKit::RWMol> parsed(RDKit::MolBlockToMol(std::string(description), false, false, true));
    if (!parsed) throw std::runtime_error("RDKit failed to parse MOL block");
    return convert_rdkit_molecule<dynamic_molecule>(*parsed);
  }
} // namespace mudock
