#include <GraphMol/FileParsers/FileParsers.h>
#include <mudock/format/rdkit_wrapper.hpp>
#include <stdexcept>
#include <string_view>

namespace mudock {

  static void compute_and_sanitize(const rw_mol_wrapper& molecule) {
    // NOTE: use partial specialization to avoid errors on kekule ambiguities
    auto failed_ops = unsigned{0};
    RDKit::MolOps::sanitizeMol(*molecule,
                               failed_ops,
                               RDKit::MolOps::SanitizeFlags::SANITIZE_FINDRADICALS |
                                   RDKit::MolOps::SanitizeFlags::SANITIZE_SETAROMATICITY |
                                   RDKit::MolOps::SanitizeFlags::SANITIZE_SETCONJUGATION |
                                   RDKit::MolOps::SanitizeFlags::SANITIZE_SETHYBRIDIZATION |
                                   RDKit::MolOps::SanitizeFlags::SANITIZE_SYMMRINGS);
    if (failed_ops != RDKit::MolOps::SANITIZE_NONE) {
      throw std::runtime_error("RDKit sanitize failed!");
    }
  }

  rw_mol_wrapper parse_mol2(const std::string_view description) {
    static constexpr auto sanitaze  = false;
    static constexpr auto remove_hs = false;
    auto molecule =
        std::unique_ptr<RDKit::RWMol>{RDKit::Mol2BlockToMol(std::string{description}, sanitaze, remove_hs)};
    compute_and_sanitize(molecule);
    return molecule;
  }

  bond_type parse_rdkit_bond_type(const RDKit::Bond::BondType bond_type) {
    switch (bond_type) {
      case RDKit::Bond::BondType::SINGLE: return bond_type::SINGLE;
      case RDKit::Bond::BondType::DOUBLE: return bond_type::DOUBLE;
      case RDKit::Bond::BondType::TRIPLE: return bond_type::TRIPLE;
      case RDKit::Bond::BondType::AROMATIC: return bond_type::AROMATIC;
      default: throw std::runtime_error("Unsopported RDKit bond type");
    }
  }

} // namespace mudock
