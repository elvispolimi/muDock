#include <memory>
#include <mudock/format/ob_wrapper.hpp>
#include <stdexcept>
#include <string_view>

namespace mudock {

  ob_mol_wrapper parse_pdbqt(const std::string_view description) {
    OpenBabel::OBConversion conv;
    conv.SetInFormat("pdbqt");

    auto mol = std::unique_ptr<OpenBabel::OBMol>{};
    std::istringstream desc{std::string(description)};
    if (!conv.Read(mol.get(), &desc)) {
      mudock::error("Couldn't open PDBQT file");
      throw std::runtime_error("PDBQT Parser failed, look to logs for details");
    }
    return mol;
  }

  bond_type parse_ob_bond_type(const OpenBabel::OBBond& bond) {
    if (bond.IsAromatic())
      return bond_type::AROMATIC;
    else if (const_cast<OpenBabel::OBBond&>(bond).IsAmide())
      return bond_type::AMIDE;
    else
      switch (bond.GetBondOrder()) {
        case 1: return bond_type::SINGLE;
        case 2: return bond_type::DOUBLE;
        case 3: return bond_type::TRIPLE;
        default: throw std::runtime_error("Unsopported OpenBabel bond type");
      }
  }

} // namespace mudock
