#pragma once

#include <cctype>
#include <openbabel/atom.h>
#include <openbabel/elements.h>
#include <openbabel/residue.h>
#include <stdexcept>
#include <string>

namespace mudock {

  inline std::string trim_copy(std::string value) {
    while (!value.empty() && std::isspace(static_cast<unsigned char>(value.front()))) {
      value.erase(value.begin());
    }

    while (!value.empty() && std::isspace(static_cast<unsigned char>(value.back()))) { value.pop_back(); }

    return value;
  }

  inline std::string get_ob_atom_name(OpenBabel::OBAtom* atom) {
    if (atom == nullptr) {
      return {};
    }

    if (auto* residue = atom->GetResidue()) {
      auto name = trim_copy(residue->GetAtomID(atom));

      if (!name.empty()) {
        return name;
      }
    }

    return std::string(OpenBabel::OBElements::GetSymbol(atom->GetAtomicNum())) +
           std::to_string(atom->GetIdx());
  }

  inline int get_ob_residue_id(OpenBabel::OBAtom* atom) {
    if (atom == nullptr) {
      return 0;
    }

    if (auto* residue = atom->GetResidue()) {
      return residue->GetNum();
    }

    return 0;
  }

  inline std::string get_ob_residue_name(OpenBabel::OBAtom* atom) {
    if (atom == nullptr) {
      return {};
    }

    if (auto* residue = atom->GetResidue()) {
      auto name = trim_copy(residue->GetName());

      if (!name.empty()) {
        return name;
      }
    }

    return {};
  }
} // namespace mudock
