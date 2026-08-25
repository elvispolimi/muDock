#pragma once

#include <array>
#include <cassert>
#include <mudock/chem/x_score_hb.hpp>
#include <mudock/chem/x_score_xtool_types.hpp>
#include <mudock/chem/residue.hpp>
#include <mudock/type_alias.hpp>
#include <span>
#include <string_view>
#include <unordered_map>

//===------------------------------------------------------------------------------------------------------
// WARNING: This file has been automatically generated from chem/x_score_residue_xtool_types.json and subsequently modified manually
//===------------------------------------------------------------------------------------------------------

namespace mudock {

  // reference residue enum for types

  // Description of an atom within a residue
  struct xtool_residue_atom_description {
    std::string_view name;
    xtool_ff basic_atom_type;
    xtool_ff x_tool_atom_type;
    fp_type vdw_radius;
    fp_type vdw_potential;
    fp_type par_charge;
    x_score_hb hbond;
    fp_type hydrophobic_scale;
    fp_type sas_parameter;
    int ring_indicator;
    std::string_view pmf_atom_type;
  };

  // Description of a bond within a residue
  struct xtool_residue_bond_description {
    std::string_view atom_1;
    std::string_view atom_2;
    int bond_type;
  };

  // Description of a residue
  struct xtool_residue_description {
    residue value;
    std::string_view name;
    fp_type total_charge;
    std::string_view description;
    std::span<const xtool_residue_atom_description> atoms;
    std::span<const xtool_residue_bond_description> bonds;
  };

  extern const std::array<xtool_residue_description, num_residues()> XTOOL_RESIDUE_DICTIONARY;
  extern const std::unordered_map<std::string_view, residue> XTOOL_RESIDUE_LOOKUP;

  // Utility function to get the description by enum
  inline const xtool_residue_description& get_description(const residue r) {
    // todo: check the asserts below
    assert(static_cast<int>(r) >= 0 && static_cast<int>(r) < num_residues());
    assert(XTOOL_RESIDUE_DICTIONARY[static_cast<int>(r)].value == r);
    return XTOOL_RESIDUE_DICTIONARY[static_cast<int>(r)];
  }

  // todo: function is probaby unnecessary once the residue is initialized as enum in the protein.
  // Utility function to get the residue by name (e.g., "ALA")
  inline const xtool_residue_description* get_residue_by_name(const std::string_view name) {
    if (auto it = XTOOL_RESIDUE_LOOKUP.find(name); it != XTOOL_RESIDUE_LOOKUP.end()) {
      return &get_description(it->second);
    }
    return nullptr;
  }

} // namespace mudock