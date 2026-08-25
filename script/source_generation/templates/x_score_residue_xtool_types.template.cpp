#include <mudock/chem/x_score_residue_xtool_types.hpp>

namespace mudock {

  // --- Static sub-arrays for atoms and bonds ---
{% for residue in residue_type %}
  static const xtool_residue_atom_description {@ residue.name @}_atoms[] = {
{% for atom in residue.atoms %}
    {
      "{@ atom.name @}",
      xtool_ff::{@ atom.basic_atom_type @},
      xtool_ff::{@ atom.x_tool_atom_type @},
      {@ atom.vdw_radius @}f,
      {@ atom.vdw_potential @}f,
      {@ atom.par_charge @}f,
      "{@ atom.hbond @}",
      {@ atom.hydrophobic_scale @}f,
      {@ atom.sas_parameter @}f,
      {@ atom.ring_indicator @},
      "{@ atom.pmf_atom_type @}"
    }{@ "," if not loop.last @}
{% endfor %}
  };

  static const xtool_residue_bond_description {@ residue.name @}_bonds[] = {
{% for bond in residue.bonds %}
    { "{@ bond.atom_1 @}", "{@ bond.atom_2 @}", {@ bond.bond_type @} }{@ "," if not loop.last @}
{% endfor %}
  };
{% endfor %}

  // --- Main Dictionary ---
  const std::array<xtool_residue_description, num_residues()> XTOOL_RESIDUE_DICTIONARY = {{
{% for residue in residue_type %}
    {
      xtool_residue::{@ residue.name @},
      "{@ residue.name @}",
      {@ residue.total_charge @}f,
      "{@ residue.description @}",
      {@ residue.name @}_atoms,
      {@ residue.name @}_bonds
    }{@ "," if not loop.last @}
{% endfor %}
  }};

  // --- Lookup Map ---
  const std::unordered_map<std::string_view, xtool_residue> XTOOL_RESIDUE_LOOKUP = {
{% for residue in residue_type %}
    { "{@ residue.name @}", xtool_residue::{@ residue.name @} }{@ "," if not loop.last @}
{% endfor %}
  };

} // namespace mudock
