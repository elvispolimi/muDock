#pragma once

#include <mudock/chem/autodock_molecule.hpp>
#include <mudock/chem/autodock_protein.hpp>
#include <mudock/molecule/fragments.hpp>

namespace mudock {
  struct autodock_ligand: public autodock_static_molecule {
    [[nodiscard]] inline auto non_bond_size() const { return non_bond_list_a1.size(); }
    [[nodiscard]] inline auto* non_bond_A() const { return non_bond_list_a1.data(); }
    [[nodiscard]] inline auto* non_bond_B() const { return non_bond_list_a2.data(); }
    [[nodiscard]] inline auto* non_bond_cA() const { return cA_v.data(); }
    [[nodiscard]] inline auto* non_bond_cB() const { return cB_v.data(); }
    [[nodiscard]] inline auto* non_bond_xB() const { return xB_v.data(); }

    [[nodiscard]] inline auto* fragments_masks() const { return frag_masks.data(); }
    [[nodiscard]] inline auto* fragmets_starts() const { return frag_start_indexes.data(); }
    [[nodiscard]] inline auto* fragments_stops() const { return frag_stop_indexes.data(); }

    [[nodiscard]] inline auto* atom_map_offsets() const { return map_offset_per_atom.data(); }
    [[nodiscard]] inline auto* atom_map_index() const { return map_index_per_atom.data(); }

    inline void update_offsets(const autodock_protein& adt_protein) {
      const auto atom_map_size = adt_protein.get_map_flat_size();
      for (int i = 0; i < (*this).num_atoms(); i++)
        map_offset_per_atom[i] = static_cast<int>(map_index_per_atom[i]) * atom_map_size;
    };

    // TODO missing remove atom etc...

    void prepare() {
      autodock_static_molecule::prepare();
      non_bond_list();
      precompute_lennard_jones();
      for (int i = 0; i < (*this).num_atoms(); i++)
        map_index_per_atom[i] = autodock_grid_from_ff((*this).autodock_type(i));
      get_linearized_fragments_mask((*this).num_atoms(),
                                    (*this).num_rotamers(),
                                    frag_masks,
                                    frag_start_indexes,
                                    frag_stop_indexes,
                                    (*this));
    };

  private:
    std::vector<int> non_bond_list_a1, non_bond_list_a2;
    std::vector<fp_type> cA_v, cB_v;
    std::vector<int> xB_v;
    std::vector<int> frag_masks;
    std::vector<int> frag_start_indexes;
    std::vector<int> frag_stop_indexes;
    static_containers::atoms_size<autodock_grid_type> map_index_per_atom;
    static_containers::atoms_size<int> map_offset_per_atom;

    void nonbonds(md_container<std::vector<uint_fast8_t>, 2>& nbmatrix,
                  const std::span<const bond> ligand_bond,
                  const int num_atoms);

    // weedbonds.cc for nonbondlist only for the first group
    /*___________________________________________________________________________
  |    ENDBRANCH---TORS---BRANCH---R O O T---BRANCH---ENDBRANCH                |
  |                                  /              \                          |
  |                                BRANCH            BRANCH--TORS---ENDBRANCH  |
  |                                /                  \                        |
  |                               ENDBRANCH            ENDBRANCH               |
  |____________________________________________________________________________|
  |  Eliminate all rigidly bonded atoms:                                       |
  |                                     * those atoms which are at the ROOT;   |
  |                                     * atoms between TORS and BRANCH;       |
  |                                     * atoms between BRANCH and ENDBRANCH.  |
  |  This is necessary for internal energy calculations.                       |
  |____________________________________________________________________________|
  | Weed out bonds in rigid pieces,                                            |
  |____________________________________________________________________________|
  */
    void weed_bonds(md_container<std::vector<uint_fast8_t>, 2>& nbmatrix,
                    const int num_atoms,
                    const fragments<static_containers>& ligand_fragments);

    void non_bond_list();

    void precompute_lennard_jones();
  };
} // namespace mudock
