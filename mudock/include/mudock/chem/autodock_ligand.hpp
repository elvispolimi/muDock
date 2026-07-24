#pragma once

#include <functional>
#include <mudock/chem/autodock_layer.hpp>
#include <mudock/chem/pdbqt_forcefield_override.hpp>
#include <mudock/molecule/fragments.hpp>

namespace mudock {
  struct autodock_ligand: public autodock_static_layer {
    autodock_ligand(static_molecule& _molecule, std::function<void(static_molecule&)> f = {})
        : autodock_static_layer(_molecule, f) {
      check_source_path_pdbqt_forcefield_override(*this, _molecule);
      prepare();
    };
    autodock_ligand(static_molecule& _molecule, std::function<void(autodock_static_layer&)> f)
        : autodock_static_layer(_molecule, f) {
      prepare();
    };
    [[nodiscard]] inline auto non_bond_size() const { return non_bond_list_a1.size(); }
    [[nodiscard]] inline auto* non_bond_A() const { return non_bond_list_a1.data(); }
    [[nodiscard]] inline auto* non_bond_B() const { return non_bond_list_a2.data(); }
    [[nodiscard]] inline auto* non_bond_cA() const { return cA_v.data(); }
    [[nodiscard]] inline auto* non_bond_cB() const { return cB_v.data(); }
    [[nodiscard]] inline auto* non_bond_xB() const { return xB_v.data(); }

    [[nodiscard]] inline auto* atom_map_offsets() const { return map_offset_per_atom.data(); }
    [[nodiscard]] inline auto* atom_map_index() const { return map_index_per_atom.data(); }

    inline void update_offsets(const int map_flat_size) {
      for (int i = 0; i < (*this).num_atoms(); i++)
        map_offset_per_atom[i] = static_cast<int>(map_index_per_atom[i]) * map_flat_size;
    };

    // TODO missing remove atom etc...

  private:
    std::vector<int> non_bond_list_a1, non_bond_list_a2;
    std::vector<fp_type> cA_v, cB_v;
    std::vector<int> xB_v;
    static_containers::atoms_size<autodock_grid_type> map_index_per_atom;
    static_containers::atoms_size<int> map_offset_per_atom;

    void prepare() {
      // autodock_static_layer::prepare(f);
      non_bond_list();
      precompute_lennard_jones();
      for (int i = 0; i < (*this).num_atoms(); i++)
        map_index_per_atom[i] = autodock_grid_from_ff((*this)().autodock_type(i));
    };

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
