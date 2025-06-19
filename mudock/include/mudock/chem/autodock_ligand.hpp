#pragma once

#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/chem/autodock_protein.hpp>
#include <mudock/cpp_implementation/weed_bonds.hpp>
#include <mudock/molecule.hpp>
#include <mudock/molecule/containers.hpp>
#include <mudock/type_alias.hpp>
#include <span>
#include <vector>

namespace mudock {
  struct autodock_ligand {
    autodock_ligand(static_molecule& _ligand): ligand(_ligand) {
      non_bond_list(ligand, non_bond_list_a1, non_bond_list_a2);
      precompute_lennard_jones(non_bond_list_a1.size(),
                               cA_v,
                               cB_v,
                               xB_v,
                               ligand,
                               non_bond_list_a1,
                               non_bond_list_a2);
      for (int i = 0; i < ligand.num_atoms(); i++)
        map_index_per_atom[i] = autodock_grid_from_ff(ligand.autodock_type(i));
      get_linearized_fragments_mask(ligand.num_atoms(),
                                    ligand.num_rotamers(),
                                    frag_masks,
                                    frag_start_indexes,
                                    frag_stop_indexes,
                                    ligand);
    };

    [[nodiscard]] inline auto get_num_atoms() const { return ligand.num_atoms(); }
    [[nodiscard]] inline auto get_num_rotatable_bonds() const { return ligand.num_rotamers(); }
    [[nodiscard]] inline auto& get_ligand() const { return ligand; }
    [[nodiscard]] inline auto get_non_bond_size() const { return non_bond_list_a1.size(); }
    [[nodiscard]] inline auto* get_non_bond_A() const { return non_bond_list_a1.data(); }
    [[nodiscard]] inline auto* get_non_bond_B() const { return non_bond_list_a2.data(); }
    [[nodiscard]] inline auto* get_non_bond_cA() const { return cA_v.data(); }
    [[nodiscard]] inline auto* get_non_bond_cB() const { return cB_v.data(); }
    [[nodiscard]] inline auto* get_non_bond_xB() const { return xB_v.data(); }

    [[nodiscard]] inline auto* get_ligand_x_p() { return ligand.get_x().data(); }
    [[nodiscard]] inline auto* get_ligand_y_p() { return ligand.get_y().data(); }
    [[nodiscard]] inline auto* get_ligand_z_p() { return ligand.get_z().data(); }
    [[nodiscard]] inline auto get_ligand_x() const { return ligand.get_x(); }
    [[nodiscard]] inline auto get_ligand_y() const { return ligand.get_y(); }
    [[nodiscard]] inline auto get_ligand_z() const { return ligand.get_z(); }
    [[nodiscard]] inline auto* get_ligand_vol() const { return ligand.get_vol().data(); }
    [[nodiscard]] inline auto* get_ligand_solpar() const { return ligand.get_solpar().data(); }
    [[nodiscard]] inline auto* get_ligand_charge() const { return ligand.get_charge().data(); }
    [[nodiscard]] inline auto* get_fragments_masks() const { return frag_masks.data(); }
    [[nodiscard]] inline auto* get_fragmets_starts() const { return frag_start_indexes.data(); }
    [[nodiscard]] inline auto* get_fragments_stops() const { return frag_stop_indexes.data(); }

    [[nodiscard]] inline auto* get_atom_map_offsets() const { return map_offset_per_atom.data(); }
    [[nodiscard]] inline auto* get_atom_map_index() const { return map_index_per_atom.data(); }

    inline void update_offsets(const autodock_protein& adt_protein) {
      const auto atom_map_size = adt_protein.get_map_flat_size();
      for (int i = 0; i < ligand.num_atoms(); i++)
        map_offset_per_atom[i] = static_cast<int>(map_index_per_atom[i]) * atom_map_size;
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

    static_molecule& ligand;
  };
} // namespace mudock
