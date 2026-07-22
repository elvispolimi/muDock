#include <mudock/compute/vina_score.hpp>

using namespace std;

namespace mudock {

  unordered_map<int, vector<int>> get_atoms_in_frag(
      const span<const bond>& bonds, 
      const size_t num_atom
      ){
    auto graph = make_graph(bonds, num_atom);
    const auto ligand_fragments =
      make_unique<fragments<static_containers>>(graph,
          bonds,
          num_atom);

    auto rigid_pieces = ligand_fragments.get()->get_rigid_pieces();

    unordered_map<int, vector<int>> atoms_in_fragment;

    for (size_t i = 0; i < num_atom; ++i) {
      atoms_in_fragment[rigid_pieces[i]].emplace_back(i);
    }

    return atoms_in_fragment;
  }


  pair<vector<int>, vector<int>> get_interactive_pairs(const static_molecule& ligand){

    pair<vector<int>, vector<int>> out;

    const span<const bond>& bonds = ligand.get_bonds(); 
    const size_t num_atom = ligand.num_atoms();

    unordered_map<int, vector<int>> atoms_in_fragment = get_atoms_in_frag(bonds, num_atom);

    const auto num_rotamers = atoms_in_fragment.size();

    for (size_t rot1 = 0; rot1 < num_rotamers; ++rot1) {

      /// const auto* bitmask_rot1 = frag_masks + rot1 * num_atoms;
      vector<int> atoms_rot1 = atoms_in_fragment[rot1];
      for(size_t i = 0; i < atoms_rot1.size(); ++i){

        int atom1 = atoms_rot1[i];

        for(size_t rot2 = rot1 + 1; rot2 < num_rotamers; ++rot2){

          /// const auto* bitmask_rot2 = frag_masks + rot2 * num_atoms;
          std::vector<int> atoms_rot2 = atoms_in_fragment[rot2];

          for(size_t j = 0; j < atoms_rot2.size(); ++j){

            int atom2 = atoms_rot2[j];

            /// Search for the atom2 in the neighbors of atom1
            bool found = false;
            for(size_t nb = 0; nb < max_static_neighbors() && !found; nb++) {
              int atom1_nb = ligand.neighbors(atom1, nb);
              if (atom1_nb == atom2) found = true;
              else if (atom1_nb == -1) break;
            }
            if (found) continue;

            /// Opendock: if [i, j] in self.torsion_bond_index or [j, i] in self.torsion_bond_index: continue

            //int mask_1 = bitmask_rot1[atom1];
            //int mask_2 = bitmask_rot2[atom2];
            ///if(mask_1 == 3 || mask_2 == 3 || mask_1 == 2 || mask_2 == 2) continue;

            /// Order pair before adding it to the output
            int a = std::min(atom1, atom2);
            int b = std::max(atom1, atom2);
            out.first.emplace_back(a);
            out.second.emplace_back(b);
          }
        }
      }
    }

    // info("Total interacting pairs: ", out.first.size());
    return out;
  }

}
