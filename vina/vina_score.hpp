#pragma once

#include <cstring>
#include <mudock/batch.hpp>
#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/chem/autodock_ligand.hpp>
#include <mudock/chem/autodock_protein.hpp>
#ifndef __CUDACC__
#include <mudock/compute/buffer_utils.hpp>
#include <mudock/compute/scoring.hpp>
#include <mudock/compute/scratchpad.hpp>
#endif
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {

  std::pair<std::vector<int>, std::vector<int>> get_interactive_pairs(const mudock::static_molecule& ligand){

    std::pair<std::vector<int>, std::vector<int>> out;

    const std::span<const mudock::bond>& bonds = ligand.get_bonds(); 
    const size_t num_atom = ligand.num_atoms();

    std::unordered_map<int, std::vector<int>> atoms_in_fragment = get_atoms_in_frag(bonds, num_atom);

    const auto num_rotamers = atoms_in_fragment.size();

    for (size_t rot1 = 0; rot1 < num_rotamers; ++rot1) {

      /// const auto* bitmask_rot1 = frag_masks + rot1 * num_atoms;
      std::vector<int> atoms_rot1 = atoms_in_fragment[rot1];
      for(size_t i = 0; i < atoms_rot1.size(); ++i){

        int atom1 = atoms_rot1[i];

        for(size_t rot2 = rot1 + 1; rot2 < num_rotamers; ++rot2){

          /// const auto* bitmask_rot2 = frag_masks + rot2 * num_atoms;
          std::vector<int> atoms_rot2 = atoms_in_fragment[rot2];

          for(size_t j = 0; j < atoms_rot2.size(); ++j){

            int atom2 = atoms_rot2[j];

            /// Search for the atom2 in the neighbors of atom1
            bool found = false;
            for(size_t nb = 0; nb < mudock::max_static_neighbors() && !found; nb++) {
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

    info("Total interacting pairs: ", out.first.size());
    return out;
  }


  template<typename queue_type>
    // requires std::derived_from<queue_type, queue>
    int get_vina_score_batch(const int, std::shared_ptr<queue_type>);

#ifndef __CUDACC__
  // TODO check that the object type and the kernel impl are the same
  template<typename queue_type>
    struct vina_score: public scoring<queue_type> {
      vina_score(std::shared_ptr<scratchpad<queue_type>> _scratch,
          std::shared_ptr<scratchpad<queue_type>> _device_scratch,
          dynamic_molecule &protein)
        : scoring<queue_type>(_scratch),
        device_scratch(_device_scratch) {
          //TODO: Initialize fields related to the protein
        }

      void prepare(batch<static_molecule> &batch) {
        batch_ligands                = batch.num_ligands;
        batch_atoms                  = batch.batch_max_atoms;
        const int batch_non_bonds    = batch_ligands * batch_atoms * batch_atoms;
        const int tot_atoms_in_batch = batch_ligands * batch_atoms;
        auto q                       = (*this->scratch).get_queue();

        load_num_rotamers(batch, this->scratch);
        load_num_atoms(batch, this->scratch);

        auto &score_b = (*this->scratch).template get<buffer_data_type::SCORES>();
        if (load_scratchs<queue_type>(batch, this->scratch)) {
          score_b.alloc(batch_ligands);
          score_b.set_valid();
        }

        //TODO: Initialize fields related to the ligand

      
        // TODO ask Davide about this performance
        // Bind to the kernel function
        const auto scores_per_ligand = score_b.num_elements() / batch_ligands;
        assert((*this->scratch).template get<buffer_data_type::X_SCRATCH>().num_elements() !=
            (scores_per_ligand * batch_ligands) &&
            "Number of scores per ligands does not match the allocated coordinates space");

        const int *num_atoms_b = (*this->scratch).template get<buffer_data_type::NUM_ATOMS>().dev_pointer();
        const int *num_rotamers_b =
          (*this->scratch).template get<buffer_data_type::NUM_ROTAMERS>().dev_pointer();
        const int *num_nonbonds_b  = num_nonbond.dev_pointer();
        const fp_type *x_scratch_b = (*this->scratch).template get<buffer_data_type::X_SCRATCH>().dev_pointer();
        const fp_type *y_scratch_b = (*this->scratch).template get<buffer_data_type::Y_SCRATCH>().dev_pointer();
        const fp_type *z_scratch_b = (*this->scratch).template get<buffer_data_type::Z_SCRATCH>().dev_pointer();

        kernel = std::make_unique<vina_score_kernel<queue_type>>();
      }

      void operator()() {
        assert(
            (((*this->scratch).template get<buffer_data_type::SCORES>().num_elements() % batch_ligands) == 0) &&
            "Number of scores is not a multiple of ligands in the batch");
        assert(kernel && "Kernel method not yet prepared");
        (*kernel)();
      }

      private:
      int batch_ligands;
      int batch_atoms;

      std::shared_ptr<scratchpad<queue_type>> device_scratch;
      std::unique_ptr<vina_score_kernel<queue_type>> kernel;

      void teardown_impl(batch<static_molecule> &batch) override {
        assert(batch.num_ligands == batch_ligands && "Scoring algorithm received different batch for teardown");

        //TODO this could be an issue if the scores per population would be equal to 1;
        if (batch_ligands ==
            static_cast<int>((*this->scratch).template get<buffer_data_type::SCORES>().num_elements())) {
          auto &scores_b = (*this->scratch).template get<buffer_data_type::SCORES>();
          scores_b.copy_device2host();
          (*this->scratch).get_queue()->synchronize();
          for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
            auto &ligand = *batch.molecules[ligand_index];

            ligand.properties.assign(property_type::SCORE, std::to_string(scores_b()[ligand_index]));
          }
        }
        // Otherwise the upper stage gave the responsibility to do so
      };
    };
#endif
}
