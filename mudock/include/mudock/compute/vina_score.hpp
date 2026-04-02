#pragma once

#include <boost/range/size.hpp>
#include <cassert>
#include <cstddef>
#include <cstring>
#include <mudock/batch.hpp>
#include <mudock/chem/autodock_ligand.hpp>
#include <mudock/compute/vina_score_kernel.hpp>
#ifndef __CUDACC__
#include <mudock/compute/buffer_utils.hpp>
#include <mudock/compute/scoring.hpp>
#include <mudock/compute/scratchpad.hpp>
#endif
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {
  
  #define MAX_INTERACTING_PAIRS_IN_BATCH 10000

  template<typename molecule_type>
    size_t remove_hydrogens(molecule_type& molecule) {
      size_t out = 0;
        for(int atom = 0; atom < molecule.num_atoms();) {                   
                element type = molecule.elements(atom);                     
                if(type == element::H) {
                  molecule.remove_atom(atom);
                  out++;
                }
                else atom++;                                                
      }    
      return out;
  }

  std::pair<std::vector<int>, std::vector<int>> get_interactive_pairs(const static_molecule& ligand);

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
        l_is_hbond_acceptor(_scratch->get_queue()),
        l_is_hbond_donor(_scratch->get_queue()),
        l_is_hydrophobic(_scratch->get_queue()),
        l_vdw_radius(_scratch->get_queue()),
        l_interacting_pairs_first(_scratch->get_queue()),
        l_interacting_pairs_second(_scratch->get_queue()),
        l_num_interacting_pairs(_scratch->get_queue()),
        l_interacting_pairs_offset(_scratch->get_queue()),
        device_scratch(_device_scratch) {
          if (!(*device_scratch).template exists<buffer_data_type::PROT_HYDROPHOBICS>()) {
            info("Removing hydrogens from ligand (", protein.num_atoms(), ")...");
            remove_hydrogens(protein);           
            info("Protein size reduced to ", protein.num_atoms());

            int num_atoms_protein = protein.num_atoms();

            auto &num_atoms           = (*device_scratch).template get<buffer_data_type::NUM_ATOMS>();
            auto &protein_x           = (*device_scratch).template get<buffer_data_type::X_COORDS>();
            auto &protein_y           = (*device_scratch).template get<buffer_data_type::Y_COORDS>();
            auto &protein_z           = (*device_scratch).template get<buffer_data_type::Z_COORDS>();
            auto &p_is_hbond_acceptor = (*device_scratch).template get<buffer_data_type::PROT_H_ACCETORS>();
            auto &p_is_hbond_donor    = (*device_scratch).template get<buffer_data_type::PROT_H_DONORS>();
            auto &p_is_hydrophobic    = (*device_scratch).template get<buffer_data_type::PROT_HYDROPHOBICS>();
            auto &p_vdw_radius        = (*device_scratch).template get<buffer_data_type::PROT_VDW_RADS>();

            num_atoms.alloc(1);
            protein_x.alloc(num_atoms_protein);
            protein_y.alloc(num_atoms_protein);
            protein_z.alloc(num_atoms_protein);
            p_is_hbond_acceptor.alloc(num_atoms_protein);
            p_is_hbond_donor.alloc(num_atoms_protein);
            p_is_hydrophobic.alloc(num_atoms_protein);
            p_vdw_radius.alloc(num_atoms_protein);

            num_atoms()[0] = num_atoms_protein;

            std::memcpy(
                protein_x(), 
                protein.get_x().data(), 
                num_atoms_protein * sizeof(fp_type)
                );

            std::memcpy(
                protein_y(), 
                protein.get_y().data(), 
                num_atoms_protein * sizeof(fp_type)
                );

            std::memcpy(
                protein_z(), 
                protein.get_z().data(),
                num_atoms_protein * sizeof(fp_type)
                );

            std::memcpy(
                p_is_hbond_acceptor(), 
                protein.get_is_hbond_acceptor().data(), 
                num_atoms_protein * sizeof(int)
                );

            std::memcpy(
                p_is_hbond_donor(), 
                protein.get_is_hbond_donor().data(), 
                num_atoms_protein * sizeof(int)
                );

            std::memcpy(
                p_is_hydrophobic(), 
                protein.get_is_hydrophobic().data(), 
                num_atoms_protein * sizeof(int)
                );

            std::memcpy(
                p_vdw_radius(), 
                protein.get_vdw_radius().data(), 
                num_atoms_protein * sizeof(fp_type)
                );

            num_atoms.copy_host2device();
            protein_x.copy_host2device();
            protein_y.copy_host2device();
            protein_z.copy_host2device();
            p_is_hbond_acceptor.copy_host2device();
            p_is_hbond_donor.copy_host2device();
            p_is_hydrophobic.copy_host2device();
            p_vdw_radius.copy_host2device();
          }
        }

      void prepare(batch<static_molecule> &batch) {
        batch_ligands                = batch.num_ligands;
        batch_atoms                  = batch.batch_max_atoms;
        auto q                       = (*this->scratch).get_queue();
        const int tot_atoms_in_batch = batch_ligands * batch_atoms;

        load_num_rotamers(batch, this->scratch);
        load_num_atoms(batch, this->scratch);

        auto &score_b = (*this->scratch).template get<buffer_data_type::SCORES>();
        if (load_scratchs<queue_type>(batch, this->scratch)) {
          score_b.alloc(batch_ligands);
          score_b.set_valid();
        }

        l_is_hbond_acceptor.alloc(tot_atoms_in_batch);
        l_is_hbond_donor.alloc(tot_atoms_in_batch);
        l_is_hydrophobic.alloc(tot_atoms_in_batch);
        l_vdw_radius.alloc(tot_atoms_in_batch);
        l_interacting_pairs_first.alloc(MAX_INTERACTING_PAIRS_IN_BATCH);
        l_interacting_pairs_second.alloc(MAX_INTERACTING_PAIRS_IN_BATCH);
        l_num_interacting_pairs.alloc(batch_ligands);
        l_interacting_pairs_offset.alloc(batch_ligands);

        int offset_interacting_pairs = 0;
        for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
          auto &ligand = *batch.molecules[ligand_index];

          info("Removing hydrogens from ligand of size ", ligand.num_atoms());
          remove_hydrogens(ligand);           
          info("Ligand size reduced to", ligand.num_atoms());

          const int stride_atoms = ligand_index * batch_atoms; 
          const int num_atoms = ligand.num_atoms();
          std::memcpy(
              (void *) (l_is_hbond_acceptor() + stride_atoms), 
                ligand.get_is_hbond_acceptor().data(), 
                num_atoms * sizeof(int)
          );

          std::memcpy(
              (void *) (l_is_hbond_donor() + stride_atoms), 
                ligand.get_is_hbond_donor().data(), 
                num_atoms * sizeof(int)
          );
          
          std::memcpy(
              (void *) (l_is_hydrophobic() + stride_atoms), 
                ligand.get_is_hydrophobic().data(), 
                num_atoms * sizeof(int)
          );
          
          std::memcpy(
              (void *) (l_vdw_radius() + stride_atoms), 
                ligand.get_vdw_radius().data(), 
                num_atoms * sizeof(fp_type)
          );

          /// Parse interacting pairs of the ligand
          auto [ip_first, ip_second] = get_interactive_pairs(ligand);
          int num_interacting_pairs = ip_first.size();

          assert(offset_interacting_pairs + num_interacting_pairs <= MAX_INTERACTING_PAIRS_IN_BATCH && "Number of interacting pairs exceeded the limit\n");

          std::memcpy(
              (void *) (l_interacting_pairs_first() + offset_interacting_pairs), 
                ip_first.data(), 
                num_interacting_pairs * sizeof(int)
          );
 
          std::memcpy(
              (void *) (l_interacting_pairs_second() + offset_interacting_pairs), 
                ip_second.data(), 
                num_interacting_pairs * sizeof(int)
          );

          std::memcpy((void *) (l_num_interacting_pairs() + ligand_index), &num_interacting_pairs, sizeof(int));
          std::memcpy((void *) (l_interacting_pairs_offset() + ligand_index), &offset_interacting_pairs, sizeof(int));
          offset_interacting_pairs += num_interacting_pairs;
        }

        l_is_hbond_acceptor.copy_host2device();
        l_is_hbond_donor.copy_host2device();
        l_is_hydrophobic.copy_host2device();
        l_vdw_radius.copy_host2device();
        l_interacting_pairs_first.copy_host2device();
        l_interacting_pairs_second.copy_host2device();
        l_num_interacting_pairs.copy_host2device();
        l_interacting_pairs_offset.copy_host2device();

        // TODO ask Davide about this performance
        // Bind to the kernel function
        const auto scores_per_ligand = score_b.num_elements() / batch_ligands;
        assert((*this->scratch).template get<buffer_data_type::X_SCRATCH>().num_elements() !=
            (scores_per_ligand * batch_ligands) &&
            "Number of scores per ligands does not match the allocated coordinates space");

        const int *num_atoms_b = (*this->scratch).template get<buffer_data_type::NUM_ATOMS>().dev_pointer();
        const int *num_rotamers_b =
          (*this->scratch).template get<buffer_data_type::NUM_ROTAMERS>().dev_pointer();
        const fp_type *x_scratch_b = (*this->scratch).template get<buffer_data_type::X_SCRATCH>().dev_pointer();
        const fp_type *y_scratch_b = (*this->scratch).template get<buffer_data_type::Y_SCRATCH>().dev_pointer();
        const fp_type *z_scratch_b = (*this->scratch).template get<buffer_data_type::Z_SCRATCH>().dev_pointer();

        /// Pointers to protein data
        const int num_atoms_protein = (*device_scratch).template get<buffer_data_type::NUM_ATOMS>().host_pointer()[0];
        const fp_type *protein_x_p =  (*device_scratch).template get<buffer_data_type::X_COORDS>().dev_pointer();
        const fp_type *protein_y_p =  (*device_scratch).template get<buffer_data_type::Y_COORDS>().dev_pointer();
        const fp_type *protein_z_p =  (*device_scratch).template get<buffer_data_type::Z_COORDS>().dev_pointer();
        const int *p_is_hbond_acceptor_p =(*device_scratch).template get<buffer_data_type::PROT_H_ACCETORS>().dev_pointer();
        const int *p_is_hbond_donor_p =  (*device_scratch).template get<buffer_data_type::PROT_H_DONORS>().dev_pointer();
        const int *p_is_hydrophobic_p =  (*device_scratch).template get<buffer_data_type::PROT_HYDROPHOBICS>().dev_pointer();
        const fp_type *p_vdw_radius_p = (*device_scratch).template get<buffer_data_type::PROT_VDW_RADS>().dev_pointer();

        // Ligand batch data
        const int *l_is_hbond_acceptor_b = l_is_hbond_acceptor.dev_pointer();
        const int *l_is_hbond_donor_b = l_is_hbond_donor.dev_pointer();
        const int *l_is_hydrophobic_b = l_is_hydrophobic.dev_pointer();
        const fp_type *l_vdw_radius_b = l_vdw_radius.dev_pointer();
        const int *l_interacting_pairs_first_b = l_interacting_pairs_first.dev_pointer();
        const int *l_interacting_pairs_second_b = l_interacting_pairs_second.dev_pointer();
        const int *l_num_interacting_pairs_b = l_num_interacting_pairs.dev_pointer();
        const int *l_interacting_pairs_offset_b = l_interacting_pairs_offset.dev_pointer();

        fp_type *scores_b = score_b.dev_pointer();

        kernel = std::make_unique<vina_score_kernel<queue_type>>(
            scores_per_ligand,
            batch_ligands,
            batch_atoms,
            num_atoms_protein,
            protein_x_p,
            protein_y_p,
            protein_z_p,
            p_is_hbond_acceptor_p,
            p_is_hbond_donor_p,
            p_is_hydrophobic_p,
            p_vdw_radius_p,
            num_atoms_b, 
            x_scratch_b,
            y_scratch_b,
            z_scratch_b,
            l_is_hbond_acceptor_b,
            l_is_hbond_donor_b,
            l_is_hydrophobic_b,
            l_vdw_radius_b,
            num_rotamers_b,
            l_interacting_pairs_first_b, 
            l_interacting_pairs_second_b,
            l_num_interacting_pairs_b, 
            l_interacting_pairs_offset_b,
            scores_b,
            q
        );
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
  
      buffer_vector<int, queue_type> l_is_hbond_acceptor; 
      buffer_vector<int, queue_type> l_is_hbond_donor; 
      buffer_vector<int, queue_type> l_is_hydrophobic;
      buffer_vector<fp_type, queue_type> l_vdw_radius;
      buffer_vector<int, queue_type> l_interacting_pairs_first;
      buffer_vector<int, queue_type> l_interacting_pairs_second;
      buffer_vector<int, queue_type> l_num_interacting_pairs;
      buffer_vector<int, queue_type> l_interacting_pairs_offset;
      
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
