#pragma once

#include <cstddef>
#include <cstring>
#include <mudock/batch.hpp>
#include <mudock/chem/vinardo_ligand.hpp>
#include <mudock/chem/vinardo_protein.hpp>
#include <mudock/compute/batch_multiple.hpp>
#include <mudock/compute/vinardo_score_kernel.hpp>
#if !defined(__CUDACC__) && !defined(__HIPCC__)
  #include <mudock/compute/buffer_utils.hpp>
  #include <mudock/compute/scoring.hpp>
  #include <mudock/compute/scratchpad.hpp>
#endif
#include <mudock/log.hpp>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {

#if !defined(__CUDACC__) && !defined(__HIPCC__)
  // TODO check that the object type and the kernel impl are the same
  template<typename queue_type>
  struct vinardo_score: public scoring<queue_type> {
    static constexpr const char stage_name[] = "VINARDO";
    //When I create buffers I pass them the queue that knows how to manage memory
    //then for all operations I talk directly with the buffer
    //in the initialization list I have to prepare all buffers that I will need for the ligands
    vinardo_score(std::shared_ptr<scratchpad<queue_type>> _scratch,        //Here I should insert the data that changes between batch
                  std::shared_ptr<scratchpad<queue_type>> _device_scratch, //This should be the container for the stable data
                  dynamic_molecule& protein)
        : scoring<queue_type>(_scratch),
          pl_offsets(_scratch->get_queue()),
          pl_counts(_scratch->get_queue()),
          pl_protein_atom_idx(_scratch->get_queue()),
          pl_ligand_atom_idx(_scratch->get_queue()),
          pl_radius_sum(_scratch->get_queue()),
          pl_hydrophobic_possible(_scratch->get_queue()),
          pl_hbond_possible(_scratch->get_queue()),
          ll_offsets(_scratch->get_queue()),
          ll_counts(_scratch->get_queue()),
          ll_atom_i_idx(_scratch->get_queue()),
          ll_atom_j_idx(_scratch->get_queue()),
          ll_radius_sum(_scratch->get_queue()),
          ll_hydrophobic_possible(_scratch->get_queue()),
          ll_hbond_possible(_scratch->get_queue()),
          vinardo_num_tors(_scratch->get_queue()),
          inter_scores(_scratch->get_queue()),
          intra_scores(_scratch->get_queue()),
          protein_vinardo(protein),
          device_scratch(_device_scratch) {

      //If protein data are not yet on the device we must put them there
      if (!device_scratch->template exists<buffer_data_type::PROT_X>() ||
          !device_scratch->template exists<buffer_data_type::PROT_Y>() ||
          !device_scratch->template exists<buffer_data_type::PROT_Z>()) {

        //retrieve the buffers
        auto& prot_x = device_scratch->template get<buffer_data_type::PROT_X>();
        auto& prot_y = device_scratch->template get<buffer_data_type::PROT_Y>();
        auto& prot_z = device_scratch->template get<buffer_data_type::PROT_Z>();
        //Allocate memory
        const auto protein_atoms = protein.num_atoms();
        prot_x.alloc(protein_atoms);
        prot_y.alloc(protein_atoms);
        prot_z.alloc(protein_atoms);
        //Fill them
        std::memcpy(prot_x(), protein.x(), protein_atoms * sizeof(fp_type));
        std::memcpy(prot_y(), protein.y(), protein_atoms * sizeof(fp_type));
        std::memcpy(prot_z(), protein.z(), protein_atoms * sizeof(fp_type));
        //Move data to device
        prot_x.copy_host2device();
        prot_y.copy_host2device();
        prot_z.copy_host2device();

      }

    }

    void prepare(batch<static_molecule> &batch) {
      //Here we should prepare the data needed for the scoring, so basically all the pre processing
      //Basically the idea is that we want to convert easy to use c++ objects into
      //faster buffer that are used in the kernel , so we do the pre processing here and then we just call the kernel in the operator
      batch_ligands                = batch.num_ligands;
      batch_atoms                   = batch.batch_max_atoms;

      //Loading of standard data needed in multiple stages, to avoid boilerplate I have useful helpers
      load_num_rotamers(batch, this->scratch);
      load_num_atoms(batch, this->scratch);

      //Allocation of vinardo specific data
      const int protein_atoms      = static_cast<int>(protein_vinardo.num_atoms());
      const int max_pl_pairs       = batch_ligands * protein_atoms * batch_atoms;
      const int max_ll_pairs       = batch_ligands * batch_atoms * batch_atoms;
      const int scores_per_ligand = std::max(1, static_cast<int>((*this->scratch).configuration.population_number));
      //I need to retrieve the coordinates, that may be modified due to the genetics
      load_scratchs<queue_type>(batch, this->scratch, scores_per_ligand);
      auto& score_b = (*this->scratch).template get<buffer_data_type::SCORES>();

      if (!score_b.is_valid() || score_b.num_elements() != static_cast<std::size_t>(batch_ligands * scores_per_ligand)) {
        score_b.alloc(batch_ligands * scores_per_ligand);
        score_b.set_valid();
      }

      vinardo_num_tors.alloc(batch_ligands);
      inter_scores.alloc(batch_ligands * scores_per_ligand);
      intra_scores.alloc(batch_ligands * scores_per_ligand);

      pl_offsets.alloc(batch_ligands);             //It tells me for ligand k where the pl pairs start
      pl_counts.alloc(batch_ligands);              //It tells me for ligand k how many pl pairs there are starting from offset
      pl_protein_atom_idx.alloc(max_pl_pairs);      //For every pl pair it tells me the index of the protein atom involved
      pl_ligand_atom_idx.alloc(max_pl_pairs);       //For every pl pair it tells me the index of the ligand atom involved
      pl_radius_sum.alloc(max_pl_pairs);            //For every pl pair it tells me the sum of the radii of the atoms involved
      pl_hydrophobic_possible.alloc(max_pl_pairs); //For every pl pair it tells me if an hydrophobic interaction is possible
      pl_hbond_possible.alloc(max_pl_pairs);       //For every pl pair it tells me if an hbond interaction is possible

      ll_offsets.alloc(batch_ligands);            //It tells me for ligand k where the ll pairs start
      ll_counts.alloc(batch_ligands);             //It tells me for ligand k how many ll pairs there are starting from offset
      ll_atom_i_idx.alloc(max_ll_pairs);           //For every ll pair it tells me the index of the first ligand atom involved
      ll_atom_j_idx.alloc(max_ll_pairs);           //For every ll pair it tells me the index of the second ligand atom involved
      ll_radius_sum.alloc(max_ll_pairs);           //For every ll pair it tells me the sum of the radii of the atoms involved
      ll_hydrophobic_possible.alloc(max_ll_pairs); //For every ll pair it tells me if an hydrophobic interaction is possible
      ll_hbond_possible.alloc(max_ll_pairs);       //For every ll pair it tells me if an hbond interaction is possible

      //I need to fill all ligand data
      int current_pl_offset = 0;
      int current_ll_offset = 0;
      for (int ligand_index = 0; ligand_index < batch_ligands; ++ligand_index) {
        //For every ligand I create the ligand layer
        auto& ligand = *batch.molecules[ligand_index];
        vinardo_ligand ligand_vinardo{ligand}; //I create the vinardo ligand that does pre processing in the constructor

        //I write the num tors for the given ligand
        vinardo_num_tors()[ligand_index] = ligand_vinardo.get_num_tors();

        //I need to load all the pl pairs of the ligand, to do it I take the data from vinardo ligand and put them in the buffers
        const auto pl_pairs  = preprocess_protein_ligand_vinardo(ligand_vinardo, protein_vinardo);
        const auto& ll_pairs = ligand_vinardo.get_ligand_ligand_pairs();
        //I update offset and count for the pl pairs
        pl_offsets()[ligand_index] = current_pl_offset;
        pl_counts()[ligand_index]  = static_cast<int>(pl_pairs.size());

        ll_offsets()[ligand_index] = current_ll_offset;
        ll_counts()[ligand_index]  = static_cast<int>(ll_pairs.size());

        for (const auto& pair: pl_pairs) {
          pl_protein_atom_idx()[current_pl_offset] = pair.protein_atom_idx;
          pl_ligand_atom_idx()[current_pl_offset]  = pair.ligand_atom_idx;
          pl_radius_sum()[current_pl_offset]       = pair.radius_sum;
          pl_hydrophobic_possible()[current_pl_offset] = pair.hydrophobic_possible;
          pl_hbond_possible()[current_pl_offset]        = pair.hbond_possible;
          ++current_pl_offset;
        }
        for (const auto& pair: ll_pairs) {
          ll_atom_i_idx()[current_ll_offset] = pair.ligand_atom_i_idx;
          ll_atom_j_idx()[current_ll_offset] = pair.ligand_atom_j_idx;
          ll_radius_sum()[current_ll_offset] = pair.radius_sum;
          ll_hydrophobic_possible()[current_ll_offset] = pair.hydrophobic_possible;
          ll_hbond_possible()[current_ll_offset]        = pair.hbond_possible;
          ++current_ll_offset;
        }
      }

      //move to device
      vinardo_num_tors.copy_host2device();

      pl_offsets.copy_host2device();
      pl_counts.copy_host2device();
      pl_protein_atom_idx.copy_host2device();
      pl_ligand_atom_idx.copy_host2device();
      pl_radius_sum.copy_host2device();
      pl_hydrophobic_possible.copy_host2device();
      pl_hbond_possible.copy_host2device();

      ll_offsets.copy_host2device();
      ll_counts.copy_host2device();
      ll_atom_i_idx.copy_host2device();
      ll_atom_j_idx.copy_host2device();
      ll_radius_sum.copy_host2device();
      ll_hydrophobic_possible.copy_host2device();
      ll_hbond_possible.copy_host2device();

      //create the kernel that will use these data on the device
      //retrieve all the pointers
      const fp_type* x_scratch_b = (*this->scratch).template get<buffer_data_type::X_SCRATCH>().dev_pointer();
      const fp_type* y_scratch_b = (*this->scratch).template get<buffer_data_type::Y_SCRATCH>().dev_pointer();
      const fp_type* z_scratch_b = (*this->scratch).template get<buffer_data_type::Z_SCRATCH>().dev_pointer();

      const fp_type* prot_x_b = device_scratch->template get<buffer_data_type::PROT_X>().dev_pointer();
      const fp_type* prot_y_b = device_scratch->template get<buffer_data_type::PROT_Y>().dev_pointer();
      const fp_type* prot_z_b = device_scratch->template get<buffer_data_type::PROT_Z>().dev_pointer();

      kernel = std::make_unique<vinardo_score_kernel<queue_type>>(
          scores_per_ligand,
          batch_ligands,
          batch_atoms,
          x_scratch_b,
          y_scratch_b,
          z_scratch_b,
          prot_x_b,
          prot_y_b,
          prot_z_b,
          vinardo_num_tors.dev_pointer(),
          pl_offsets.dev_pointer(),
          pl_counts.dev_pointer(),
          pl_protein_atom_idx.dev_pointer(),
          pl_ligand_atom_idx.dev_pointer(),
          pl_radius_sum.dev_pointer(),
          pl_hydrophobic_possible.dev_pointer(),
          pl_hbond_possible.dev_pointer(),
          ll_offsets.dev_pointer(),
          ll_counts.dev_pointer(),
          ll_atom_i_idx.dev_pointer(),
          ll_atom_j_idx.dev_pointer(),
          ll_radius_sum.dev_pointer(),
          ll_hydrophobic_possible.dev_pointer(),
          ll_hbond_possible.dev_pointer(),
          inter_scores.dev_pointer(),
          intra_scores.dev_pointer(),
          score_b.dev_pointer(),
          (*this->scratch).get_queue());
    }

    void operator()() {
      assert(
          (((*this->scratch).template get<buffer_data_type::SCORES>().num_elements() % batch_ligands) == 0) &&
          "Number of scores is not a multiple of ligands in the batch");
      assert(kernel && "Kernel method not yet prepared");
      (*kernel)();
    }

    static std::size_t get_shared_ligand_mem(const int max_atoms, const knobs conf) {
      const int scores_per_ligand = std::max(1, static_cast<int>(conf.population_number));
      return sizeof(int) + sizeof(int) + sizeof(fp_type) * scores_per_ligand +
             3 * sizeof(fp_type) * max_atoms * scores_per_ligand;
    }

    //Now we assume the worst case of a protein with 1024 atoms in order to estimate the memory needed
    static std::size_t get_private_ligand_mem(const int max_atoms, const knobs conf) {
      const int scores_per_ligand = std::max(1, static_cast<int>(conf.population_number));
      constexpr int max_protein_atoms = 1024;
      const int max_pl_pairs = max_protein_atoms * max_atoms;
      const int max_ll_pairs = max_atoms * max_atoms;
      std::size_t mem{0};

      mem += sizeof(int); // vinardo_num_tors

      mem += sizeof(int);                        // pl_offsets
      mem += sizeof(int);                        // pl_counts
      mem += sizeof(int) * max_pl_pairs;         // pl_protein_atom_idx
      mem += sizeof(int) * max_pl_pairs;         // pl_ligand_atom_idx
      mem += sizeof(fp_type) * max_pl_pairs;     // pl_radius_sum
      mem += sizeof(std::uint8_t) * max_pl_pairs; // pl_hydrophobic_possible
      mem += sizeof(std::uint8_t) * max_pl_pairs; // pl_hbond_possible

      mem += sizeof(int);                            // ll_offsets
      mem += sizeof(int);                            // ll_counts
      mem += sizeof(int) * max_ll_pairs;             // ll_atom_i_idx
      mem += sizeof(int) * max_ll_pairs;             // ll_atom_j_idx
      mem += sizeof(fp_type) * max_ll_pairs;         // ll_radius_sum
      mem += sizeof(std::uint8_t) * max_ll_pairs; // ll_hydrophobic_possible
      mem += sizeof(std::uint8_t) * max_ll_pairs;    // ll_hbond_possible

      mem += 2 * sizeof(fp_type) * scores_per_ligand; // inter_scores, intra_scores

      return mem;
    }

    static int get_ligand_mem(const int max_atoms, const knobs conf) {
      return static_cast<int>(get_shared_ligand_mem(max_atoms, conf) +
                              get_private_ligand_mem(max_atoms, conf));
    }
    static batch_multiple get_batch_size(const int max_atoms,
                                         std::shared_ptr<queue_type> q,
                                         const knobs& conf,
                                         const size_t max_bucket_size) {
      (void) max_atoms;
      (void) q;
      (void) conf;
      return batch_multiple{static_cast<int>(std::max<std::size_t>(1, max_bucket_size)), 1};
    }
  private:
    int batch_ligands = 0; //Ligands in the batch
    int batch_atoms   = 0; //Atoms in the batch

    //Protein-Ligand Batch Data
    buffer_vector<int, queue_type> pl_offsets; //pl_offsets[k] is the offset in the pl arrays for the k-th ligand
    buffer_vector<int, queue_type> pl_counts;  //pl_counts[k] is the number of pair of the kth ligand
    //Using the same index on both return the index of the atoms of the couple, and for the others, the properties
    buffer_vector<int, queue_type> pl_protein_atom_idx;
    buffer_vector<int, queue_type> pl_ligand_atom_idx;
    buffer_vector<fp_type, queue_type> pl_radius_sum;
    buffer_vector<std::uint8_t, queue_type> pl_hydrophobic_possible;
    buffer_vector<std::uint8_t, queue_type> pl_hbond_possible;
    //Ligand-Ligand Data
    buffer_vector<int, queue_type> ll_offsets;
    buffer_vector<int, queue_type> ll_counts;
    buffer_vector<int, queue_type> ll_atom_i_idx;
    buffer_vector<int, queue_type> ll_atom_j_idx;
    buffer_vector<fp_type, queue_type> ll_radius_sum;
    buffer_vector<std::uint8_t, queue_type> ll_hydrophobic_possible;
    buffer_vector<std::uint8_t, queue_type> ll_hbond_possible;
    buffer_vector<int, queue_type> vinardo_num_tors;
    buffer_vector<fp_type, queue_type> inter_scores;
    buffer_vector<fp_type, queue_type> intra_scores;

    vinardo_protein protein_vinardo;
    std::shared_ptr<scratchpad<queue_type>> device_scratch;
    std::unique_ptr<vinardo_score_kernel<queue_type>> kernel;

    void teardown_impl(batch<static_molecule> &batch) override {
      assert(batch.num_ligands == batch_ligands && "Scoring algorithm received different batch for teardown");

      auto &scores_b              = (*this->scratch).template get<buffer_data_type::SCORES>();
      const int scores_per_ligand = static_cast<int>(scores_b.num_elements() / batch_ligands);
      scores_b.copy_device2host();
      (*this->scratch).get_queue()->synchronize();
      for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
        auto &ligand          = *batch.molecules[ligand_index];
        const int score_index = ligand_index * scores_per_ligand;
        ligand.properties.assign(property_type::SCORE, std::to_string(scores_b()[score_index]));
      }
    };
  };
#endif
} // namespace mudock
