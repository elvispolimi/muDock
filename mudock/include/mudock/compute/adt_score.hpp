#pragma once

#include <cstddef>
#include <cstring>
#include <mudock/batch.hpp>
#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/chem/autodock_ligand.hpp>
#include <mudock/chem/autodock_protein.hpp>
#include <mudock/compute/adt_score_kernel.hpp>
#include <mudock/compute/batch_multiple.hpp>
#if !defined(__CUDACC__) && !defined(__HIPCC__)
  #include <mudock/compute/buffer_utils.hpp>
  #include <mudock/compute/scoring.hpp>
  #include <mudock/compute/scratchpad.hpp>
#endif
#include <mudock/log.hpp>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {

  template<typename queue_type>
  batch_multiple get_adt_score_batch_multiple(const int, std::shared_ptr<queue_type>) {
    return {};
  }

#if !defined(__CUDACC__) && !defined(__HIPCC__)
  // TODO check that the object type and the kernel impl are the same
  template<typename queue_type>
  struct adt_score: public scoring<queue_type> {
    static constexpr const char stage_name[] = "ADT";

    adt_score(std::shared_ptr<scratchpad<queue_type>> _scratch,
              std::shared_ptr<scratchpad<queue_type>> _device_scratch,
              dynamic_molecule &protein)
        : scoring<queue_type>(_scratch),
          vols(_scratch->get_queue()),
          solpars(_scratch->get_queue()),
          charges(_scratch->get_queue()),
          map_offsets(_scratch->get_queue()),
          num_nonbond(_scratch->get_queue()),
          nonbond_a1(_scratch->get_queue()),
          nonbond_a2(_scratch->get_queue()),
          nonbond_cA(_scratch->get_queue()),
          nonbond_cB(_scratch->get_queue()),
          nonbond_xB(_scratch->get_queue()),
          device_scratch(_device_scratch) {
      if (!(*device_scratch).template exists<buffer_data_type::PROT_GRID_MAPS>()) {
        autodock_protein adt_prot(protein);

        auto &prot_min       = (*device_scratch).template get<buffer_data_type::PROT_MIN>();
        auto &prot_max       = (*device_scratch).template get<buffer_data_type::PROT_MAX>();
        auto &prot_center    = (*device_scratch).template get<buffer_data_type::PROT_CENTER>();
        auto &prot_index_x   = (*device_scratch).template get<buffer_data_type::PROT_SIZE_X>();
        auto &prot_index_xy  = (*device_scratch).template get<buffer_data_type::PROT_SIZE_XY>();
        auto &prot_index_xyz = (*device_scratch).template get<buffer_data_type::PROT_SIZE_XYZ>();
        auto &prot_grid_maps = (*device_scratch).template get<buffer_data_type::PROT_GRID_MAPS>();

        prot_min.alloc(3);
        prot_max.alloc(3);
        prot_center.alloc(3);
        prot_index_x.alloc(1);
        prot_index_xy.alloc(1);
        prot_index_xyz.alloc(1);
        prot_grid_maps.alloc(adt_prot.get_size_xyz() * num_autodock_grids());

        std::memcpy(prot_min(), adt_prot.get_min_p(), 3 * sizeof(fp_type));
        std::memcpy(prot_max(), adt_prot.get_max_p(), 3 * sizeof(fp_type));
        std::memcpy(prot_center(), adt_prot.get_center_p(), 3 * sizeof(fp_type));
        prot_index_x()[0]   = static_cast<int>(adt_prot.get_size_x());
        prot_index_xy()[0]  = static_cast<int>(adt_prot.get_size_xy());
        prot_index_xyz()[0] = static_cast<int>(adt_prot.get_size_xyz());
        std::memcpy(prot_grid_maps(),
                    adt_prot.get_maps_pointer(),
                    adt_prot.get_map_flat_size() * num_autodock_grids() * sizeof(fp_type));

        prot_min.copy_host2device();
        prot_max.copy_host2device();
        prot_center.copy_host2device();
        prot_index_x.copy_host2device();
        prot_index_xy.copy_host2device();
        prot_index_xyz.copy_host2device();
        // On CPU is not required and on GPUS we have probably to laod texture memory etc...
        // prot_grid_maps.copy_host2device();
      }
    }

    void prepare(batch<static_molecule> &batch) {
      batch_ligands                = batch.num_ligands;
      batch_atoms                  = batch.batch_max_atoms;
      const int batch_non_bonds    = batch_ligands * batch_atoms * batch_atoms;
      const int tot_atoms_in_batch = batch_ligands * batch_atoms;
      auto q                       = (*this->scratch).get_queue();

      load_num_rotamers(batch, this->scratch);
      load_num_atoms(batch, this->scratch);

      const int scores_per_ligand =
          std::max(1, static_cast<int>((*this->scratch).configuration.population_number));

      auto &score_b = (*this->scratch).template get<buffer_data_type::SCORES>();
      if (load_scratchs<queue_type>(batch, this->scratch, scores_per_ligand)) {
        if (!score_b.is_valid() ||
            score_b.num_elements() != static_cast<size_t>(batch_ligands * scores_per_ligand)) {
          score_b.alloc(batch_ligands * scores_per_ligand);
          score_b.set_valid();
        }
      } else if (!score_b.is_valid() ||
                 score_b.num_elements() != static_cast<size_t>(batch_ligands * scores_per_ligand)) {
        score_b.alloc(batch_ligands * scores_per_ligand);
        score_b.set_valid();
      }

      vols.alloc(tot_atoms_in_batch);
      solpars.alloc(tot_atoms_in_batch);
      charges.alloc(tot_atoms_in_batch);
      map_offsets.alloc(tot_atoms_in_batch);
      num_nonbond.alloc(batch_ligands + 1);
      num_nonbond()[0] = 0;
      nonbond_a1.alloc(batch_non_bonds);
      nonbond_a2.alloc(batch_non_bonds);
      nonbond_cA.alloc(batch_non_bonds);
      nonbond_cB.alloc(batch_non_bonds);
      nonbond_xB.alloc(batch_non_bonds);

      const int map_flat_size =
          (*device_scratch).template get<buffer_data_type::PROT_SIZE_XYZ>().host_pointer()[0];
      for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
        auto &ligand = *batch.molecules[ligand_index];

        autodock_ligand adt_ligand{ligand};

        adt_ligand.update_offsets(map_flat_size);
        const int stride_atoms = ligand_index * batch_atoms;
        // Atoms and bonds
        const int num_atoms = ligand.num_atoms();

        const auto non_bond_size = adt_ligand.non_bond_size();
        std::memcpy((void *) (nonbond_a1() + num_nonbond()[ligand_index]),
                    adt_ligand.non_bond_A(),
                    non_bond_size * sizeof(int));
        std::memcpy((void *) (nonbond_a2() + num_nonbond()[ligand_index]),
                    adt_ligand.non_bond_B(),
                    non_bond_size * sizeof(int));
        std::memcpy((void *) (nonbond_cA() + num_nonbond()[ligand_index]),
                    adt_ligand.non_bond_cA(),
                    non_bond_size * sizeof(fp_type));
        std::memcpy((void *) (nonbond_cB() + num_nonbond()[ligand_index]),
                    adt_ligand.non_bond_cB(),
                    non_bond_size * sizeof(fp_type));
        std::memcpy((void *) (nonbond_xB() + num_nonbond()[ligand_index]),
                    adt_ligand.non_bond_xB(),
                    non_bond_size * sizeof(int));
        num_nonbond()[ligand_index + 1] = static_cast<int>(num_nonbond()[ligand_index] + non_bond_size);

        // Autodock typing
        std::memcpy((void *) (vols() + stride_atoms), adt_ligand.vol(), num_atoms * sizeof(fp_type));
        std::memcpy((void *) (solpars() + stride_atoms), adt_ligand.solpar(), num_atoms * sizeof(fp_type));
        std::memcpy((void *) (charges() + stride_atoms), ligand.charge(), num_atoms * sizeof(fp_type));

        std::memcpy((void *) (map_offsets() + stride_atoms),
                    adt_ligand.atom_map_offsets(),
                    num_atoms * sizeof(int));
      }

      vols.copy_host2device();
      solpars.copy_host2device();
      charges.copy_host2device();
      map_offsets.copy_host2device();
      num_nonbond.copy_host2device();
      num_nonbond()[0] = 0;
      nonbond_a1.copy_host2device();
      nonbond_a2.copy_host2device();
      nonbond_cA.copy_host2device();
      nonbond_cB.copy_host2device();
      nonbond_xB.copy_host2device();

      // TODO ask Davide about this performance
      // Bind to the kernel function
      const auto expected_scratch_elements =
          static_cast<std::size_t>(scores_per_ligand) * static_cast<std::size_t>(tot_atoms_in_batch);
      (void) expected_scratch_elements;
      assert((*this->scratch).template get<buffer_data_type::X_SCRATCH>().num_elements() ==
                 expected_scratch_elements &&
             "Number of scores per ligand does not match the allocated coordinates space");

      const int *num_atoms_b = (*this->scratch).template get<buffer_data_type::NUM_ATOMS>().dev_pointer();
      const int *num_rotamers_b =
          (*this->scratch).template get<buffer_data_type::NUM_ROTAMERS>().dev_pointer();
      const int *num_nonbonds_b  = num_nonbond.dev_pointer();
      const fp_type *x_scratch_b = (*this->scratch).template get<buffer_data_type::X_SCRATCH>().dev_pointer();
      const fp_type *y_scratch_b = (*this->scratch).template get<buffer_data_type::Y_SCRATCH>().dev_pointer();
      const fp_type *z_scratch_b = (*this->scratch).template get<buffer_data_type::Z_SCRATCH>().dev_pointer();

      const fp_type *vols_b       = vols.dev_pointer();
      const fp_type *solpars_b    = solpars.dev_pointer();
      const fp_type *charges_b    = charges.dev_pointer();
      const int *map_offsets_b    = map_offsets.dev_pointer();
      const int *nonbond_a1_b     = nonbond_a1.dev_pointer();
      const int *nonbond_a2_b     = nonbond_a2.dev_pointer();
      const fp_type *nonbond_cA_b = nonbond_cA.dev_pointer();
      const fp_type *nonbond_cB_b = nonbond_cB.dev_pointer();
      const int *nonbond_xB_b     = nonbond_xB.dev_pointer();

      // Use host pointer as on CPP you can use it, on GPU they will load their own memory
      const fp_type *grid_maps =
          (*device_scratch).template get<buffer_data_type::PROT_GRID_MAPS>().host_pointer();
      const fp_type *minimum = (*device_scratch).template get<buffer_data_type::PROT_MIN>().dev_pointer();
      const fp_type *maximum = (*device_scratch).template get<buffer_data_type::PROT_MAX>().dev_pointer();
      const fp_type *center  = (*device_scratch).template get<buffer_data_type::PROT_CENTER>().dev_pointer();
      const int map_index_x =
          (*device_scratch).template get<buffer_data_type::PROT_SIZE_X>().host_pointer()[0];
      const int map_index_xy =
          (*device_scratch).template get<buffer_data_type::PROT_SIZE_XY>().host_pointer()[0];
      const int map_index_xyz =
          (*device_scratch).template get<buffer_data_type::PROT_SIZE_XYZ>().host_pointer()[0];

      fp_type *scores_b = score_b.dev_pointer();

      kernel = std::make_unique<adt_score_kernel<queue_type>>(scores_per_ligand,
                                                              batch_ligands,
                                                              batch_atoms,
                                                              num_atoms_b,
                                                              num_rotamers_b,
                                                              num_nonbonds_b,
                                                              x_scratch_b,
                                                              y_scratch_b,
                                                              z_scratch_b,
                                                              vols_b,
                                                              solpars_b,
                                                              charges_b,
                                                              map_offsets_b,
                                                              nonbond_a1_b,
                                                              nonbond_a2_b,
                                                              nonbond_cA_b,
                                                              nonbond_cB_b,
                                                              nonbond_xB_b,
                                                              grid_maps,
                                                              minimum,
                                                              maximum,
                                                              center,
                                                              map_index_x,
                                                              map_index_xy,
                                                              map_index_xyz,
                                                              scores_b,
                                                              q);
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

    static std::size_t get_private_ligand_mem(const int max_atoms, const knobs) {
      const int non_bonds_atoms = max_atoms * max_atoms;
      std::size_t mem{0};
      mem += sizeof(fp_type) * max_atoms;       // vols
      mem += sizeof(fp_type) * max_atoms;       // solpars
      mem += sizeof(fp_type) * max_atoms;       // charges
      mem += sizeof(int) * max_atoms;           // map_offsets
      mem += sizeof(int);                       // num_nonbond
      mem += sizeof(int) * non_bonds_atoms;     // nonbond_a1
      mem += sizeof(int) * non_bonds_atoms;     // nonbond_a2
      mem += sizeof(fp_type) * non_bonds_atoms; // nonbond_cA
      mem += sizeof(fp_type) * non_bonds_atoms; // nonbond_cB
      mem += sizeof(int) * non_bonds_atoms;     // nonbond_xB
      return mem;
    }

    static int get_ligand_mem(const int max_atoms, const knobs conf) {
      return static_cast<int>(get_shared_ligand_mem(max_atoms, conf) +
                              get_private_ligand_mem(max_atoms, conf));
    }

    static batch_multiple get_batch_size(const int atoms,
                                         std::shared_ptr<queue_type> q,
                                         const knobs &conf,
                                         const size_t max_bucket_size) {
      (void) conf;
      const auto plain_multiple_info =
          normalize_batch_multiple(get_adt_score_batch_multiple<queue_type>(atoms, q));
      mudock::stage_bucket_trace("ADT stage plain multiple for ",
                                 atoms,
                                 " atoms -> total=",
                                 plain_multiple_info.total_multiple(),
                                 " (active_blocks_per_sm=",
                                 plain_multiple_info.active_blocks_per_sm,
                                 ", num_sms=",
                                 plain_multiple_info.num_sms,
                                 ")",
                                 " (max_bucket_size hint=",
                                 max_bucket_size,
                                 ")");
      return plain_multiple_info;
    }

  private:
    int batch_ligands;
    int batch_atoms;

    buffer_vector<fp_type, queue_type> vols;
    buffer_vector<fp_type, queue_type> solpars;
    buffer_vector<fp_type, queue_type> charges;
    buffer_vector<int, queue_type> map_offsets;
    buffer_vector<int, queue_type> num_nonbond;
    buffer_vector<int, queue_type> nonbond_a1;
    buffer_vector<int, queue_type> nonbond_a2;
    buffer_vector<fp_type, queue_type> nonbond_cA;
    buffer_vector<fp_type, queue_type> nonbond_cB;
    buffer_vector<int, queue_type> nonbond_xB;

    std::shared_ptr<scratchpad<queue_type>> device_scratch;
    std::unique_ptr<adt_score_kernel<queue_type>> kernel;

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
