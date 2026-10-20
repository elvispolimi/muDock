#pragma once

#include <cstddef>
#include <cstring>
#include <mudock/batch.hpp>
#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/chem/autodock_ligand.hpp>
#include <mudock/chem/autodock_protein.hpp>
#include <mudock/compute/adt_score_kernel.hpp>
#include <mudock/compute/adt_score_soa.hpp>
#include <mudock/compute/batch_multiple.hpp>
#include <mudock/compute/scratchpad_layout_sets.hpp>
#if !defined(__CUDACC__) && !defined(__HIPCC__)
  #include <mudock/compute/buffer_utils.hpp>
  #include <mudock/compute/scoring.hpp>
  #include <mudock/compute/scratchpad.hpp>
#endif
#include <mudock/log.hpp>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {

  template<typename queue_type, typename layout_set = all_scratch_layouts>
  batch_multiple get_adt_score_batch_multiple(const int, std::shared_ptr<queue_type>) {
    return {};
  }

#if !defined(__CUDACC__) && !defined(__HIPCC__)
  // TODO check that the object type and the kernel impl are the same
  template<typename queue_type, typename layout_set = all_scratch_layouts>
  struct adt_score: public scoring<queue_type, layout_set> {
    static constexpr const char stage_name[] = "ADT";

    adt_score(std::shared_ptr<scratchpad<queue_type, layout_set>> _scratch,
          std::shared_ptr<scratchpad<queue_type, layout_set>> _device_scratch,
              dynamic_molecule &protein)
        : scoring<queue_type, layout_set>(_scratch),
          atom_data(_scratch->get_queue()),
          num_nonbond_data(_scratch->get_queue()),
          nonbond_data(_scratch->get_queue()),
          device_scratch(_device_scratch) {
      if (!device_scratch->protein_maps().is_valid()) {
        autodock_protein adt_prot(protein);

        auto &prot_bounds = device_scratch->protein_bounds();
        auto &prot_dimensions = device_scratch->protein_dimensions();
        auto &prot_grid_maps = device_scratch->protein_maps();

        prot_bounds.alloc(3);
        prot_dimensions.alloc(1);
        prot_grid_maps.alloc(adt_prot.get_size_xyz() * num_autodock_grids());
        auto bounds = prot_bounds.host_view();
        auto dimensions = prot_dimensions.host_view();

        std::memcpy(bounds.minimum().data(), adt_prot.get_min_p(), 3 * sizeof(fp_type));
        std::memcpy(bounds.maximum().data(), adt_prot.get_max_p(), 3 * sizeof(fp_type));
        std::memcpy(bounds.center().data(), adt_prot.get_center_p(), 3 * sizeof(fp_type));
        dimensions.size_x()[0]   = static_cast<int>(adt_prot.get_size_x());
        dimensions.size_xy()[0]  = static_cast<int>(adt_prot.get_size_xy());
        dimensions.size_xyz()[0] = static_cast<int>(adt_prot.get_size_xyz());
        std::memcpy(prot_grid_maps.host_view().maps().data(),
                    adt_prot.get_maps_pointer(),
                    adt_prot.get_map_flat_size() * num_autodock_grids() * sizeof(fp_type));

        prot_bounds.copy_host2device();
        prot_dimensions.copy_host2device();
        prot_grid_maps.copy_host2device();
        // Workers use this device-shared pointer from their own queues.
        device_scratch->get_queue()->synchronize();
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

      auto &score_data = this->scratch->scores();
      if (load_scratchs<queue_type>(batch, this->scratch, scores_per_ligand)) {
        if (!score_data.is_valid() ||
            score_data.num_elements() != static_cast<size_t>(batch_ligands * scores_per_ligand)) {
          score_data.alloc(batch_ligands * scores_per_ligand);
          score_data.set_valid();
        }
      } else if (!score_data.is_valid() ||
                 score_data.num_elements() != static_cast<size_t>(batch_ligands * scores_per_ligand)) {
        score_data.alloc(batch_ligands * scores_per_ligand);
        score_data.set_valid();
      }

      atom_data.alloc(tot_atoms_in_batch);
      num_nonbond_data.alloc(batch_ligands + 1);
      nonbond_data.alloc(batch_non_bonds);
      auto atom_view = atom_data.host_view();
      auto num_nonbond_view = num_nonbond_data.host_view();
      auto nonbond_view = nonbond_data.host_view();
      num_nonbond_view.num_nonbond()[0] = 0;

        const int map_flat_size = device_scratch->protein_dimensions().host_view().size_xyz()[0];
      for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
        auto &ligand = *batch.molecules[ligand_index];

        autodock_ligand adt_ligand{ligand};

        adt_ligand.update_offsets(map_flat_size);
        const int stride_atoms = ligand_index * batch_atoms;
        // Atoms and bonds
        const int num_atoms = ligand.num_atoms();

        const auto non_bond_size = adt_ligand.non_bond_size();
        std::memcpy(nonbond_view.a1().data() + num_nonbond_view.num_nonbond()[ligand_index],
                    adt_ligand.non_bond_A(),
                    non_bond_size * sizeof(int));
        std::memcpy(nonbond_view.a2().data() + num_nonbond_view.num_nonbond()[ligand_index],
                    adt_ligand.non_bond_B(),
                    non_bond_size * sizeof(int));
        std::memcpy(nonbond_view.cA().data() + num_nonbond_view.num_nonbond()[ligand_index],
                    adt_ligand.non_bond_cA(),
                    non_bond_size * sizeof(fp_type));
        std::memcpy(nonbond_view.cB().data() + num_nonbond_view.num_nonbond()[ligand_index],
                    adt_ligand.non_bond_cB(),
                    non_bond_size * sizeof(fp_type));
        std::memcpy(nonbond_view.xB().data() + num_nonbond_view.num_nonbond()[ligand_index],
                    adt_ligand.non_bond_xB(),
                    non_bond_size * sizeof(int));
        num_nonbond_view.num_nonbond()[ligand_index + 1] =
          static_cast<int>(num_nonbond_view.num_nonbond()[ligand_index] + non_bond_size);

        // Autodock typing
        std::memcpy(atom_view.vols().data() + stride_atoms, adt_ligand.vol(), num_atoms * sizeof(fp_type));
        std::memcpy(atom_view.solpars().data() + stride_atoms, adt_ligand.solpar(), num_atoms * sizeof(fp_type));
        std::memcpy(atom_view.charges().data() + stride_atoms, ligand.charge(), num_atoms * sizeof(fp_type));

        std::memcpy(atom_view.map_offsets().data() + stride_atoms,
                    adt_ligand.atom_map_offsets(),
                    num_atoms * sizeof(int));
      }

      atom_data.copy_host2device();
      num_nonbond_data.copy_host2device();
      nonbond_data.copy_host2device();

      // TODO ask Davide about this performance
      // Bind to the kernel function
      const auto expected_scratch_elements =
          static_cast<std::size_t>(scores_per_ligand) * static_cast<std::size_t>(tot_atoms_in_batch);
      (void) expected_scratch_elements;
      assert(this->scratch->scratch_coordinates().num_elements() ==
                 expected_scratch_elements &&
             "Number of scores per ligand does not match the allocated coordinates space");

        auto batch_view = this->scratch->batch_info().device_view();
        auto scratch_view = this->scratch->scratch_coordinates().device_view();
        auto atom_device_view = atom_data.device_view();
        auto num_nonbond_device_view = num_nonbond_data.device_view();
        auto nonbond_device_view = nonbond_data.device_view();
        const int *num_atoms_b = batch_view.num_atoms().data();
        const int *num_rotamers_b = batch_view.num_rotamers().data();
        const int *num_nonbonds_b  = num_nonbond_device_view.num_nonbond().data();
        const fp_type *x_scratch_b = scratch_view.x().data();
        const fp_type *y_scratch_b = scratch_view.y().data();
        const fp_type *z_scratch_b = scratch_view.z().data();

        const fp_type *vols_b       = atom_device_view.vols().data();
        const fp_type *solpars_b    = atom_device_view.solpars().data();
        const fp_type *charges_b    = atom_device_view.charges().data();
        const int *map_offsets_b    = atom_device_view.map_offsets().data();
        const int *nonbond_a1_b     = nonbond_device_view.a1().data();
        const int *nonbond_a2_b     = nonbond_device_view.a2().data();
        const fp_type *nonbond_cA_b = nonbond_device_view.cA().data();
        const fp_type *nonbond_cB_b = nonbond_device_view.cB().data();
        const int *nonbond_xB_b     = nonbond_device_view.xB().data();

        auto &prot_grid_maps = device_scratch->protein_maps();
      const fp_type *grid_maps = (*device_scratch).get_queue()->obj_required()
                       ? prot_grid_maps.device_view().maps().data()
                       : prot_grid_maps.host_view().maps().data();
        auto bounds_view = device_scratch->protein_bounds().device_view();
        auto dimensions_view = device_scratch->protein_dimensions().host_view();
        const fp_type *minimum = bounds_view.minimum().data();
        const fp_type *maximum = bounds_view.maximum().data();
        const fp_type *center  = bounds_view.center().data();
        const int map_index_x = dimensions_view.size_x()[0];
        const int map_index_xy = dimensions_view.size_xy()[0];
        const int map_index_xyz = dimensions_view.size_xyz()[0];

      fp_type *scores_b = score_data.device_view().scores().data();

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
      mudock::stage_bucket_trace("ADT memory device=",
                                 q->get_id(),
                                 " ligands=",
                                 batch_ligands,
                                 " atoms=",
                                 batch_atoms,
                                 " theoretical_bytes_per_ligand=",
                                 get_ligand_mem(batch_atoms, (*this->scratch).configuration),
                                 " queue_allocated_bytes=",
                                 q->allocated_bytes(),
                                 " queue_peak_allocated_bytes=",
                                 q->peak_allocated_bytes(),
                                 " device_shared_queue_bytes=",
                                 device_scratch->get_queue()->allocated_bytes(),
                                 " device_current_bytes=",
                                 q->device_allocated_bytes(),
                                 " device_peak_allocated_bytes=",
                                 q->device_peak_allocated_bytes());
    }

    [[nodiscard]] auto get_device_shared_queue() const { return device_scratch->get_queue(); }

    void operator()() {
      assert(
          ((this->scratch->scores().num_elements() % batch_ligands) == 0) &&
          "Number of scores is not a multiple of ligands in the batch");
      assert(kernel && "Kernel method not yet prepared");
      (*kernel)();
    }

    static std::size_t get_shared_ligand_mem(const int max_atoms, const knobs conf) {
      const int scores_per_ligand = std::max(1, static_cast<int>(conf.population_number));
      return sizeof(int) * 2 + sizeof(fp_type) * scores_per_ligand +
             sizeof(fp_type) * 3 * max_atoms * scores_per_ligand;
    }

    static std::size_t get_private_ligand_mem(const int max_atoms, const knobs) {
      const auto non_bonds_atoms = static_cast<std::size_t>(max_atoms) * max_atoms;
      std::size_t mem{0};
      mem += sizeof(fp_type) * 3 * max_atoms;                         // vols, solpars, charges
      mem += sizeof(int) * max_atoms;                                 // map_offsets
      mem += sizeof(int);                                           // num_nonbond (plus one fixed batch sentinel)
      mem += sizeof(int) * 2 * non_bonds_atoms;                       // nonbond_a1, nonbond_a2
      mem += sizeof(fp_type) * 2 * non_bonds_atoms;                   // nonbond_cA, nonbond_cB
      mem += sizeof(int) * non_bonds_atoms;                           // nonbond_xB
      return mem;
    }

    static std::size_t get_ligand_mem(const int max_atoms, const knobs conf) {
      return get_shared_ligand_mem(max_atoms, conf) + get_private_ligand_mem(max_atoms, conf);
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

    soa_buffer<AdtAtomData<>, queue_type> atom_data;
    soa_buffer<AdtNumNonbondData<>, queue_type> num_nonbond_data;
    soa_buffer<AdtNonbondData<>, queue_type> nonbond_data;

    std::shared_ptr<scratchpad<queue_type, layout_set>> device_scratch;
    std::unique_ptr<adt_score_kernel<queue_type>> kernel;

    void teardown_impl(batch<static_molecule> &batch) override {
      assert(batch.num_ligands == batch_ligands && "Scoring algorithm received different batch for teardown");

      auto &score_data = this->scratch->scores();
      const int scores_per_ligand = static_cast<int>(score_data.num_elements() / batch_ligands);
      score_data.copy_device2host();
      (*this->scratch).get_queue()->synchronize();
      auto score_view = score_data.host_view();
      for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
        auto &ligand          = *batch.molecules[ligand_index];
        const int score_index = ligand_index * scores_per_ligand;
        ligand.properties.assign(property_type::SCORE, std::to_string(score_view.scores()[score_index]));
      }
    };
  };
#endif
} // namespace mudock
