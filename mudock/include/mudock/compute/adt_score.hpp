#pragma once

#include <cstddef>
#include <cstring>
#include <type_traits>
#include <mudock/batch.hpp>
#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/chem/autodock_ligand.hpp>
#include <mudock/chem/autodock_protein.hpp>
#include <mudock/compute/adt_score_kernel.hpp>
#include <mudock/data_formats/adt_score_soa.hpp>
#include <mudock/compute/batch_multiple.hpp>
#include <mudock/data_formats/data_layout_sets.hpp>
#include <mudock/compute/ligand_placement.hpp>
#if !defined(__CUDACC__) && !defined(__HIPCC__)
  #include <mudock/compute/buffer_utils.hpp>
  #include <mudock/compute/scoring.hpp>
  #include <mudock/data_formats/data_store.hpp>
#endif
#include <mudock/log.hpp>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {

  template<typename queue_type, typename layout_set = scoring_data_layouts>
  batch_multiple get_adt_score_batch_multiple(const int, std::shared_ptr<queue_type>) {
    return {};
  }

#if !defined(__CUDACC__) && !defined(__HIPCC__)
  // TODO check that the object type and the kernel impl are the same
  template<typename queue_type, typename layout_set = scoring_data_layouts,
           typename protein_layout_set = protein_data_layouts>
  struct adt_score: public scoring<queue_type, layout_set> {
    static constexpr const char stage_name[] = "ADT";

    adt_score(std::shared_ptr<data_store<queue_type, layout_set>> _data,
              std::shared_ptr<data_store<queue_type, protein_layout_set>> _protein_data,
              dynamic_molecule &protein,
              ligand_placement _placement = {})
        : scoring<queue_type, layout_set>(_data),
          adt_data(_data->get_queue()),
          protein_data(_protein_data),
          protein_center(protein.get_center()),
          placement(_placement) {
      if (!protein_data->is_valid()) {
        autodock_protein adt_prot(protein);

        protein_data->alloc(3,
                              1,
                              static_cast<int>(adt_prot.get_size_xyz() * num_autodock_grids()));
        auto protein_view = protein_data->host_view();
        auto bounds = protein_view.protein_bounds();
        auto dimensions = protein_view.protein_dimensions();
        auto maps = protein_view.protein_maps();

        std::memcpy(bounds.minimum().data(), adt_prot.get_min_p(), 3 * sizeof(fp_type));
        std::memcpy(bounds.maximum().data(), adt_prot.get_max_p(), 3 * sizeof(fp_type));
        std::memcpy(bounds.center().data(), adt_prot.get_center_p(), 3 * sizeof(fp_type));
        dimensions.size_x()[0]   = static_cast<int>(adt_prot.get_size_x());
        dimensions.size_xy()[0]  = static_cast<int>(adt_prot.get_size_xy());
        dimensions.size_xyz()[0] = static_cast<int>(adt_prot.get_size_xyz());
        std::memcpy(maps.maps().data(),
                    adt_prot.get_maps_pointer(),
                    adt_prot.get_map_flat_size() * num_autodock_grids() * sizeof(fp_type));

        protein_data->copy_host2device();
        // Workers use this device-shared pointer from their own queues.
        protein_data->get_queue()->synchronize();
      }
    }

    void prepare(batch<static_molecule> &batch) {
      batch_ligands                = batch.num_ligands;
      batch_atoms                  = batch.batch_max_atoms;
      const int batch_non_bonds    = batch_ligands * batch_atoms * batch_atoms;
      const int tot_atoms_in_batch = batch_ligands * batch_atoms;
      auto q                       = (*this->data).get_queue();

      for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index)
        // Correct protein center putting the real data structure
        apply_ligand_placement(*batch.molecules[ligand_index], placement, protein_center);

      const int scores_per_ligand =
          std::max(1, static_cast<int>((*this->data).configuration.population_number));
      if constexpr (std::is_same_v<layout_set, scoring_data_layouts>) {
        this->data->alloc(tot_atoms_in_batch * scores_per_ligand,
                             batch_ligands,
                             batch_ligands * scores_per_ligand);
      }

      load_num_rotamers(batch, this->data);
      load_num_atoms(batch, this->data);
      load_scoring_coordinates<queue_type>(batch, this->data, scores_per_ligand);
      this->data->copy_host2device();

      adt_data.alloc(tot_atoms_in_batch, batch_ligands + 1, batch_non_bonds);
      auto adt_view = adt_data.host_view();
      auto atom_view = adt_view.atom_data();
      auto num_nonbond_view = adt_view.num_nonbond_data();
      auto nonbond_view = adt_view.nonbond_data();
      num_nonbond_view.num_nonbond()[0] = 0;

      auto protein_host_view = protein_data->host_view();
      const int map_flat_size = protein_host_view.protein_dimensions().size_xyz()[0];
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

      adt_data.copy_host2device();

      // TODO ask Davide about this performance
      // Bind to the kernel function
      const auto expected_data_elements =
          static_cast<std::size_t>(scores_per_ligand) * static_cast<std::size_t>(tot_atoms_in_batch);
      (void) expected_data_elements;
      assert(this->data->host_view().scoring_coordinates().metadata().size() ==
                 expected_data_elements &&
             "Number of scores per ligand does not match the allocated coordinates space");

        auto batch_view = this->data->device_view().batch_info();
        auto  data_view = this->data->device_view().scoring_coordinates();
        auto adt_device_view = adt_data.device_view();
        auto atom_device_view = adt_device_view.atom_data();
        auto num_nonbond_device_view = adt_device_view.num_nonbond_data();
        auto nonbond_device_view = adt_device_view.nonbond_data();
        const int *num_atoms_b = batch_view.num_atoms().data();
        const int *num_rotamers_b = batch_view.num_rotamers().data();
        const int *num_nonbonds_b  = num_nonbond_device_view.num_nonbond().data();
        const fp_type *x_data_b =  data_view.x().data();
        const fp_type *y_data_b =  data_view.y().data();
        const fp_type *z_data_b =  data_view.z().data();

        const fp_type *vols_b       = atom_device_view.vols().data();
        const fp_type *solpars_b    = atom_device_view.solpars().data();
        const fp_type *charges_b    = atom_device_view.charges().data();
        const int *map_offsets_b    = atom_device_view.map_offsets().data();
        const int *nonbond_a1_b     = nonbond_device_view.a1().data();
        const int *nonbond_a2_b     = nonbond_device_view.a2().data();
        const fp_type *nonbond_cA_b = nonbond_device_view.cA().data();
        const fp_type *nonbond_cB_b = nonbond_device_view.cB().data();
        const int *nonbond_xB_b     = nonbond_device_view.xB().data();

      auto protein_device_view = protein_data->device_view();
      const fp_type *grid_maps = (*protein_data).get_queue()->obj_required()
                       ? protein_device_view.protein_maps().maps().data()
                       : protein_host_view.protein_maps().maps().data();
        auto bounds_view = protein_device_view.protein_bounds();
        auto dimensions_view = protein_host_view.protein_dimensions();
        const fp_type *minimum = bounds_view.minimum().data();
        const fp_type *maximum = bounds_view.maximum().data();
        const fp_type *center  = bounds_view.center().data();
        const int map_index_x = dimensions_view.size_x()[0];
        const int map_index_xy = dimensions_view.size_xy()[0];
        const int map_index_xyz = dimensions_view.size_xyz()[0];

      fp_type *scores_b = this->data->device_view().scores().scores().data();

      kernel = std::make_unique<adt_score_kernel<queue_type>>(scores_per_ligand,
                                                              batch_ligands,
                                                              batch_atoms,
                                                              num_atoms_b,
                                                              num_rotamers_b,
                                                              num_nonbonds_b,
                                                              x_data_b,
                                                              y_data_b,
                                                              z_data_b,
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

    [[nodiscard]] auto get_device_shared_queue() const { return protein_data->get_queue(); }

    void operator()() {
      assert(
          ((this->data->host_view().scores().metadata().size() % batch_ligands) == 0) &&
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
      const std::size_t non_bonds_atoms = static_cast<std::size_t>(max_atoms) *
                                          static_cast<std::size_t>(max_atoms);
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

    soa_buffer<AdtScoreData<>, queue_type> adt_data;

    std::shared_ptr<data_store<queue_type, protein_layout_set>> protein_data;
    point<fp_type, 3> protein_center;
    ligand_placement placement;
    std::unique_ptr<adt_score_kernel<queue_type>> kernel;

    void teardown_impl(batch<static_molecule> &batch) override {
      assert(batch.num_ligands == batch_ligands && "Scoring algorithm received different batch for teardown");

      const int scores_per_ligand = static_cast<int>(
          this->data->host_view().scores().metadata().size() / batch_ligands);
      // TODO is it copying everything?? We need only scores
      this->data->copy_device2host_scores();
      (*this->data).get_queue()->synchronize();
      auto score_view = this->data->host_view().scores();
      for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
        auto &ligand          = *batch.molecules[ligand_index];
        const int score_index = ligand_index * scores_per_ligand;
        ligand.properties.assign(property_type::SCORE, std::to_string(score_view.scores()[score_index]));
      }
    };
  };
#endif
} // namespace mudock
