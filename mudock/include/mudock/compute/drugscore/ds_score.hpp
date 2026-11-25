#pragma once

#include <concepts>
#include <cstring>
#include <functional>
#include <memory>
#include <mudock/batch.hpp>
#include <mudock/chem/drugscore/drugscore_ligand.hpp>
#include <mudock/chem/drugscore/drugscore_protein.hpp>
#include <mudock/chem/drugscore/drugscore_types.hpp>
#include <mudock/compute/batch_multiple.hpp>
#include <mudock/compute/drugscore/ds_score_kernel.hpp>
#include <mudock/compute/ligand_placement.hpp>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>
#include <variant>

#ifndef __CUDACC__
  #include <mudock/compute/buffer_utils.hpp>
  #include <mudock/compute/scoring.hpp>
  #include <mudock/compute/scratchpad.hpp>
#endif

namespace mudock {

  template<typename queue_type>
  batch_multiple get_ds_score_batch_multiple(const int, std::shared_ptr<queue_type>) {
    return {};
  }

#ifndef __CUDACC__
  template<typename queue_type, bool PRECOMPUTE>
  struct ds_score_impl: public scoring<queue_type> {
    static constexpr const char stage_name[] = "DS";

    ds_score_impl(std::shared_ptr<scratchpad<queue_type>> scratch_,
             std::shared_ptr<scratchpad<queue_type>> device_scratch_,
             dynamic_molecule& protein,
             ligand_placement placement_ = {})
        : scoring<queue_type>(scratch_),
          lig_ds_types(scratch_->get_queue()),
          device_scratch(device_scratch_),
          placement(placement_),
          protein_center(protein.get_center()) {
      if constexpr (PRECOMPUTE) {
        if (!(*device_scratch).template exists<buffer_data_type::PROT_DS_GRID_MAPS>()) {
          drugscore_precomp_protein ds_prot(protein);

          auto& protein_min          = (*device_scratch).template get<buffer_data_type::PROT_MIN>();
          auto& protein_max          = (*device_scratch).template get<buffer_data_type::PROT_MAX>();
          auto& protein_index_x      = (*device_scratch).template get<buffer_data_type::PROT_SIZE_X>();
          auto& protein_index_xy     = (*device_scratch).template get<buffer_data_type::PROT_SIZE_XY>();
          auto& protein_index_xyz    = (*device_scratch).template get<buffer_data_type::PROT_SIZE_XYZ>();
          auto& protein_ds_grid_maps = (*device_scratch).template get<buffer_data_type::PROT_DS_GRID_MAPS>();

          protein_min.alloc(3);
          protein_max.alloc(3);
          protein_index_x.alloc(1);
          protein_index_xy.alloc(1);
          protein_index_xyz.alloc(1);
          protein_ds_grid_maps.alloc(ds_prot.get_map_flat_size() * num_valued_ds_types());

          std::memcpy(protein_min(), ds_prot.get_min_p(), 3 * sizeof(fp_type));
          std::memcpy(protein_max(), ds_prot.get_max_p(), 3 * sizeof(fp_type));
          protein_index_x()[0]   = static_cast<int>(ds_prot.get_size_x());
          protein_index_xy()[0]  = static_cast<int>(ds_prot.get_size_xy());
          protein_index_xyz()[0] = static_cast<int>(ds_prot.get_size_xyz());
          std::memcpy(protein_ds_grid_maps(),
                      ds_prot.get_maps_pointer(),
                      ds_prot.get_map_flat_size() * num_valued_ds_types() * sizeof(fp_type));

          protein_min.copy_host2device();
          protein_max.copy_host2device();
          protein_index_x.copy_host2device();
          protein_index_xy.copy_host2device();
          protein_index_xyz.copy_host2device();
          protein_ds_grid_maps.copy_host2device();
        }
      } else {
        if (!(*device_scratch).template exists<buffer_data_type::PROT_X_COORDS>()) {
          drugscore_protein ds_prot(protein);
          const int n_atoms = protein.num_atoms();

          // sort indices on drugscore types,
          // N.B. FEATURE TO TEST
          // std::vector<unsigned int> indices(n_atoms);
          // std::iota(indices.begin(), indices.end(), 0);
          // std::sort(indices.begin(), indices.end(), [&](const int a, const int b) {
          //   return static_cast<int>(ds_prot.drugscore_type(a)) < static_cast<int>(ds_prot.drugscore_type(b));
          // });
          //ds_prot.permute(indices.data(), n_atoms);

          auto& prot_x_coords  = (*device_scratch).template get<buffer_data_type::PROT_X_COORDS>();
          auto& prot_y_coords  = (*device_scratch).template get<buffer_data_type::PROT_Y_COORDS>();
          auto& prot_z_coords  = (*device_scratch).template get<buffer_data_type::PROT_Z_COORDS>();
          auto& prot_num_atoms = (*device_scratch).template get<buffer_data_type::PROT_NUM_ATOMS>();
          auto& prot_ds_types  = (*device_scratch).template get<buffer_data_type::PROT_DS_TYPES>();

          prot_num_atoms.alloc(1);
          prot_x_coords.alloc(n_atoms);
          prot_y_coords.alloc(n_atoms);
          prot_z_coords.alloc(n_atoms);
          prot_ds_types.alloc(n_atoms);

          std::memcpy(prot_x_coords(), protein.x(), n_atoms * sizeof(fp_type));
          std::memcpy(prot_y_coords(), protein.y(), n_atoms * sizeof(fp_type));
          std::memcpy(prot_z_coords(), protein.z(), n_atoms * sizeof(fp_type));
          std::memcpy(prot_ds_types(), ds_prot.drugscore_type(), n_atoms * sizeof(ds_type));
          prot_num_atoms()[0] = n_atoms;
          prot_x_coords.copy_host2device();
          prot_y_coords.copy_host2device();
          prot_z_coords.copy_host2device();
          prot_ds_types.copy_host2device();
          prot_num_atoms.copy_host2device();
        }

        if (!(*device_scratch).template exists<buffer_data_type::DS_MODEL>()) {
          auto& ds_model = (*device_scratch).template get<buffer_data_type::DS_MODEL>();

          const size_t total_size = num_valued_ds_types() * num_valued_ds_types() * num_buckets();
          ds_model.alloc(total_size);
          std::memcpy(ds_model(), DRUGSCORE_MODEL_FLAT(), total_size * sizeof(fp_type));
          ds_model.copy_host2device();
        }
      }

      device_scratch->get_queue()->synchronize();
    }

    void prepare(batch<static_molecule>& batch) {
      batch_ligands                = batch.num_ligands;
      batch_atoms                  = batch.batch_max_atoms;
      const int tot_atoms_in_batch = batch_ligands * batch_atoms;

      for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index)
        apply_ligand_placement(*batch.molecules[ligand_index], placement, protein_center);


      load_num_atoms(batch, this->scratch);

      // Prepare the scores buffer only if scratches are NOT already loaded
      auto& scores_b = (*this->scratch).template get<buffer_data_type::SCORES>();
      if (load_scratchs<queue_type>(batch, this->scratch)) {
        scores_b.alloc(batch_ligands);
        scores_b.set_valid();
      }

      //-- Transfer ligand drugscore types
      lig_ds_types.alloc(tot_atoms_in_batch);
      for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
        auto& ligand = *batch.molecules[ligand_index];

        drugscore_ligand ds_ligand{ligand};

        const int stride_atoms = ligand_index * batch_atoms;
        const int num_atoms    = ds_ligand.get_base_molecule().num_atoms();

        // sort indices on drugscore types,
        // N.B. FEATURE TO TEST
        // std::vector<unsigned int> indices(num_atoms);
        // std::iota(indices.begin(), indices.end(), 0);
        // std::sort(indices.begin(), indices.end(), [&](const int a, const int b) {
        //   return static_cast<int>(ds_ligand.drugscore_type(a)) <
        //          static_cast<int>(ds_ligand.drugscore_type(b));
        // });
        //ds_ligand.permute(indices.data(), num_atoms);

        std::memcpy((void*) (lig_ds_types() + stride_atoms),
                    ds_ligand.drugscore_type(),
                    num_atoms * sizeof(ds_type));
      }
      lig_ds_types.copy_host2device();

      //-- Bind to the kernel function
      const auto scores_per_ligand = scores_b.num_elements() / batch_ligands;
      assert((*this->scratch).template get<buffer_data_type::X_SCRATCH>().num_elements() !=
                 (scores_per_ligand * batch_ligands) &&
             "Number of scores per ligands does not match the allocated coordinates space");

      // ligand data pointers
      const fp_type* x_lig_p = (*this->scratch).template get<buffer_data_type::X_SCRATCH>().dev_pointer();
      const fp_type* y_lig_p = (*this->scratch).template get<buffer_data_type::Y_SCRATCH>().dev_pointer();
      const fp_type* z_lig_p = (*this->scratch).template get<buffer_data_type::Z_SCRATCH>().dev_pointer();
      const ds_type* lig_ds_types_p = lig_ds_types.dev_pointer();
      const int* lig_n_atm_p = (*this->scratch).template get<buffer_data_type::NUM_ATOMS>().dev_pointer();

      // result scores pointer
      fp_type* scores_p = (*this->scratch).template get<buffer_data_type::SCORES>().dev_pointer();
      // Device queue
      auto q = this->scratch->get_queue();

      // protein data pointers
      if constexpr (PRECOMPUTE) {
        const fp_type* prot_min = (*device_scratch).template get<buffer_data_type::PROT_MIN>().dev_pointer();
        const fp_type* prot_max = (*device_scratch).template get<buffer_data_type::PROT_MAX>().dev_pointer();
        const int prot_index_x =
            (*device_scratch).template get<buffer_data_type::PROT_SIZE_X>().host_pointer()[0];
        const int prot_index_xy =
            (*device_scratch).template get<buffer_data_type::PROT_SIZE_XY>().host_pointer()[0];
        const int prot_index_xyz =
            (*device_scratch).template get<buffer_data_type::PROT_SIZE_XYZ>().host_pointer()[0];
        const fp_type* prot_ds_grid_maps =
            (*device_scratch).template get<buffer_data_type::PROT_DS_GRID_MAPS>().dev_pointer();

        kernel = std::make_unique<ds_score_precomp_kernel<queue_type>>(batch_atoms,
                                                                      batch_ligands,
                                                                      scores_per_ligand,
                                                                      x_lig_p,
                                                                      y_lig_p,
                                                                      z_lig_p,
                                                                      lig_ds_types_p,
                                                                      lig_n_atm_p,
                                                                      prot_min,
                                                                      prot_max,
                                                                      prot_index_x,
                                                                      prot_index_xy,
                                                                      prot_index_xyz,
                                                                      prot_ds_grid_maps,
                                                                      scores_p,
                                                                      q);
      } else {
        const fp_type* prot_x_p =
            (*device_scratch).template get<buffer_data_type::PROT_X_COORDS>().dev_pointer();
        const fp_type* prot_y_p =
            (*device_scratch).template get<buffer_data_type::PROT_Y_COORDS>().dev_pointer();
        const fp_type* prot_z_p =
            (*device_scratch).template get<buffer_data_type::PROT_Z_COORDS>().dev_pointer();
        const ds_type* prot_ds_types_p =
            (*device_scratch).template get<buffer_data_type::PROT_DS_TYPES>().dev_pointer();
        const int prot_n_atm =
            (*device_scratch).template get<buffer_data_type::PROT_NUM_ATOMS>().host_pointer()[0];
        // drugscore scores pointer
        const fp_type* ds_scores = (*device_scratch).template get<buffer_data_type::DS_MODEL>().dev_pointer();

        kernel = std::make_unique<ds_score_kernel<queue_type>>(batch_atoms,
                                                              batch_ligands,
                                                              scores_per_ligand,
                                                              x_lig_p,
                                                              y_lig_p,
                                                              z_lig_p,
                                                              lig_ds_types_p,
                                                              lig_n_atm_p,
                                                              prot_x_p,
                                                              prot_y_p,
                                                              prot_z_p,
                                                              prot_ds_types_p,
                                                              prot_n_atm,
                                                              ds_scores,
                                                              scores_p,
                                                              q);
      }
    }

    void operator()() {
      assert(
          (((*this->scratch).template get<buffer_data_type::SCORES>().num_elements() % batch_ligands) == 0) &&
          "Number of scores is not a multiple of ligands in the batch");

      if constexpr (PRECOMPUTE) {
        const auto& ker_ptr = std::get<std::unique_ptr<ds_score_precomp_kernel<queue_type>>>(kernel);
        assert(ker_ptr && "Kernel method not yet prepared");
        (*ker_ptr)();
      } else {
        const auto& ker_ptr = std::get<std::unique_ptr<ds_score_kernel<queue_type>>>(kernel);
        assert(ker_ptr && "Kernel method not yet prepared");
        (*ker_ptr)();
      };
    }

    static batch_multiple get_batch_size(const int atoms,
                                         std::shared_ptr<queue_type> q,
                                         const knobs &conf,
                                         const size_t max_bucket_size) {
      (void) conf;
      const auto plain_multiple_info =
          normalize_batch_multiple(get_ds_score_batch_multiple<queue_type>(atoms, q));
      mudock::stage_bucket_trace("DS stage plain multiple for ",
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

    buffer_vector<ds_type, queue_type> lig_ds_types;

    std::shared_ptr<scratchpad<queue_type>> device_scratch;
    ligand_placement placement;
    point3D protein_center;

    std::variant<
      std::unique_ptr<ds_score_kernel<queue_type>>,
      std::unique_ptr<ds_score_precomp_kernel<queue_type>>
    > kernel;

    void teardown_impl(batch<static_molecule>& batch) override {
      assert(batch.num_ligands == batch_ligands && "Scoring algorithm received different batch for teardown");
      assert(batch_ligands ==
                 static_cast<int>((*this->scratch).template get<buffer_data_type::SCORES>().num_elements()) &&
             "Number of scores and ligands in batch are different");

      auto& scores_b = (*this->scratch).template get<buffer_data_type::SCORES>();
      scores_b.copy_device2host();
      for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
        auto& ligand = *batch.molecules[ligand_index];

        ligand.properties.assign(property_type::SCORE, std::to_string(scores_b()[ligand_index]));
      }
    };
  };

  template<typename queue_type>
  struct ds_score : ds_score_impl<queue_type, false> {
    using ds_score_impl<queue_type,false>::ds_score_impl;
  };

  template<typename queue_type>
  struct ds_score_precomp : ds_score_impl<queue_type, true> {
    using ds_score_impl<queue_type,true>::ds_score_impl;
  };

#endif
} // namespace mudock
