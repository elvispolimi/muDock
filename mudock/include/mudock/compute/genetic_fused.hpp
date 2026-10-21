#pragma once

#include <algorithm>
#include <chrono>
#include <concepts>
#include <cstring>
#include <memory>
#include <mudock/chem/autodock_ligand.hpp>
#include <mudock/chem/autodock_protein.hpp>
#include <mudock/chem/geom_ligand.hpp>
#include <mudock/compute/docking.hpp>
#include <mudock/compute/bucket_size.hpp>
#include <mudock/cpp_implementation/center_of_mass.hpp>
#include <mudock/cuda_implementation/genetic_fused_cuda.cuh>
#include <mudock/log.hpp>
#include <mudock/utils.hpp>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#ifdef MUDOCK_USE_CUDA
  #include <mudock/cuda_implementation/queue_cuda.cuh>
#endif

namespace mudock {
  template<typename queue_t>
  struct genetic_fused: docking<queue_t> {
    static constexpr const char stage_name[] = "GENETIC_FUSED";

    genetic_fused(std::shared_ptr<scratchpad<queue_t>> scratch,
                  std::shared_ptr<scratchpad<queue_t>> device_scratch_,
                  dynamic_molecule& protein)
        : docking<queue_t>(scratch),
          device_scratch(std::move(device_scratch_)),
          protein_center(protein.get_center()),
          ligand_device(std::make_unique<object<std::byte>>(scratch->get_queue())) {}

    static void prepare_shared_protein(std::shared_ptr<scratchpad<queue_t>> device_scratch,
                                       dynamic_molecule& protein) {
      if (device_scratch->template exists<buffer_data_type::PROT_GRID_MAPS>()) return;

      autodock_protein adt_protein(protein);
      auto& prot_min = device_scratch->template get<buffer_data_type::PROT_MIN>();
      auto& prot_max = device_scratch->template get<buffer_data_type::PROT_MAX>();
      auto& prot_center = device_scratch->template get<buffer_data_type::PROT_CENTER>();
      auto& prot_index_x = device_scratch->template get<buffer_data_type::PROT_SIZE_X>();
      auto& prot_index_xy = device_scratch->template get<buffer_data_type::PROT_SIZE_XY>();
      auto& prot_index_xyz = device_scratch->template get<buffer_data_type::PROT_SIZE_XYZ>();
      auto& prot_grid_maps = device_scratch->template get<buffer_data_type::PROT_GRID_MAPS>();

      prot_min.alloc(3);
      prot_max.alloc(3);
      prot_center.alloc(3);
      prot_index_x.alloc(1);
      prot_index_xy.alloc(1);
      prot_index_xyz.alloc(1);
      prot_grid_maps.alloc(adt_protein.get_size_xyz() * num_autodock_grids());

      std::memcpy(prot_min(), adt_protein.get_min_p(), 3 * sizeof(fp_type));
      std::memcpy(prot_max(), adt_protein.get_max_p(), 3 * sizeof(fp_type));
      std::memcpy(prot_center(), adt_protein.get_center_p(), 3 * sizeof(fp_type));
      prot_index_x()[0] = static_cast<int>(adt_protein.get_size_x());
      prot_index_xy()[0] = static_cast<int>(adt_protein.get_size_xy());
      prot_index_xyz()[0] = static_cast<int>(adt_protein.get_size_xyz());
      std::memcpy(prot_grid_maps(),
                  adt_protein.get_maps_pointer(),
                  adt_protein.get_size_xyz() * num_autodock_grids() * sizeof(fp_type));

      prot_min.copy_host2device();
      prot_max.copy_host2device();
      prot_center.copy_host2device();
      prot_index_x.copy_host2device();
      prot_index_xy.copy_host2device();
      prot_index_xyz.copy_host2device();
      prot_grid_maps.copy_host2device();
      device_scratch->get_queue()->synchronize();
    }

    void prepare(batch<static_molecule>& batch) override {
#if defined(MUDOCK_USE_CUDA)
      if constexpr (!std::same_as<queue_t, queue_cuda>)
        throw std::runtime_error("Fused search is supported only with CUDA:GPU");
#else
      throw std::runtime_error("Fused search is unavailable because this build has no CUDA backend");
#endif
      auto q = this->get_queue();
      const auto& conf = this->scratch->configuration;
      const int ligands = batch.num_ligands;
      const int atoms = batch.batch_max_atoms;
      const int population = static_cast<int>(conf.population_number);
      const auto ligand_count = static_cast<std::size_t>(ligands);
      const auto atom_stride = static_cast<std::size_t>(atoms);
      const auto atom_capacity = ligand_count * atom_stride;
      const auto max_rotamers = static_cast<std::size_t>(std::max(0, atoms - 3));
      const auto fragment_capacity = atom_capacity * max_rotamers;
      const auto rotamer_capacity = ligand_count * max_rotamers;
      const auto nonbond_capacity = atom_capacity * atom_stride;
      const auto population_capacity = ligand_count * static_cast<std::size_t>(population);

      if (!device_scratch->template get<buffer_data_type::PROT_GRID_MAPS>().is_valid())
        throw std::runtime_error("Fused protein data was not initialized before worker startup");
      const auto map_flat_size = device_scratch->template get<buffer_data_type::PROT_SIZE_XYZ>().host_pointer()[0];

      const auto input_bytes = ligand_input_bytes(ligands, atoms);
      const auto total_bytes = ligand_slab_bytes(ligands, atoms, conf);
      if (host_ligand_scratch.size() < input_bytes) host_ligand_scratch.resize(input_bytes);

      auto* host_cursor = host_ligand_scratch.data();
      const auto take_host = [&]<typename T>(const std::size_t count) {
        auto* field = reinterpret_cast<T*>(host_cursor);
        host_cursor += count * sizeof(T);
        return field;
      };
      auto* num_atoms = take_host.template operator()<int>(ligand_count);
      auto* num_rotamers = take_host.template operator()<int>(ligand_count);
      auto* x_coords = take_host.template operator()<fp_type>(atom_capacity);
      auto* y_coords = take_host.template operator()<fp_type>(atom_capacity);
      auto* z_coords = take_host.template operator()<fp_type>(atom_capacity);
      auto* fragments = take_host.template operator()<int>(fragment_capacity);
      auto* fragments_start = take_host.template operator()<int>(ligand_count + 1);
      auto* frag_indices_start = take_host.template operator()<int>(ligand_count + 1);
      auto* frag_start = take_host.template operator()<int>(rotamer_capacity);
      auto* frag_stop = take_host.template operator()<int>(rotamer_capacity);
      auto* num_nonbond = take_host.template operator()<int>(ligand_count + 1);
      auto* vols = take_host.template operator()<fp_type>(atom_capacity);
      auto* solpars = take_host.template operator()<fp_type>(atom_capacity);
      auto* charges = take_host.template operator()<fp_type>(atom_capacity);
      auto* map_offsets = take_host.template operator()<int>(atom_capacity);
      auto* nonbond_a1 = take_host.template operator()<int>(nonbond_capacity);
      auto* nonbond_a2 = take_host.template operator()<int>(nonbond_capacity);
      auto* nonbond_cA = take_host.template operator()<fp_type>(nonbond_capacity);
      auto* nonbond_cB = take_host.template operator()<fp_type>(nonbond_capacity);
      auto* nonbond_xB = take_host.template operator()<int>(nonbond_capacity);

      std::size_t fragment_cursor = 0;
      std::size_t rotamer_cursor = 0;
      std::size_t nonbond_cursor = 0;
      fragments_start[0] = frag_indices_start[0] = num_nonbond[0] = 0;
      for (int ligand_index = 0; ligand_index < ligands; ++ligand_index) {
        auto& ligand = *batch.molecules[ligand_index];
        const int ligand_atoms = ligand.num_atoms();
        const int ligand_rotamers = static_cast<int>(ligand.num_rotamers());
        const auto atom_offset = static_cast<std::size_t>(ligand_index) * atom_stride;
        num_atoms[ligand_index] = ligand_atoms;
        num_rotamers[ligand_index] = ligand_rotamers;

        const auto centroid = compute_center_of_mass(ligand.x(), ligand.y(), ligand.z(), ligand_atoms);
        const auto offset = protein_center - centroid;
        for (int atom = 0; atom < ligand_atoms; ++atom) {
          ligand.x()[atom] += offset.x();
          ligand.y()[atom] += offset.y();
          ligand.z()[atom] += offset.z();
        }
        std::memcpy(x_coords + atom_offset, ligand.x(), ligand_atoms * sizeof(fp_type));
        std::memcpy(y_coords + atom_offset, ligand.y(), ligand_atoms * sizeof(fp_type));
        std::memcpy(z_coords + atom_offset, ligand.z(), ligand_atoms * sizeof(fp_type));

        geom_ligand geometry(ligand);
        const auto fragment_count = static_cast<std::size_t>(ligand_atoms) * ligand_rotamers;
        std::memcpy(fragments + fragment_cursor,
                    geometry.fragments_masks(),
                    fragment_count * sizeof(int));
        std::memcpy(frag_start + rotamer_cursor,
                    geometry.fragmets_starts(),
                    ligand_rotamers * sizeof(int));
        std::memcpy(frag_stop + rotamer_cursor,
                    geometry.fragments_stops(),
                    ligand_rotamers * sizeof(int));
        fragment_cursor += fragment_count;
        rotamer_cursor += static_cast<std::size_t>(ligand_rotamers);
        fragments_start[ligand_index + 1] = static_cast<int>(fragment_cursor);
        frag_indices_start[ligand_index + 1] = static_cast<int>(rotamer_cursor);

        autodock_ligand adt_ligand{ligand};
        adt_ligand.update_offsets(map_flat_size);
        std::memcpy(vols + atom_offset, adt_ligand.vol(), ligand_atoms * sizeof(fp_type));
        std::memcpy(solpars + atom_offset, adt_ligand.solpar(), ligand_atoms * sizeof(fp_type));
        std::memcpy(charges + atom_offset, ligand.charge(), ligand_atoms * sizeof(fp_type));
        std::memcpy(map_offsets + atom_offset,
                    adt_ligand.atom_map_offsets(),
                    ligand_atoms * sizeof(int));

        const auto nonbond_count = adt_ligand.non_bond_size();
        std::memcpy(nonbond_a1 + nonbond_cursor, adt_ligand.non_bond_A(), nonbond_count * sizeof(int));
        std::memcpy(nonbond_a2 + nonbond_cursor, adt_ligand.non_bond_B(), nonbond_count * sizeof(int));
        std::memcpy(nonbond_cA + nonbond_cursor,
                    adt_ligand.non_bond_cA(),
                    nonbond_count * sizeof(fp_type));
        std::memcpy(nonbond_cB + nonbond_cursor,
                    adt_ligand.non_bond_cB(),
                    nonbond_count * sizeof(fp_type));
        std::memcpy(nonbond_xB + nonbond_cursor, adt_ligand.non_bond_xB(), nonbond_count * sizeof(int));
        nonbond_cursor += nonbond_count;
        num_nonbond[ligand_index + 1] = static_cast<int>(nonbond_cursor);
      }

      ligand_device->alloc(total_bytes);
      q->copy_host2device(host_ligand_scratch.data(), ligand_device->dev_pointer(), input_bytes);
      auto* device_cursor = reinterpret_cast<std::byte*>(ligand_device->dev_pointer());
      const auto take_device = [&]<typename T>(const std::size_t count) {
        auto* field = reinterpret_cast<T*>(device_cursor);
        device_cursor += count * sizeof(T);
        return field;
      };
      const auto take_device_input = [&]<typename T>(const std::size_t count) {
        auto* field = take_device.template operator()<T>(count);
        return static_cast<const T*>(field);
      };
      const int* num_atoms_device = take_device_input.template operator()<int>(ligand_count);
      const int* num_rotamers_device = take_device_input.template operator()<int>(ligand_count);
      const fp_type* x_coords_device = take_device_input.template operator()<fp_type>(atom_capacity);
      const fp_type* y_coords_device = take_device_input.template operator()<fp_type>(atom_capacity);
      const fp_type* z_coords_device = take_device_input.template operator()<fp_type>(atom_capacity);
      const int* fragments_device = take_device_input.template operator()<int>(fragment_capacity);
      const int* fragments_start_device = take_device_input.template operator()<int>(ligand_count + 1);
      const int* frag_indices_start_device = take_device_input.template operator()<int>(ligand_count + 1);
      const int* frag_start_device = take_device_input.template operator()<int>(rotamer_capacity);
      const int* frag_stop_device = take_device_input.template operator()<int>(rotamer_capacity);
      const int* num_nonbond_device = take_device_input.template operator()<int>(ligand_count + 1);
      const fp_type* vols_device = take_device_input.template operator()<fp_type>(atom_capacity);
      const fp_type* solpars_device = take_device_input.template operator()<fp_type>(atom_capacity);
      const fp_type* charges_device = take_device_input.template operator()<fp_type>(atom_capacity);
      const int* map_offsets_device = take_device_input.template operator()<int>(atom_capacity);
      const int* nonbond_a1_device = take_device_input.template operator()<int>(nonbond_capacity);
      const int* nonbond_a2_device = take_device_input.template operator()<int>(nonbond_capacity);
      const fp_type* nonbond_cA_device = take_device_input.template operator()<fp_type>(nonbond_capacity);
      const fp_type* nonbond_cB_device = take_device_input.template operator()<fp_type>(nonbond_capacity);
      const int* nonbond_xB_device = take_device_input.template operator()<int>(nonbond_capacity);

      auto* current_device = take_device.template operator()<chromosome>(population_capacity);
      auto* next_device = take_device.template operator()<chromosome>(population_capacity);
      auto* scores_device = take_device.template operator()<fp_type>(population_capacity);
      auto* best_scores_device = take_device.template operator()<fp_type>(ligand_count);
      auto* best_chromosomes_device = take_device.template operator()<chromosome>(ligand_count);
      auto* x_scratch_device = take_device.template operator()<fp_type>(atom_capacity);
      auto* y_scratch_device = take_device.template operator()<fp_type>(atom_capacity);
      auto* z_scratch_device = take_device.template operator()<fp_type>(atom_capacity);

      const auto seed = conf.seed.value_or(static_cast<std::size_t>(
          std::chrono::high_resolution_clock::now().time_since_epoch().count()));
      genetic_fused_data data{ligands,
                              atoms,
                              population,
                              static_cast<int>(conf.num_generations),
                              static_cast<int>(conf.tournament_length),
                              conf.mutation_prob,
                              seed,
                              current_device,
                              next_device,
                              scores_device,
                              best_scores_device,
                              best_chromosomes_device,
                              num_atoms_device,
                              num_rotamers_device,
                              x_coords_device,
                              y_coords_device,
                              z_coords_device,
                              x_scratch_device,
                              y_scratch_device,
                              z_scratch_device,
                              fragments_device,
                              fragments_start_device,
                              frag_indices_start_device,
                              frag_start_device,
                              frag_stop_device,
                              num_nonbond_device,
                              vols_device,
                              solpars_device,
                              charges_device,
                              map_offsets_device,
                              nonbond_a1_device,
                              nonbond_a2_device,
                              nonbond_cA_device,
                              nonbond_cB_device,
                              nonbond_xB_device,
                              device_scratch->template get<buffer_data_type::PROT_GRID_MAPS>().dev_pointer(),
                              device_scratch->template get<buffer_data_type::PROT_MIN>().dev_pointer(),
                              device_scratch->template get<buffer_data_type::PROT_MAX>().dev_pointer(),
                              device_scratch->template get<buffer_data_type::PROT_CENTER>().dev_pointer(),
                              device_scratch->template get<buffer_data_type::PROT_SIZE_X>().host_pointer()[0],
                              device_scratch->template get<buffer_data_type::PROT_SIZE_XY>().host_pointer()[0],
                              map_flat_size};
      fused_kernel = std::make_unique<genetic_fused_launcher<queue_t>>(data, q);
      best_scores_device_ptr = best_scores_device;
      best_chromosomes_device_ptr = best_chromosomes_device;
      batch_ligands = ligands;
      best_scores_result.resize(ligand_count);
      best_chromosomes_result.resize(ligand_count);
      mudock::stage_bucket_trace(stage_name,
                                 " prepared device=",
                                 q->get_id(),
                                 " ligands=",
                                 ligands,
                                 " atoms=",
                                 atoms,
                                 " slab_bytes=",
                                 total_bytes,
                                 " h2d_bytes=",
                                 input_bytes,
                                 " population=",
                                 population,
                                 " generations=",
                                 conf.num_generations,
                                 " execution=one CUDA kernel launch per batch");
    }

    void operator()() override {
      if (!fused_kernel) throw std::runtime_error("Fused stage was not prepared");
      (*fused_kernel)();
    }

    [[nodiscard]] const std::vector<fp_type>& best_scores() const { return best_scores_result; }
    [[nodiscard]] const std::vector<chromosome>& best_chromosomes() const { return best_chromosomes_result; }

    static std::size_t get_ligand_mem(const int atoms, const knobs& conf) {
      return ligand_slab_bytes(1, atoms, conf);
    }

    template<typename queue_type>
    static int get_batch_size(const int atoms,
                              std::shared_ptr<queue_type> q,
                              const knobs& conf,
                              const std::size_t max_mem) {
      const auto mem_per_ligand = get_ligand_mem(atoms, conf);
      const auto max_bucket_size = std::max<std::size_t>(1, max_mem / mem_per_ligand);
      const auto fallback = [&]() -> batch_multiple {
#if defined(MUDOCK_USE_CUDA)
        if constexpr (std::same_as<queue_type, queue_cuda>)
          return get_genetic_fused_batch_multiple_cuda(q->get_id(), atoms);
#endif
        return {};
      };
      mudock::stage_bucket_trace("PIPELINE(",
                                 stage_name,
                                 ") pre-resolve for ",
                                 atoms,
                                 " atoms: mem_budget=",
                                 max_mem,
                                 " B, mem_per_ligand=",
                                 mem_per_ligand,
                                 " B, max_bucket_size=",
                                 max_bucket_size);
      return resolve_stage_bucket_size(stage_name,
                                       atoms,
                                       max_bucket_size,
                                       mem_per_ligand,
                                       q->honors_stage_bucket_policy(),
                                       fallback);
    }

  protected:
    void teardown_impl(batch<static_molecule>& batch) override {
      if (batch.num_ligands != batch_ligands)
        throw std::runtime_error("Fused stage received a different batch at teardown");
      auto q = this->get_queue();
      q->copy_device2host(best_scores_device_ptr,
                          best_scores_result.data(),
                          best_scores_result.size() * sizeof(fp_type));
      q->copy_device2host(best_chromosomes_device_ptr,
                          best_chromosomes_result.data(),
                          best_chromosomes_result.size() * sizeof(chromosome));
      q->synchronize();
      for (int ligand = 0; ligand < batch_ligands; ++ligand)
        batch.molecules[ligand]->properties.assign(property_type::SCORE,
                                                    std::to_string(best_scores_result[static_cast<std::size_t>(ligand)]));
    }

  private:
    static std::size_t ligand_input_bytes(const std::size_t ligands, const int atoms) {
      const auto atom_capacity = ligands * static_cast<std::size_t>(atoms);
      const auto max_rotamers = static_cast<std::size_t>(std::max(0, atoms - 3));
      const auto fragment_capacity = atom_capacity * max_rotamers;
      const auto rotamer_capacity = ligands * max_rotamers;
      const auto nonbond_capacity = atom_capacity * static_cast<std::size_t>(atoms);
      return 2 * ligands * sizeof(int) + 3 * atom_capacity * sizeof(fp_type) +
             fragment_capacity * sizeof(int) + 3 * (ligands + 1) * sizeof(int) +
             2 * rotamer_capacity * sizeof(int) +
             atom_capacity * (3 * sizeof(fp_type) + sizeof(int)) +
             nonbond_capacity * (3 * sizeof(int) + 2 * sizeof(fp_type));
    }

    static std::size_t ligand_slab_bytes(const std::size_t ligands, const int atoms, const knobs& conf) {
      const auto atom_capacity = ligands * static_cast<std::size_t>(atoms);
      const auto population_capacity = ligands * std::max<std::size_t>(1, conf.population_number);
      return ligand_input_bytes(ligands, atoms) +
             2 * population_capacity * sizeof(chromosome) +
             population_capacity * sizeof(fp_type) +
             ligands * (sizeof(fp_type) + sizeof(chromosome)) +
             3 * atom_capacity * sizeof(fp_type);
    }

    std::shared_ptr<scratchpad<queue_t>> device_scratch;
    point3D protein_center;
    std::unique_ptr<object<std::byte>> ligand_device;
    std::vector<std::byte> host_ligand_scratch;
    std::unique_ptr<genetic_fused_launcher<queue_t>> fused_kernel;
    std::vector<fp_type> best_scores_result;
    std::vector<chromosome> best_chromosomes_result;
    fp_type* best_scores_device_ptr{};
    chromosome* best_chromosomes_device_ptr{};
    int batch_ligands{};
  };
} // namespace mudock
