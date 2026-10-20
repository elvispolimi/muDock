#pragma once

#include "mudock/type_alias.hpp"

#include <concepts>
#include <memory>
#include <mudock/chem/geom_ligand.hpp>
#include <mudock/compute/batch_multiple.hpp>
#include <mudock/compute/queue.hpp>
#include <mudock/log.hpp>
#if !defined(__CUDACC__) && !defined(__HIPCC__)
  #include <mudock/compute/buffer_utils.hpp>
  #include <mudock/compute/scratchpad.hpp>
  #include <mudock/compute/transform.hpp>
#endif
#include <mudock/cpp_implementation/center_of_mass.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/cpp_implementation/mutate_cpp.hpp>
#include <mudock/grid/point3D.hpp>
#include <mudock/molecule.hpp>

namespace mudock {
  template<typename queue_type>
  batch_multiple get_geom_transform_batch_multiple(const int, std::shared_ptr<queue_type>) {
    return {};
  }

  template<typename queue_type>
    requires std::derived_from<queue_type, queue>
  struct geom_kernel {
    static constexpr char geom_region_name[] = "geometric_trasformation";
    geom_kernel(const int batch_ligands_,
                const int batch_atoms_,
                const int chromsomes_per_ligand_,
                const chromosome* __restrict__ chromosomes_b_,
                const int* __restrict__ num_atoms_b_,
                const int* __restrict__ num_rotamers_b_,
                const fp_type* __restrict__ x_coords_b_,
                const fp_type* __restrict__ y_coords_b_,
                const fp_type* __restrict__ z_coords_b_,
                fp_type* __restrict__ x_scratch_b_,
                fp_type* __restrict__ y_scratch_b_,
                fp_type* __restrict__ z_scratch_b_,
                const int* __restrict__ ligand_fragments_b_,
                const int* __restrict__ ligand_fragments_start_b_,
                const int* __restrict__ frag_indices_start_b_,
                const int* __restrict__ frag_start_indices_b_,
                const int* __restrict__ frag_stop_indices_b_,
                std::shared_ptr<queue_type> q_)
        : batch_ligands(batch_ligands_),
          batch_atoms(batch_atoms_),
          chromsomes_per_ligand(chromsomes_per_ligand_),
          chromosomes_b(chromosomes_b_),
          num_atoms_b(num_atoms_b_),
          num_rotamers_b(num_rotamers_b_),
          x_coords_b(x_coords_b_),
          y_coords_b(y_coords_b_),
          z_coords_b(z_coords_b_),
          x_scratch_b(x_scratch_b_),
          y_scratch_b(y_scratch_b_),
          z_scratch_b(z_scratch_b_),
          ligand_fragments_b(ligand_fragments_b_),
          ligand_fragments_start_b(ligand_fragments_start_b_),
          frag_indices_start_b(frag_indices_start_b_),
          frag_start_indices_b(frag_start_indices_b_),
          frag_stop_indices_b(frag_stop_indices_b_),
          q(q_) {}

    void operator()();
    inline void set_chromosomes_buffer(const chromosome* chromosomes_b_) { chromosomes_b = chromosomes_b_; }

    geom_kernel(const geom_kernel&)            = default;
    geom_kernel(geom_kernel&&)                 = default;
    geom_kernel& operator=(const geom_kernel&) = delete;
    geom_kernel& operator=(geom_kernel&&)      = delete;

    ~geom_kernel() = default;

  private:
    const int batch_ligands;
    const int batch_atoms;
    const int chromsomes_per_ligand;
    const chromosome* __restrict__ chromosomes_b;
    const int* __restrict__ num_atoms_b;
    const int* __restrict__ num_rotamers_b;
    const fp_type* __restrict__ x_coords_b;
    const fp_type* __restrict__ y_coords_b;
    const fp_type* __restrict__ z_coords_b;
    fp_type* __restrict__ x_scratch_b;
    fp_type* __restrict__ y_scratch_b;
    fp_type* __restrict__ z_scratch_b;
    const int* __restrict__ ligand_fragments_b;
    const int* __restrict__ ligand_fragments_start_b;
    const int* __restrict__ frag_indices_start_b;
    const int* __restrict__ frag_start_indices_b;
    const int* __restrict__ frag_stop_indices_b;
    std::shared_ptr<queue_type> q;
  };

#if !defined(__CUDACC__) && !defined(__HIPCC__)
  template<typename queue_t, typename layout_set = all_scratch_layouts>
    requires std::derived_from<queue_t, queue>
  struct geometric: public transform<queue_t, layout_set> {
    geometric(std::shared_ptr<scratchpad<queue_t, layout_set>> _scratch, dynamic_molecule& protein)
        : transform<queue_t, layout_set>(_scratch),
          geometry_data(_scratch->get_queue()),
          protein_center(protein.get_center()) {};

    void prepare(batch<static_molecule>& batch) {
      batch_ligands                           = batch.num_ligands;
      batch_atoms                             = batch.batch_max_atoms;
      const int batch_rotamers                = batch_atoms - 3;
      const int tot_atoms_in_batch            = batch_ligands * batch_atoms;
      const int tot_rotamers_atoms_in_batch   = tot_atoms_in_batch * batch_rotamers;
      const std::size_t tot_rotamers_in_batch = batch_ligands * batch_rotamers;
      auto q                                  = (*this->scratch).get_queue();
      const auto chromsomes_per_ligand =
          std::max(1, static_cast<int>((*this->scratch).configuration.population_number));
      const auto expected_chromosomes =
          static_cast<std::size_t>(batch_ligands) * static_cast<std::size_t>(chromsomes_per_ligand);
      (void) expected_chromosomes;
      assert(this->scratch->genetic().num_elements() ==
                 expected_chromosomes &&
             "Chromosomes buffer not allocated before geom transform construction");

      geometry_data.alloc(tot_rotamers_atoms_in_batch,
              batch_ligands + 1,
              tot_rotamers_in_batch,
              tot_rotamers_in_batch,
              batch_ligands + 1);
      auto geometry_host_view = geometry_data.host_view();
      geometry_host_view.fragments_start().values()[0] = 0;
      geometry_host_view.fragments_indices_start().values()[0] = 0;

      // TODO check what happens if coordinates are already available
      for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
        auto& ligand = *batch.molecules[ligand_index];
        geom_ligand geom_lig{ligand};
        // Atoms and bonds
        const int num_atoms = ligand.num_atoms();
        // Place the molecule to the center of the target protein
        const auto x = ligand.x(), y = ligand.y(), z = ligand.z();

        const auto ligand_center_of_mass = compute_center_of_mass(x, y, z, num_atoms);
        const auto offset                = protein_center - ligand_center_of_mass;
        translate_molecule<cpu_vectorization::AUTO>(x, y, z, num_atoms, offset.x(), offset.y(), offset.z());

        const auto num_rotamers = ligand.num_rotamers();
        assert(batch_rotamers > num_rotamers);

        std::memcpy(geometry_host_view.fragments().values().data() +
            geometry_host_view.fragments_start().values()[ligand_index],
                    geom_lig.fragments_masks(),
                    num_atoms * num_rotamers * sizeof(int));
        geometry_host_view.fragments_start().values()[ligand_index + 1] =
          static_cast<int>(geometry_host_view.fragments_start().values()[ligand_index] +
                   (num_atoms * num_rotamers));
        std::memcpy(geometry_host_view.fragments_start_atoms().values().data() +
                geometry_host_view.fragments_indices_start().values()[ligand_index],
                    geom_lig.fragmets_starts(),
                    num_rotamers * sizeof(int));
        std::memcpy(geometry_host_view.fragments_stop_atoms().values().data() +
            geometry_host_view.fragments_indices_start().values()[ligand_index],
                    geom_lig.fragments_stops(),
                    num_rotamers * sizeof(int));
        geometry_host_view.fragments_indices_start().values()[ligand_index + 1] =
            static_cast<int>(geometry_host_view.fragments_indices_start().values()[ligand_index] + num_rotamers);
      }
      geometry_data.copy_host2device();

      load_num_rotamers<queue_t>(batch, this->scratch);
      load_num_atoms<queue_t>(batch, this->scratch);
      auto& x_scratch_b = this->scratch->scratch_coordinates();
      if (load_coords<queue_t>(batch, this->scratch)) {
        x_scratch_b.alloc(tot_atoms_in_batch * chromsomes_per_ligand);
        x_scratch_b.set_valid();
      }

      // Binding
        auto batch_view = this->scratch->batch_info().device_view();
        auto genetic_view = this->scratch->genetic().device_view();
        auto coordinates_view = this->scratch->coordinates().device_view();
        auto scratch_view = x_scratch_b.device_view();
        auto geometry_view = geometry_data.device_view();
        chromosome* chromosomes_p = genetic_view.chromosomes().data();
        int* num_atoms_p = batch_view.num_atoms().data();
        int* num_rotamers_p = batch_view.num_rotamers().data();
        const fp_type* x_coords_p = coordinates_view.x().data();
        const fp_type* y_coords_p = coordinates_view.y().data();
        const fp_type* z_coords_p = coordinates_view.z().data();
        fp_type* x_scratch_p = scratch_view.x().data();
        fp_type* y_scratch_p = scratch_view.y().data();
        fp_type* z_scratch_p = scratch_view.z().data();
        int* ligand_fragments_p = geometry_view.fragments().values().data();
        int* ligand_fragments_start_p = geometry_view.fragments_start().values().data();
        int* frag_indices_start_p = geometry_view.fragments_indices_start().values().data();
        int* frag_start_indices_p = geometry_view.fragments_start_atoms().values().data();
        int* frag_stop_indices_p = geometry_view.fragments_stop_atoms().values().data();

      kernel = std::make_unique<geom_kernel<queue_t>>(batch_ligands,
                                                      batch_atoms,
                                                      chromsomes_per_ligand,
                                                      chromosomes_p,
                                                      num_atoms_p,
                                                      num_rotamers_p,
                                                      x_coords_p,
                                                      y_coords_p,
                                                      z_coords_p,
                                                      x_scratch_p,
                                                      y_scratch_p,
                                                      z_scratch_p,
                                                      ligand_fragments_p,
                                                      ligand_fragments_start_p,
                                                      frag_indices_start_p,
                                                      frag_start_indices_p,
                                                      frag_stop_indices_p,
                                                      q);
    }

    void operator()() {
      assert(
          (this->scratch->genetic().num_elements() % batch_ligands ==
           0) &&
          "Number of chromosomes per ligand is not a multiple of the number of ligands expected to be docked");

      assert(kernel && "Kernel method not yet prepared");
      (*kernel)();
    };

    inline void set_chromosomes_buffer(const chromosome* chromosomes_b_) {
      if (kernel) {
        kernel->set_chromosomes_buffer(chromosomes_b_);
      }
    }

    static std::size_t get_shared_ligand_mem(const int max_atoms, const knobs conf) {
      const int chromosomes_per_ligand = std::max(1, static_cast<int>(conf.population_number));
      return sizeof(int) + sizeof(int) + 3 * sizeof(fp_type) * max_atoms +
             3 * sizeof(fp_type) * max_atoms * chromosomes_per_ligand;
    }

    static std::size_t get_private_ligand_mem(const int max_atoms, const knobs) {
      const int batch_rotamers              = max_atoms - 3;
      const int tot_rotamers_atoms_in_batch = max_atoms * batch_rotamers;
      std::size_t mem{0};
      mem += sizeof(int) * tot_rotamers_atoms_in_batch; // ligand fragments
      mem += sizeof(int);                               // ligand_fragments_start
      mem += sizeof(int) * batch_rotamers;              // frag_start_atom_indices
      mem += sizeof(int) * batch_rotamers;              // frag_stop_atom_indices
      mem += sizeof(int);                               // frag_indices_start
      return mem;
    }

    static std::size_t get_ligand_mem(const int max_atoms, const knobs conf) {
      return get_shared_ligand_mem(max_atoms, conf) + get_private_ligand_mem(max_atoms, conf);
    }

    static batch_multiple get_batch_size(const int atoms,
                                         std::shared_ptr<queue_t> q,
                                         const knobs& conf,
                                         const size_t max_bucket_size) {
      (void) conf;
      const auto plain_multiple_info =
          normalize_batch_multiple(get_geom_transform_batch_multiple<queue_t>(atoms, q));
      mudock::stage_bucket_trace("GEOM stage plain multiple for ",
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
    soa_buffer<GeometryData<>, queue_t> geometry_data;

    int batch_ligands;
    int batch_atoms;
    point<fp_type, 3> protein_center;
    std::unique_ptr<geom_kernel<queue_t>> kernel;

    void teardown_impl(batch<static_molecule>& batch) {
      const auto chromsomes_per_ligand =
          std::max(1, static_cast<int>((*this->scratch).configuration.population_number));
      auto& scratch_coordinates = this->scratch->scratch_coordinates();
      scratch_coordinates.copy_device2host();
      (*this->scratch).get_queue()->synchronize();
      auto scratch_view = scratch_coordinates.host_view();
      for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
        auto& ligand           = *batch.molecules[ligand_index];
        const int num_atoms    = ligand.num_atoms();
        const int stride_atoms = ligand_index * batch_atoms * chromsomes_per_ligand;
        std::memcpy(ligand.x(), scratch_view.x().data() + stride_atoms, num_atoms * sizeof(fp_type));
        std::memcpy(ligand.y(), scratch_view.y().data() + stride_atoms, num_atoms * sizeof(fp_type));
        std::memcpy(ligand.z(), scratch_view.z().data() + stride_atoms, num_atoms * sizeof(fp_type));
      }
    };
  };
#endif
} // namespace mudock
