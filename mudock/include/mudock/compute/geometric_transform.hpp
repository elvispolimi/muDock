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
  template<typename queue_t>
    requires std::derived_from<queue_t, queue>
  struct geometric: public transform<queue_t> {
    geometric(std::shared_ptr<scratchpad<queue_t>> _scratch, dynamic_molecule& protein)
        : transform<queue_t>(_scratch),
          ligand_fragments(_scratch->get_queue()),
          ligand_fragments_start(_scratch->get_queue()),
          frag_start_atom_indices(_scratch->get_queue()),
          frag_stop_atom_indices(_scratch->get_queue()),
          frag_indices_start(_scratch->get_queue()),
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
      assert((*this->scratch).template get<buffer_data_type::CHROMOSOMES>().num_elements() ==
                 expected_chromosomes &&
             "Chromosomes buffer not allocated before geom transform construction");

      ligand_fragments.alloc(tot_rotamers_atoms_in_batch);
      ligand_fragments_start.alloc(batch_ligands + 1);
      ligand_fragments_start()[0] = 0;
      frag_start_atom_indices.alloc(tot_rotamers_in_batch);
      frag_stop_atom_indices.alloc(tot_rotamers_in_batch);
      frag_indices_start.alloc(batch_ligands + 1);
      frag_indices_start()[0] = 0;

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

        std::memcpy((void*) (ligand_fragments() + ligand_fragments_start()[ligand_index]),
                    geom_lig.fragments_masks(),
                    num_atoms * num_rotamers * sizeof(int));
        ligand_fragments_start()[ligand_index + 1] =
            static_cast<int>(ligand_fragments_start()[ligand_index] + (num_atoms * num_rotamers));
        std::memcpy((void*) (frag_start_atom_indices() + frag_indices_start()[ligand_index]),
                    geom_lig.fragmets_starts(),
                    num_rotamers * sizeof(int));
        std::memcpy((void*) (frag_stop_atom_indices() + frag_indices_start()[ligand_index]),
                    geom_lig.fragments_stops(),
                    num_rotamers * sizeof(int));
        frag_indices_start()[ligand_index + 1] =
            static_cast<int>(frag_indices_start()[ligand_index] + num_rotamers);
      }
      ligand_fragments.copy_host2device();
      ligand_fragments_start.copy_host2device();
      frag_start_atom_indices.copy_host2device();
      frag_stop_atom_indices.copy_host2device();
      frag_indices_start.copy_host2device();

      load_num_rotamers<queue_t>(batch, this->scratch);
      load_num_atoms<queue_t>(batch, this->scratch);
      auto& x_scratch_b = (*this->scratch).template get<buffer_data_type::X_SCRATCH>();
      auto& y_scratch_b = (*this->scratch).template get<buffer_data_type::Y_SCRATCH>();
      auto& z_scratch_b = (*this->scratch).template get<buffer_data_type::Z_SCRATCH>();
      if (load_coords<queue_t>(batch, this->scratch)) {
        x_scratch_b.alloc(tot_atoms_in_batch * chromsomes_per_ligand);
        y_scratch_b.alloc(tot_atoms_in_batch * chromsomes_per_ligand);
        z_scratch_b.alloc(tot_atoms_in_batch * chromsomes_per_ligand);
        x_scratch_b.set_valid();
        y_scratch_b.set_valid();
        z_scratch_b.set_valid();
      }

      // Binding
      chromosome* chromosomes_p =
          (*this->scratch).template get<buffer_data_type::CHROMOSOMES>().dev_pointer();
      int* num_atoms_p    = (*this->scratch).template get<buffer_data_type::NUM_ATOMS>().dev_pointer();
      int* num_rotamers_p = (*this->scratch).template get<buffer_data_type::NUM_ROTAMERS>().dev_pointer();
      const fp_type* x_coords_p = (*this->scratch).template get<buffer_data_type::X_COORDS>().dev_pointer();
      const fp_type* y_coords_p = (*this->scratch).template get<buffer_data_type::Y_COORDS>().dev_pointer();
      const fp_type* z_coords_p = (*this->scratch).template get<buffer_data_type::Z_COORDS>().dev_pointer();

      fp_type* x_scratch_p = x_scratch_b.dev_pointer();
      fp_type* y_scratch_p = y_scratch_b.dev_pointer();
      fp_type* z_scratch_p = z_scratch_b.dev_pointer();

      int* ligand_fragments_p       = ligand_fragments.dev_pointer();
      int* ligand_fragments_start_p = ligand_fragments_start.dev_pointer();
      int* frag_indices_start_p     = frag_indices_start.dev_pointer();
      int* frag_start_indices_p     = frag_start_atom_indices.dev_pointer();
      int* frag_stop_indices_p      = frag_stop_atom_indices.dev_pointer();

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
          ((*this->scratch).template get<buffer_data_type::CHROMOSOMES>().num_elements() % batch_ligands ==
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

    static int get_ligand_mem(const int max_atoms, const knobs conf) {
      return static_cast<int>(get_shared_ligand_mem(max_atoms, conf) +
                              get_private_ligand_mem(max_atoms, conf));
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
    buffer_vector<int, queue_t> ligand_fragments;
    buffer_vector<int, queue_t> ligand_fragments_start;
    buffer_vector<int, queue_t> frag_start_atom_indices;
    buffer_vector<int, queue_t> frag_stop_atom_indices;
    buffer_vector<int, queue_t> frag_indices_start;

    int batch_ligands;
    int batch_atoms;
    point<fp_type, 3> protein_center;
    std::unique_ptr<geom_kernel<queue_t>> kernel;

    void teardown_impl(batch<static_molecule>& batch) {
      const auto chromsomes_per_ligand =
          std::max(1, static_cast<int>((*this->scratch).configuration.population_number));
      auto& x_scratch_b = (*this->scratch).template get<buffer_data_type::X_SCRATCH>();
      auto& y_scratch_b = (*this->scratch).template get<buffer_data_type::Y_SCRATCH>();
      auto& z_scratch_b = (*this->scratch).template get<buffer_data_type::Z_SCRATCH>();
      x_scratch_b.copy_device2host();
      y_scratch_b.copy_device2host();
      z_scratch_b.copy_device2host();
      (*this->scratch).get_queue()->synchronize();
      for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
        auto& ligand           = *batch.molecules[ligand_index];
        const int num_atoms    = ligand.num_atoms();
        const int stride_atoms = ligand_index * batch_atoms * chromsomes_per_ligand;
        std::memcpy(ligand.x(), x_scratch_b.host_pointer() + stride_atoms, num_atoms * sizeof(fp_type));
        std::memcpy(ligand.y(), y_scratch_b.host_pointer() + stride_atoms, num_atoms * sizeof(fp_type));
        std::memcpy(ligand.z(), z_scratch_b.host_pointer() + stride_atoms, num_atoms * sizeof(fp_type));
      }
    };
  };
#endif
} // namespace mudock
