#pragma once

#include "mudock/type_alias.hpp"

#include <concepts>
#include <memory>
#include <mudock/chem/geom_ligand.hpp>
#include <mudock/compute/queue.hpp>
#include <mudock/compute/geom_transform_args.hpp>
#include <mudock/compute/buffer_utils.hpp>
#include <mudock/compute/ligand_placement.hpp>
#include <mudock/compute/scratchpad.hpp>
#include <mudock/compute/transform.hpp>
#include <mudock/cpp_implementation/centroid.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/cpp_implementation/mutate_cpp.hpp>
#include <mudock/grid/point3D.hpp>
#include <mudock/molecule.hpp>

namespace mudock {
  template<typename queue_t>
    requires std::derived_from<queue_t, queue>
  struct geometric: public transform<queue_t> {
    geometric(std::shared_ptr<scratchpad<queue_t>> _scratch,
              dynamic_molecule& protein,
              ligand_placement _placement = {})
        : transform<queue_t>(_scratch),
          ligand_fragments(_scratch->get_queue()),
          ligand_fragments_start(_scratch->get_queue()),
          frag_start_atom_indices(_scratch->get_queue()),
          frag_stop_atom_indices(_scratch->get_queue()),
          frag_indices_start(_scratch->get_queue()),
          protein_center(protein.get_center()),
          placement(_placement) {};

    void prepare(batch<static_molecule>& batch) override {
      batch_ligands                           = batch.num_ligands;
      batch_atoms                             = batch.batch_max_atoms;
      const int batch_rotamers                = batch_atoms - 3;
      const std::size_t tot_atoms_in_batch = static_cast<std::size_t>(batch_ligands) *
                                             static_cast<std::size_t>(batch_atoms);
      const std::size_t tot_rotamers_atoms_in_batch =
          tot_atoms_in_batch * static_cast<std::size_t>(batch_rotamers);
      const std::size_t tot_rotamers_in_batch =
          static_cast<std::size_t>(batch_ligands) * static_cast<std::size_t>(batch_rotamers);
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
        // Place the molecule according to the selected initial placement policy.
        apply_ligand_placement(ligand, placement, protein_center);

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

      args = {batch_ligands,
              batch_atoms,
              chromsomes_per_ligand,
              chromosomes_p,
              x_coords_p,
              y_coords_p,
              z_coords_p,
              num_atoms_p,
              num_rotamers_p,
              ligand_fragments_p,
              ligand_fragments_start_p,
              frag_start_indices_p,
              frag_stop_indices_p,
              frag_indices_start_p,
              x_scratch_p,
              y_scratch_p,
              z_scratch_p};
    }

    void operator()() override;

    [[nodiscard]] geom_transform_args kernel_args() const {
      return args;
    }

    void set_chromosomes_buffer(const chromosome* chromosomes) {
      args.chromosomes = chromosomes;
    }

    static std::size_t get_shared_ligand_mem(const int max_atoms, const knobs conf) {
      const int chromosomes_per_ligand = std::max(1, static_cast<int>(conf.population_number));
      return sizeof(int) + sizeof(int) + 3 * sizeof(fp_type) * max_atoms +
             3 * sizeof(fp_type) * max_atoms * chromosomes_per_ligand;
    }

    static std::size_t get_private_ligand_mem(const int max_atoms, const knobs) {
      const std::size_t batch_rotamers              = static_cast<std::size_t>(max_atoms - 3);
      const std::size_t tot_rotamers_atoms_in_batch = static_cast<std::size_t>(max_atoms) * batch_rotamers;
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

  private:
    buffer_vector<int, queue_t> ligand_fragments;
    buffer_vector<int, queue_t> ligand_fragments_start;
    buffer_vector<int, queue_t> frag_start_atom_indices;
    buffer_vector<int, queue_t> frag_stop_atom_indices;
    buffer_vector<int, queue_t> frag_indices_start;

    int batch_ligands;
    int batch_atoms;
    point<fp_type, 3> protein_center;
    ligand_placement placement;
    geom_transform_args args{};

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
} // namespace mudock
