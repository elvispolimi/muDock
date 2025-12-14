#pragma once

#include <concepts>
#include <cstring>
#include <functional>
#include <memory>
#include <mudock/batch.hpp>
#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/chem/autodock_ligand.hpp>
#include <mudock/chem/autodock_protein.hpp>
#include <mudock/compute/buffer.hpp>
#include <mudock/compute/scoring.hpp>
#include <mudock/compute/scratchpad.hpp>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {

  template<typename queue_type>
    requires std::derived_from<queue_type, queue>
  struct adt_score_kernel {
    adt_score_kernel(int scores_per_ligand_,
                     int batch_ligands_,
                     int batch_atoms_,
                     const int *__restrict__ num_atoms_b_,
                     const int *__restrict__ num_rotamers_b_,
                     const int *__restrict__ num_nonbonds_b_,
                     const fp_type *__restrict__ x_scratch_b_,
                     const fp_type *__restrict__ y_scratch_b_,
                     const fp_type *__restrict__ z_scratch_b_,
                     const fp_type *__restrict__ vols_b_,
                     const fp_type *__restrict__ solpars_b_,
                     const fp_type *__restrict__ charges_b_,
                     const int *__restrict__ map_offsets_b_,
                     const int *__restrict__ nonbond_a1_b_,
                     const int *__restrict__ nonbond_a2_b_,
                     const fp_type *__restrict__ nonbond_cA_b_,
                     const fp_type *__restrict__ nonbond_cB_b_,
                     const int *__restrict__ nonbond_xB_b_,
                     const fp_type *__restrict__ grid_maps_,
                     const fp_type *__restrict__ minimum_,
                     const fp_type *__restrict__ maximum_,
                     const fp_type *__restrict__ center_,
                     int map_index_x_,
                     int map_index_xy_,
                     int map_index_xyz_,
                     fp_type *__restrict__ scores_b_)
        : scores_per_ligand(scores_per_ligand_),
          batch_ligands(batch_ligands_),
          batch_atoms(batch_atoms_),
          num_atoms_b(num_atoms_b_),
          num_rotamers_b(num_rotamers_b_),
          num_nonbonds_b(num_nonbonds_b_),
          x_scratch_b(x_scratch_b_),
          y_scratch_b(y_scratch_b_),
          z_scratch_b(z_scratch_b_),
          vols_b(vols_b_),
          solpars_b(solpars_b_),
          charges_b(charges_b_),
          map_offsets_b(map_offsets_b_),
          nonbond_a1_b(nonbond_a1_b_),
          nonbond_a2_b(nonbond_a2_b_),
          nonbond_cA_b(nonbond_cA_b_),
          nonbond_cB_b(nonbond_cB_b_),
          nonbond_xB_b(nonbond_xB_b_),
          grid_maps(grid_maps_),
          minimum(minimum_),
          maximum(maximum_),
          center(center_),
          map_index_x(map_index_x_),
          map_index_xy(map_index_xy_),
          map_index_xyz(map_index_xyz_),
          scores_b(scores_b_) {}

    void operator()();

    adt_score_kernel(const adt_score_kernel &)            = default;
    adt_score_kernel(adt_score_kernel &&)                 = default;
    adt_score_kernel &operator=(const adt_score_kernel &) = delete;
    adt_score_kernel &operator=(adt_score_kernel &&)      = delete;

    ~adt_score_kernel() = default;

  private:
    int scores_per_ligand;
    int batch_ligands;
    int batch_atoms;
    const int *__restrict__ num_atoms_b;
    const int *__restrict__ num_rotamers_b;
    const int *__restrict__ num_nonbonds_b;
    const fp_type *__restrict__ x_scratch_b;
    const fp_type *__restrict__ y_scratch_b;
    const fp_type *__restrict__ z_scratch_b;
    const fp_type *__restrict__ vols_b;
    const fp_type *__restrict__ solpars_b;
    const fp_type *__restrict__ charges_b;
    const int *__restrict__ map_offsets_b;
    const int *__restrict__ nonbond_a1_b;
    const int *__restrict__ nonbond_a2_b;
    const fp_type *__restrict__ nonbond_cA_b;
    const fp_type *__restrict__ nonbond_cB_b;
    const int *__restrict__ nonbond_xB_b;
    const fp_type *__restrict__ grid_maps;
    const fp_type *__restrict__ minimum;
    const fp_type *__restrict__ maximum;
    const fp_type *__restrict__ center;
    int map_index_x;
    int map_index_xy;
    int map_index_xyz;
    fp_type *__restrict__ scores_b;
  };

  // TODO check that the object type and the kernel impl are the same
  template<typename queue_type>
  struct adt_score: public scoring<queue_type> {
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
        prot_index_x()[0]   = adt_prot.get_size_x();
        prot_index_xy()[0]  = adt_prot.get_size_xy();
        prot_index_xyz()[0] = adt_prot.get_size_xyz();
        std::memcpy(prot_grid_maps(),
                    adt_prot.get_maps_pointer(),
                    adt_prot.get_map_flat_size() * num_autodock_grids() * sizeof(fp_type));
      }
    }

    void prepare(batch<static_molecule> &batch) {
      batch_ligands                = batch.num_ligands;
      batch_atoms                  = batch.batch_max_atoms;
      const int batch_non_bonds    = batch_ligands * batch_atoms * batch_atoms;
      const int tot_atoms_in_batch = batch_ligands * batch_atoms;

      // TODO issue here if another batch arrives
      // Or any of the other buffers containing ligands info
      auto &scratch_x      = (*this->scratch).template get<buffer_data_type::X_SCRATCH>();
      auto &scratch_y      = (*this->scratch).template get<buffer_data_type::Y_SCRATCH>();
      auto &scratch_z      = (*this->scratch).template get<buffer_data_type::Z_SCRATCH>();
      auto &num_atoms_b    = (*this->scratch).template get<buffer_data_type::NUM_ATOMS>();
      auto &num_rotamers_b = (*this->scratch).template get<buffer_data_type::NUM_ROTAMERS>();
      auto &score_b        = (*this->scratch).template get<buffer_data_type::SCORES>();

      if (!scratch_x.is_valid()) {
        score_b.alloc(batch_ligands);
        scratch_x.alloc(tot_atoms_in_batch);
        scratch_y.alloc(tot_atoms_in_batch);
        scratch_z.alloc(tot_atoms_in_batch);
        num_atoms_b.alloc(batch_ligands);
        num_rotamers_b.alloc(batch_ligands);
        for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
          auto &ligand = *batch.molecules[ligand_index];

          const int stride_atoms         = ligand_index * batch_atoms;
          const int num_atoms            = ligand.num_atoms();
          num_atoms_b()[ligand_index]    = num_atoms;
          num_rotamers_b()[ligand_index] = ligand.num_rotamers();

          const auto x = ligand.x(), y = ligand.y(), z = ligand.z();
          std::memcpy((void *) (scratch_x() + stride_atoms), x, num_atoms * sizeof(fp_type));
          std::memcpy((void *) (scratch_y() + stride_atoms), y, num_atoms * sizeof(fp_type));
          std::memcpy((void *) (scratch_z() + stride_atoms), z, num_atoms * sizeof(fp_type));
        }
        scratch_x.copy_host2device();
        scratch_y.copy_host2device();
        scratch_z.copy_host2device();
        num_atoms_b.copy_host2device();
        num_rotamers_b.copy_host2device();
        score_b.is_valid();
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
        num_nonbond()[ligand_index + 1] = num_nonbond()[ligand_index] + non_bond_size;

        // Autodock typing
        std::memcpy((void *) (vols() + stride_atoms), adt_ligand.vol(), num_atoms * sizeof(fp_type));
        std::memcpy((void *) (solpars() + stride_atoms), adt_ligand.solpar(), num_atoms * sizeof(fp_type));
        std::memcpy((void *) (charges() + stride_atoms), ligand.charge(), num_atoms * sizeof(fp_type));

        std::memcpy((void *) (map_offsets() + stride_atoms),
                    adt_ligand.atom_map_index(),
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
      const auto scores_per_ligand =
          (*this->scratch).template get<buffer_data_type::SCORES>().num_elements() / batch_ligands;
      assert((*this->scratch).template get<buffer_data_type::X_SCRATCH>().num_elements() !=
                 (scores_per_ligand * batch_ligands) &&
             "Number of scores per ligands does not match the allocated coordinates space");
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

      const fp_type *grid_maps =
          (*device_scratch).template get<buffer_data_type::PROT_GRID_MAPS>().dev_pointer();
      const fp_type *minimum = (*device_scratch).template get<buffer_data_type::PROT_MIN>().dev_pointer();
      const fp_type *maximum = (*device_scratch).template get<buffer_data_type::PROT_MAX>().dev_pointer();
      const fp_type *center  = (*device_scratch).template get<buffer_data_type::PROT_CENTER>().dev_pointer();
      const int map_index_x =
          (*device_scratch).template get<buffer_data_type::PROT_SIZE_X>().host_pointer()[0];
      const int map_index_xy =
          (*device_scratch).template get<buffer_data_type::PROT_SIZE_XY>().host_pointer()[0];
      const int map_index_xyz =
          (*device_scratch).template get<buffer_data_type::PROT_SIZE_XYZ>().host_pointer()[0];

      fp_type *scores_b = (*this->scratch).template get<buffer_data_type::SCORES>().dev_pointer();

      kernel = std::make_unique<adt_score_kernel<queue_type>>(scores_per_ligand,
                                                              batch_ligands,
                                                              batch_atoms,
                                                              num_atoms_b.dev_pointer(),
                                                              num_rotamers_b.dev_pointer(),
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
                                                              scores_b);
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
      assert(batch_ligands ==
                 static_cast<int>((*this->scratch).template get<buffer_data_type::SCORES>().num_elements()) &&
             "Number of scores and ligands in batch are different");

      auto &scores_b = (*this->scratch).template get<buffer_data_type::SCORES>();
      scores_b.copy_device2host();
      for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index) {
        auto &ligand = *batch.molecules[ligand_index];

        ligand.properties.assign(property_type::SCORE, std::to_string(scores_b()[ligand_index]));
      }
    };
  };
} // namespace mudock
