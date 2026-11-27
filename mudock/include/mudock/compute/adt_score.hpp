#pragma once

#include <cstddef>
#include <cstring>
#include <mudock/batch.hpp>
#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/chem/autodock_ligand.hpp>
#include <mudock/chem/autodock_protein.hpp>
#include <mudock/compute/algorithm.hpp>
#include <mudock/compute/adt_score_args.hpp>
#include <mudock/compute/ligand_placement.hpp>
#include <mudock/compute/buffer_utils.hpp>
#include <mudock/compute/scoring.hpp>
#include <mudock/compute/scratchpad.hpp>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {

  struct adt_score_device_kernel;

  // TODO check that the object type and the kernel impl are the same
  template<typename queue_type>
  struct adt_score: public scoring<queue_type> {
    static constexpr const char stage_name[] = "ADT";
    using device_kernel = adt_score_device_kernel;

    adt_score(std::shared_ptr<scratchpad<queue_type>> _scratch,
              std::shared_ptr<scratchpad<queue_type>> _device_scratch,
              dynamic_molecule &protein,
              ligand_placement _placement = {})
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
          device_scratch(_device_scratch),
          placement(_placement),
          protein_center(protein.get_center()) {
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
        prot_grid_maps.copy_host2device();
        // Workers use this device-shared pointer from their own queues.
        device_scratch->get_queue()->synchronize();
      }
    }

    void prepare(batch<static_molecule>& batch) override {
      batch_ligands                = batch.num_ligands;
      batch_atoms                  = batch.batch_max_atoms;
      const int batch_non_bonds    = batch_ligands * batch_atoms * batch_atoms;
      const int tot_atoms_in_batch = batch_ligands * batch_atoms;
      for (int ligand_index{0}; ligand_index < batch_ligands; ++ligand_index)
        apply_ligand_placement(*batch.molecules[ligand_index], placement, protein_center);

      load_num_rotamers(batch, this->scratch);
      load_num_atoms(batch, this->scratch);

      const int scores_per_ligand =
          std::max(1, static_cast<int>((*this->scratch).configuration.population_number));

      auto &score_b = (*this->scratch).template get<buffer_data_type::SCORES>();
      load_scratchs<queue_type>(batch, this->scratch, scores_per_ligand);
      if (!score_b.is_valid() ||
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

      auto &prot_grid_maps = (*device_scratch).template get<buffer_data_type::PROT_GRID_MAPS>();
      const fp_type *grid_maps = (*device_scratch).get_queue()->obj_required()
                                     ? prot_grid_maps.dev_pointer()
                                     : prot_grid_maps.host_pointer();
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

      args = {};
      args.batch_ligands = batch_ligands;
      args.scores_per_ligand = scores_per_ligand;
      args.batch_atom_stride = batch_atoms;
      args.map_index_x = map_index_x;
      args.map_index_xy = map_index_xy;
      args.map_index_xyz = map_index_xyz;
      args.vols = vols_b;
      args.solpars = solpars_b;
      args.charges = charges_b;
      args.num_atoms = num_atoms_b;
      args.num_rotamers = num_rotamers_b;
      args.num_nonbonds = num_nonbonds_b;
      args.nonbond_a1 = nonbond_a1_b;
      args.nonbond_a2 = nonbond_a2_b;
      args.nonbond_cA = nonbond_cA_b;
      args.nonbond_cB = nonbond_cB_b;
      args.nonbond_xB = nonbond_xB_b;
      args.grid_maps = grid_maps;
      args.atom_map_offsets = map_offsets_b;
      args.minimum = minimum;
      args.maximum = maximum;
      args.center = center;
      args.scores = scores_b;
      args.x = x_scratch_b;
      args.y = y_scratch_b;
      args.z = z_scratch_b;
    }


    void operator()() override;

    [[nodiscard]] adt_score_args kernel_args() const {
      return args;
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
      mem += sizeof(fp_type) * 3 * max_atoms;       // vols, solpars, charges
      mem += sizeof(int) * max_atoms;               // map_offsets
      mem += sizeof(int);                           // num_nonbond sentinel
      mem += sizeof(int) * 2 * non_bonds_atoms;     // nonbond_a1, nonbond_a2
      mem += sizeof(fp_type) * 2 * non_bonds_atoms; // nonbond_cA, nonbond_cB
      mem += sizeof(int) * non_bonds_atoms;         // nonbond_xB
      return mem;
    }

    static std::size_t get_ligand_mem(const int max_atoms, const knobs conf) {
      return get_shared_ligand_mem(max_atoms, conf) + get_private_ligand_mem(max_atoms, conf);
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
    ligand_placement placement;
    point3D protein_center;
    adt_score_args args{};

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

  template<typename queue_t>
  struct scoring_backend<queue_t, scoring_function::ADT> {
    using stage_type = adt_score<queue_t>;
    using arguments = adt_score_args;
    using device_kernel = typename stage_type::device_kernel;
  };
} // namespace mudock
