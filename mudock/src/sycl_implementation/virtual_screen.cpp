#include <mudock/cpp_implementation/center_of_mass.hpp>
#include <mudock/cpp_implementation/geometric_transformations.hpp>
#include <mudock/cpp_implementation/mutate.hpp>
#include <mudock/cpp_implementation/weed_bonds.hpp>
#include <mudock/cuda_implementation/map_textures.cuh>
#include <mudock/sycl_implementation/evaluate_fitness.hpp>
#include <mudock/sycl_implementation/virtual_screen.hpp>

namespace mudock {
  // TODO create a single conf file for all implementations
  static constexpr std::size_t max_non_bonds{1024 * 10};
  // Grid spacing fixed to 0.5 Angstrom
  static constexpr fp_type inv_spacing{2};

  void init_texture_memory(const grid_map &map, sycl_object<fp_type> &tex_obj) {
    tex_obj.alloc(map.index.size_x() * map.index.size_y() * map.index.size_z());
    tex_obj.copy_host2device(map.data());
  }

  virtual_screen_sycl::virtual_screen_sycl(const knobs k,
                                           std::shared_ptr<const grid_atom_mapper> &grid_atom_maps,
                                           std::shared_ptr<const grid_map> &electro_map,
                                           std::shared_ptr<const grid_map> &desolv_map,
                                           sycl::queue &_queue)
      : queue(_queue),
        configuration(k),
        original_ligand_x(queue),
        original_ligand_y(queue),
        original_ligand_z(queue),
        scratch_ligand_x(queue),
        scratch_ligand_y(queue),
        scratch_ligand_z(queue),
        ligand_vol(queue),
        ligand_solpar(queue),
        ligand_charge(queue),
        ligand_Rij_hb(queue),
        ligand_Rii(queue),
        ligand_epsij_hb(queue),
        ligand_epsii(queue),
        ligand_num_hbond(queue),
        ligand_num_atoms(queue),
        ligand_num_rotamers(queue),
        ligand_fragments(queue),
        frag_start_atom_indices(queue),
        frag_stop_atom_indices(queue),
        num_nonbonds(queue),
        nonbond_a1(queue),
        nonbond_a2(queue),
        map_texture_index(queue),
        ligand_scores(queue),
        chromosomes(queue),
        best_chromosomes(queue),
        center(electro_map.get()->center),
        minimum_coord(electro_map.get()->minimum_coord),
        maximum_coord(electro_map.get()->maximum_coord),
        electro_tex(queue),
        desolv_tex(queue),
        atom_texs(queue),
        random_states(queue) {
    assert(electro_map.get()->index == desolv_map.get()->index &&
           desolv_map.get()->index == grid_atom_maps.get()->get_index());
    init_texture_memory(*electro_map.get(), electro_tex);

    init_texture_memory(*desolv_map.get(), desolv_tex);

    atom_texs.wrappers_pointer.alloc(num_device_map_textures());
    atom_texs.wrappers.reserve(num_device_map_textures());
    for (int index{0}; index < num_device_map_textures(); ++index) {
      atom_texs.wrappers.emplace_back(queue);
      sycl_object<fp_type> &atom_tex = atom_texs.wrappers.back();
      const grid_map &grid_atom =
          grid_atom_maps.get()->get_atom_map(autodock_type_from_map(static_cast<device_map_textures>(index)));
      init_texture_memory(grid_atom, atom_tex);
      atom_texs.wrappers_pointer.host[index] = atom_tex.dev_pointer();
    }
    atom_texs.wrappers_pointer.copy_host2device();

    // Get the device associated with the queue
    sycl::device device = queue.get_device();

    // Get the maximum work-group size for the device
    subgroup_size = static_cast<int>(device.get_info<sycl::info::device::sub_group_sizes>().at(0));
  }

  void virtual_screen_sycl::operator()(batch &incoming_batch) {
    const std::size_t batch_atoms    = incoming_batch.batch_max_atoms;
    const std::size_t batch_rotamers = incoming_batch.batch_max_rotamers;
    const std::size_t batch_ligands  = incoming_batch.num_ligands;
    // Resize data structures
    const std::size_t tot_atoms_in_batch = batch_ligands * batch_atoms;
    // const std::size_t tot_atoms_in_population     = tot_atoms_in_batch * configuration.population_number;
    const std::size_t tot_rotamers_atoms_in_batch = tot_atoms_in_batch * batch_rotamers;
    const std::size_t tot_rotamers_in_batch       = batch_ligands * batch_rotamers;
    const std::size_t batch_nonbonds              = batch_ligands * max_non_bonds;
    // Use double buffering on the GPU for actual and next population at each iteration
    const std::size_t population_stride = configuration.population_number * 2;
    original_ligand_x.alloc(tot_atoms_in_batch);
    original_ligand_y.alloc(tot_atoms_in_batch);
    original_ligand_z.alloc(tot_atoms_in_batch);
    scratch_ligand_x.alloc(tot_atoms_in_batch);
    scratch_ligand_y.alloc(tot_atoms_in_batch);
    scratch_ligand_z.alloc(tot_atoms_in_batch);
    ligand_fragments.alloc(tot_rotamers_atoms_in_batch);
    frag_start_atom_indices.alloc(tot_rotamers_in_batch);
    frag_stop_atom_indices.alloc(tot_rotamers_in_batch);
    ligand_vol.alloc(tot_atoms_in_batch);
    ligand_solpar.alloc(tot_atoms_in_batch);
    ligand_charge.alloc(tot_atoms_in_batch);
    ligand_Rij_hb.alloc(tot_atoms_in_batch);
    ligand_Rii.alloc(tot_atoms_in_batch);
    ligand_epsij_hb.alloc(tot_atoms_in_batch);
    ligand_epsii.alloc(tot_atoms_in_batch);
    // TODO check if it is required
    ligand_num_hbond.alloc(tot_atoms_in_batch);
    ligand_num_atoms.alloc(batch_ligands);
    ligand_num_rotamers.alloc(batch_ligands);
    ligand_scores.alloc(batch_ligands);
    best_chromosomes.alloc(batch_ligands);
    // Bonds
    num_nonbonds.alloc(batch_ligands);
    nonbond_a1.alloc(batch_nonbonds);
    nonbond_a2.alloc(batch_nonbonds);
    // GA Data structures
    chromosomes.alloc(population_stride * batch_ligands);
    // Support data precomputation
    map_texture_index.alloc(tot_atoms_in_batch);

    // Copy data
    std::size_t index{0};
    // TODO pragma OMP improve performance
    // Keep the fragments for the output
    std::vector<fragments<static_containers>> batch_fragments;
    batch_fragments.reserve(batch_ligands);
    // Support data structures
    grid<uint_fast8_t, index2D> nbmatrix{{static_cast<int>(batch_atoms), static_cast<int>(batch_atoms)}};
    std::vector<non_bond_parameter> non_bond_list;
    non_bond_list.reserve(max_non_bonds);
    for (auto &ligand: std::span(incoming_batch.molecules.data(), batch_ligands)) {
      const int stride_atoms = index * batch_atoms;
      // Atoms and bonds
      const int num_atoms                    = ligand.get()->num_atoms();
      ligand_num_atoms.host_pointer()[index] = num_atoms;
      // TODO bonds
      // Place the molecule to the center of the target protein
      const auto x = ligand.get()->get_x(), y = ligand.get()->get_y(), z = ligand.get()->get_z();
      const auto ligand_center_of_mass = compute_center_of_mass(x, y, z);
      translate_molecule(x,
                         y,
                         z,
                         center.x - ligand_center_of_mass.x,
                         center.y - ligand_center_of_mass.y,
                         center.z - ligand_center_of_mass.z);

      std::memcpy((void *) (original_ligand_x.host_pointer() + stride_atoms),
                  x.data(),
                  num_atoms * sizeof(fp_type));
      std::memcpy((void *) (original_ligand_y.host_pointer() + stride_atoms),
                  y.data(),
                  num_atoms * sizeof(fp_type));
      std::memcpy((void *) (original_ligand_z.host_pointer() + stride_atoms),
                  z.data(),
                  num_atoms * sizeof(fp_type));
      // Fragments
      // Find out the rotatable bonds in the ligand
      auto graph = make_graph(ligand.get()->get_bonds());
      // TODO check this assignment
      batch_fragments[index]  = {graph, ligand.get()->get_bonds(), ligand.get()->num_atoms()};
      const auto &l_fragments = batch_fragments[index];

      // Randomly initialize the population
      const auto num_rotamers                   = l_fragments.get_num_rotatable_bonds();
      ligand_num_rotamers.host_pointer()[index] = num_rotamers;
      const int stride_masks                    = index * batch_rotamers * batch_atoms;
      const int stride_rotamers                 = index * batch_rotamers;
      assert(static_cast<int>(batch_rotamers) > ligand.get()->num_rotamers());
      for (int rot = 0; rot < num_rotamers; ++rot) {
        std::memcpy((void *) (ligand_fragments.host_pointer() + stride_masks + rot * batch_atoms),
                    l_fragments.get_mask(rot).data(),
                    num_atoms * sizeof(int));
        const auto [start_index, stop_index]                          = l_fragments.get_rotatable_atoms(rot);
        frag_start_atom_indices.host_pointer()[stride_rotamers + rot] = start_index;
        frag_stop_atom_indices.host_pointer()[stride_rotamers + rot]  = stop_index;
      }

      // Weed bonds
      // TODO can be accelerated on the GPU
      // No need to move data -> already computed on the GPU
      nbmatrix.reset();
      nonbonds(nbmatrix, ligand.get()->get_bonds(), num_atoms);
      non_bond_list.clear();
      weed_bonds(nbmatrix, non_bond_list, num_atoms, l_fragments);
      if constexpr (is_debug())
        if (non_bond_list.size() >= max_non_bonds) {
          throw std::runtime_error("Bond list size exceed maximum value " +
                                   std::to_string(non_bond_list.size()) + ".");
        }

      num_nonbonds.host_pointer()[index] = non_bond_list.size();
      const int stride_nonbonds          = index * max_non_bonds;
      int nonbond_index{0};
      for (auto &bond: non_bond_list) {
        nonbond_a1.host_pointer()[stride_nonbonds + nonbond_index] = bond.a1;
        nonbond_a2.host_pointer()[stride_nonbonds + nonbond_index] = bond.a2;
        ++nonbond_index;
      }

      // Autodock typing
      std::memcpy((void *) (ligand_vol.host_pointer() + stride_atoms),
                  ligand.get()->get_vol().data(),
                  num_atoms * sizeof(fp_type));
      std::memcpy((void *) (ligand_solpar.host_pointer() + stride_atoms),
                  ligand.get()->get_solpar().data(),
                  num_atoms * sizeof(fp_type));
      std::memcpy((void *) (ligand_charge.host_pointer() + stride_atoms),
                  ligand.get()->get_charge().data(),
                  num_atoms * sizeof(fp_type));
      std::memcpy((void *) (ligand_Rij_hb.host_pointer() + stride_atoms),
                  ligand.get()->get_Rij_hb().data(),
                  num_atoms * sizeof(fp_type));
      std::memcpy((void *) (ligand_Rii.host_pointer() + stride_atoms),
                  ligand.get()->get_Rii().data(),
                  num_atoms * sizeof(fp_type));
      std::memcpy((void *) (ligand_epsij_hb.host_pointer() + stride_atoms),
                  ligand.get()->get_epsij_hb().data(),
                  num_atoms * sizeof(fp_type));
      std::memcpy((void *) (ligand_epsii.host_pointer() + stride_atoms),
                  ligand.get()->get_epsii().data(),
                  num_atoms * sizeof(fp_type));
      std::memcpy((void *) (ligand_num_hbond.host_pointer() + stride_atoms),
                  ligand.get()->get_num_hbond().data(),
                  num_atoms * sizeof(int));

      std::size_t atom_index{0};
      for (auto &autodock_t: ligand.get()->get_autodock_type()) {
        const auto map_index = static_cast<int>(map_from_autodock_type(autodock_t));
        map_texture_index.host_pointer()[stride_atoms + atom_index] = map_index;
        ++atom_index;
      }

      ++index;
    }
    // Copy in
    ligand_num_atoms.copy_host2device();
    original_ligand_x.copy_host2device();
    original_ligand_y.copy_host2device();
    original_ligand_z.copy_host2device();
    ligand_fragments.copy_host2device();
    ligand_num_rotamers.copy_host2device();
    frag_start_atom_indices.copy_host2device();
    frag_stop_atom_indices.copy_host2device();
    ligand_vol.copy_host2device();
    ligand_solpar.copy_host2device();
    ligand_charge.copy_host2device();
    ligand_Rij_hb.copy_host2device();
    ligand_Rii.copy_host2device();
    ligand_epsij_hb.copy_host2device();
    ligand_epsii.copy_host2device();
    ligand_num_hbond.copy_host2device();
    map_texture_index.copy_host2device();
    num_nonbonds.copy_host2device();
    nonbond_a1.copy_host2device();
    nonbond_a2.copy_host2device();

    // Setup sycl random
    // A state for each thread
    const int num_threads = batch_ligands * subgroup_size;
    random_states.alloc(num_threads);
    // Simulate the population evolution for the given amount of time
    const auto num_generations = configuration.num_generations;
    // TODO checks if everything fit into shared memory
    // The shared memory contains:
    // - each chromosome's score at last population evaluation
    // Max due to the reduction at the end to find the highest scores per each chromosome
    const std::size_t min_energy_reduction_s_mem =
        std::max(configuration.population_number, static_cast<std::size_t>(subgroup_size)) * sizeof(fp_type);
    const std::size_t shared_mem = min_energy_reduction_s_mem;
    queue
        .submit([&](sycl::handler &h) {
          const auto *original_ligand_x_k       = original_ligand_x.dev_pointer();
          const auto *original_ligand_y_k       = original_ligand_y.dev_pointer();
          const auto *original_ligand_z_k       = original_ligand_z.dev_pointer();
          auto *scratch_ligand_x_k              = scratch_ligand_x.dev_pointer();
          auto *scratch_ligand_y_k              = scratch_ligand_y.dev_pointer();
          auto *scratch_ligand_z_k              = scratch_ligand_z.dev_pointer();
          const auto *ligand_vol_k              = ligand_vol.dev_pointer();
          const auto *ligand_solpar_k           = ligand_solpar.dev_pointer();
          const auto *ligand_charge_k           = ligand_charge.dev_pointer();
          const auto *ligand_num_hbond_k        = ligand_num_hbond.dev_pointer();
          const auto *ligand_Rij_hb_k           = ligand_Rij_hb.dev_pointer();
          const auto *ligand_Rii_k              = ligand_Rii.dev_pointer();
          const auto *ligand_epsij_hb_k         = ligand_epsij_hb.dev_pointer();
          const auto *ligand_epsii_k            = ligand_epsii.dev_pointer();
          const auto *num_nonbonds_k            = num_nonbonds.dev_pointer();
          const auto *nonbond_a1_k              = nonbond_a1.dev_pointer();
          const auto *nonbond_a2_k              = nonbond_a2.dev_pointer();
          const auto *ligand_num_atoms_k        = ligand_num_atoms.dev_pointer();
          const auto *ligand_num_rotamers_k     = ligand_num_rotamers.dev_pointer();
          const auto *ligand_fragments_k        = ligand_fragments.dev_pointer();
          const auto *frag_start_atom_indices_k = frag_start_atom_indices.dev_pointer();
          const auto *frag_stop_atom_indices_k  = frag_stop_atom_indices.dev_pointer();
          auto *chromosomes_k                   = chromosomes.dev_pointer();
          const auto *wrappers_pointer_k        = atom_texs.wrappers_pointer.dev_pointer();
          const auto *map_texture_index_k       = map_texture_index.dev_pointer();
          const auto *electro_tex_k             = electro_tex.dev_pointer();
          const auto *desolv_tex_k              = desolv_tex.dev_pointer();
          auto *random_states_k                 = random_states.dev_pointer();
          auto *ligand_scores_k                 = ligand_scores.dev_pointer();
          auto *best_chromosomes_k              = best_chromosomes.dev_pointer();

          sycl::local_accessor<fp_type> shm_acc(sycl::range<1>(shared_mem), h);

          const auto num_generations_l   = num_generations;
          const auto tournament_length_l = configuration.tournament_length;
          const auto mutation_prob_l     = configuration.mutation_prob;
          const auto chromosome_number_l = configuration.population_number;
          const auto chromosome_stride_l = population_stride;
          const auto atom_stride_l       = batch_atoms;
          const auto rotamers_stride_l   = batch_rotamers;
          const auto nonbond_stride_l    = max_non_bonds;
          const auto minimum_l           = minimum_coord;
          const auto maximum_l           = maximum_coord;
          const auto center_l            = center;
          const auto index_l             = index;
          const auto inv_spacing_l       = inv_spacing;

          h.parallel_for(sycl::nd_range<1>{batch_ligands * subgroup_size, subgroup_size},
                         [=](sycl::nd_item<1> it) {
                           evaluate_fitness(num_generations_l,
                                            tournament_length_l,
                                            mutation_prob_l,
                                            chromosome_number_l,
                                            chromosome_stride_l,
                                            atom_stride_l,
                                            rotamers_stride_l,
                                            nonbond_stride_l,
                                            original_ligand_x_k,
                                            original_ligand_y_k,
                                            original_ligand_z_k,
                                            scratch_ligand_x_k,
                                            scratch_ligand_y_k,
                                            scratch_ligand_z_k,
                                            ligand_vol_k,
                                            ligand_solpar_k,
                                            ligand_charge_k,
                                            ligand_num_hbond_k,
                                            ligand_Rij_hb_k,
                                            ligand_Rii_k,
                                            ligand_epsij_hb_k,
                                            ligand_epsii_k,
                                            num_nonbonds_k,
                                            nonbond_a1_k,
                                            nonbond_a2_k,
                                            ligand_num_atoms_k,
                                            ligand_num_rotamers_k,
                                            ligand_fragments_k,
                                            frag_start_atom_indices_k,
                                            frag_stop_atom_indices_k,
                                            chromosomes_k,
                                            minimum_l,
                                            maximum_l,
                                            center_l,
                                            index_l,
                                            inv_spacing_l,
                                            wrappers_pointer_k,
                                            map_texture_index_k,
                                            electro_tex_k,
                                            desolv_tex_k,
                                            random_states_k,
                                            shm_acc.get_pointer(),
                                            ligand_scores_k,
                                            best_chromosomes_k,
                                            it);
                         });
        })
        .wait();

    // Copy back chromosomes and scores
    best_chromosomes.copy_device2host();
    ligand_scores.copy_device2host();

    // update the ligand position with the best one that we found
    index = 0;
    for (auto &ligand: std::span(incoming_batch.molecules.data(), incoming_batch.num_ligands)) {
      // Reset the random number generator to improve consistency
      apply(ligand.get()->get_x(),
            ligand.get()->get_y(),
            ligand.get()->get_z(),
            *(best_chromosomes.host_pointer() + index),
            batch_fragments[index]);
      ligand->properties.assign(property_type::SCORE, std::to_string(ligand_scores.host_pointer()[index]));
      ++index;
    }
  }
} // namespace mudock
