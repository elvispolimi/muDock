#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/chem/autodock_ligand.hpp>
#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/cpp_implementation/center_of_mass.hpp>
#include <mudock/cpp_implementation/geometric_transformations.hpp>
#include <mudock/cpp_implementation/mutate.hpp>
#include <mudock/cpp_implementation/weed_bonds.hpp>
#include <mudock/sycl_implementation/evaluate_fitness.hpp>
#include <mudock/sycl_implementation/virtual_screen.hpp>

namespace mudock {
  // TODO create a single conf file for all implementations
  static constexpr std::size_t max_non_bonds{1 << 26};
  static constexpr std::size_t max_rotamers_per_ligand{64};

  virtual_screen_sycl::virtual_screen_sycl(const knobs k, const std::shared_ptr<const device> dev)
      : configuration(k),
        dev(dev),
        queue(dev->get_queue()),
        subgroup_size(dev->get_sub_group_size()),
        original_ligand_x(queue),
        original_ligand_y(queue),
        original_ligand_z(queue),
        scratch_ligand_x(queue),
        scratch_ligand_y(queue),
        scratch_ligand_z(queue),
        ligand_vol(queue),
        ligand_solpar(queue),
        ligand_charge(queue),
        ligand_num_atoms(queue),
        ligand_num_rotamers(queue),
        ligand_fragments(queue),
        frag_start_atom_indices(queue),
        frag_stop_atom_indices(queue),
        index_nonbonds(queue),
        nonbond_a1(queue),
        nonbond_a2(queue),
        nonbond_xB(queue),
        nonbond_cA(queue),
        nonbond_cB(queue),
        map_texture_index(queue),
        ligand_scores(queue),
        chromosomes(queue),
        best_chromosomes(queue),
        random_states(queue) {}

  void virtual_screen_sycl::operator()(batch &incoming_batch) {
    const std::size_t batch_atoms   = incoming_batch.batch_max_atoms;
    const std::size_t batch_ligands = incoming_batch.num_ligands;
    // Resize data structures
    const std::size_t tot_atoms_in_batch = batch_ligands * batch_atoms;
    // const std::size_t tot_atoms_in_population     = tot_atoms_in_batch * configuration.population_number;
    const std::size_t tot_rotamers_atoms_in_batch = tot_atoms_in_batch * batch_rotamers;
    const std::size_t tot_rotamers_in_batch       = batch_ligands * max_rotam\max_rotamers_per_ligand;
    const std::size_t batch_nonbonds              = batch_ligands * max_non_bonds;
    // Use double buffering on the GPU for actual and next population at each iteration
    const std::size_t population_stride = configuration.population_number * 2;
    const auto &adt_protein             = (*dev).adt_protein;
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
    // TODO check if it is required
    ligand_num_atoms.alloc(batch_ligands);
    ligand_num_rotamers.alloc(batch_ligands);
    ligand_scores.alloc(batch_ligands);
    best_chromosomes.alloc(batch_ligands);
    // Bonds
    index_nonbonds.alloc(batch_ligands + 1);
    index_nonbonds.host_pointer()[0] = 0;
    nonbond_a1.alloc(max_non_bonds);
    nonbond_a2.alloc(max_non_bonds);
    nonbond_cA.alloc(max_non_bonds);
    nonbond_cB.alloc(max_non_bonds);
    nonbond_xB.alloc(max_non_bonds);
    // GA Data structures
    chromosomes.alloc(population_stride * batch_ligands);
    // Support data precomputation
    map_texture_index.alloc(tot_atoms_in_batch);

    // Copy data
    std::vector<autodock_ligand> vector_adt_ligands;
    for (std::size_t index{0}; index < batch_ligands; ++index) {
      auto &ligand = incoming_batch.molecules[index];
      vector_adt_ligands.emplace_back(*ligand);
      auto &adt_ligand = vector_adt_ligands.back();
      adt_ligand.update_offsets(adt_protein);
      const int stride_atoms = index * batch_atoms;
      // Atoms and bonds
      const int num_atoms                    = adt_ligand.get_num_atoms();
      ligand_num_atoms.host_pointer()[index] = num_atoms;
      // TODO bonds
      // Place the molecule to the center of the target protein
      const auto x = adt_ligand.get_ligand_x(), y = adt_ligand.get_ligand_y(), z = adt_ligand.get_ligand_z();
      auto x_p = adt_ligand.get_ligand_x_p(), y_p = adt_ligand.get_ligand_y_p(),
           z_p = adt_ligand.get_ligand_z_p();

      const auto ligand_center_of_mass = compute_center_of_mass(x, y, z);
      const auto offset                = adt_protein.get_center() - ligand_center_of_mass;
      translate_molecule<cpu_vectorization::AUTO>(x_p,
                                                  y_p,
                                                  z_p,
                                                  num_atoms,
                                                  offset.x(),
                                                  offset.y(),
                                                  offset.z());

      std::memcpy((void *) (original_ligand_x.host_pointer() + stride_atoms),
                  x_p,
                  num_atoms * sizeof(fp_type));
      std::memcpy((void *) (original_ligand_y.host_pointer() + stride_atoms),
                  y_p,
                  num_atoms * sizeof(fp_type));
      std::memcpy((void *) (original_ligand_z.host_pointer() + stride_atoms),
                  z_p,
                  num_atoms * sizeof(fp_type));

      // Randomly initialize the population
      const auto num_rotamers                   = adt_ligand.get_num_rotatable_bonds();
      ligand_num_rotamers.host_pointer()[index] = num_rotamers;
      const int stride_masks                    = index * max_rotam\max_rotamers_per_ligand * batch_atoms;
      const int stride_rotamers                 = index * max_rotam\max_rotamers_per_ligand;
      assert(max_rotam\max_rotamers_per_ligand > ligand.get()->num_rotamers());

      std::memcpy((void *) (ligand_fragments.host_pointer() + stride_masks),
                  adt_ligand.get_fragments_masks(),
                  num_atoms * num_rotamers * sizeof(fp_type));
      std::memcpy((void *) (frag_start_atom_indices.host_pointer() + stride_rotamers),
                  adt_ligand.get_fragmets_starts(),
                  num_rotamers * sizeof(fp_type));
      std::memcpy((void *) (frag_stop_atom_indices.host_pointer() + stride_rotamers),
                  adt_ligand.get_fragments_stops(),
                  num_rotamers * sizeof(fp_type));

      const auto non_bond_size = adt_ligand.get_non_bond_size();
      std::memcpy((void *) (nonbond_a1.host_pointer() + index_nonbonds.host_pointer()[index]),
                  adt_ligand.get_non_bond_A(),
                  non_bond_size * sizeof(int));
      std::memcpy((void *) (nonbond_a2.host_pointer() + index_nonbonds.host_pointer()[index]),
                  adt_ligand.get_non_bond_B(),
                  non_bond_size * sizeof(int));
      std::memcpy((void *) (nonbond_cA.host_pointer() + index_nonbonds.host_pointer()[index]),
                  adt_ligand.get_non_bond_cA(),
                  non_bond_size * sizeof(fp_type));
      std::memcpy((void *) (nonbond_cB.host_pointer() + index_nonbonds.host_pointer()[index]),
                  adt_ligand.get_non_bond_cB(),
                  non_bond_size * sizeof(fp_type));
      std::memcpy((void *) (nonbond_xB.host_pointer() + index_nonbonds.host_pointer()[index]),
                  adt_ligand.get_non_bond_xB(),
                  non_bond_size * sizeof(int));
      index_nonbonds.host_pointer()[index + 1] = index_nonbonds.host_pointer()[index] + non_bond_size;

      // Autodock typing
      std::memcpy((void *) (ligand_vol.host_pointer() + stride_atoms),
                  adt_ligand.get_ligand_vol(),
                  num_atoms * sizeof(fp_type));
      std::memcpy((void *) (ligand_solpar.host_pointer() + stride_atoms),
                  adt_ligand.get_ligand_solpar(),
                  num_atoms * sizeof(fp_type));
      std::memcpy((void *) (ligand_charge.host_pointer() + stride_atoms),
                  adt_ligand.get_ligand_charge(),
                  num_atoms * sizeof(fp_type));

      std::memcpy((void *) (map_texture_index.host_pointer() + stride_atoms),
                  adt_ligand.get_atom_map_offsets(),
                  num_atoms * sizeof(int));
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
    map_texture_index.copy_host2device();
    index_nonbonds.copy_host2device();
    nonbond_a1.copy_host2device();
    nonbond_a2.copy_host2device();
    nonbond_cA.copy_host2device();
    nonbond_cB.copy_host2device();
    nonbond_xB.copy_host2device();

    // Setup sycl random
    // A state for each thread
    const int num_threads = batch_ligands * subgroup_size;
    random_states.alloc(num_threads);
    if (configuration.seed.has_value())
      random_states.alloc(num_threads, configuration.seed.value());
    else
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
    constexpr_for<0, reorder_buffer::atoms_clusters.size(), 1>([&, this](const auto atoms_index) {
      const auto n_atoms = reorder_buffer::atoms_clusters[atoms_index];
      if (batch_atoms == n_atoms)
        queue.submit([&, this](sycl::handler &h) {
          auto size_x      = adt_protein.get_size_x();
          auto size_xy     = adt_protein.get_size_xy();
          auto size_xyz    = adt_protein.get_size_xyz();
          auto original_x  = original_ligand_x.dev_pointer();
          auto original_y  = original_ligand_y.dev_pointer();
          auto original_z  = original_ligand_z.dev_pointer();
          auto scratch_x   = scratch_ligand_x.dev_pointer();
          auto scratch_y   = scratch_ligand_y.dev_pointer();
          auto scratch_z   = scratch_ligand_z.dev_pointer();
          auto vol         = ligand_vol.dev_pointer();
          auto solpar      = ligand_solpar.dev_pointer();
          auto charge      = ligand_charge.dev_pointer();
          auto nonbonds    = index_nonbonds.dev_pointer();
          auto a1          = nonbond_a1.dev_pointer();
          auto a2          = nonbond_a2.dev_pointer();
          auto ca          = nonbond_cA.dev_pointer();
          auto cb          = nonbond_cB.dev_pointer();
          auto xb          = nonbond_xB.dev_pointer();
          auto atoms       = ligand_num_atoms.dev_pointer();
          auto rotamers    = ligand_num_rotamers.dev_pointer();
          auto fragments   = ligand_fragments.dev_pointer();
          auto start       = frag_start_atom_indices.dev_pointer();
          auto stop        = frag_stop_atom_indices.dev_pointer();
          auto chromo      = chromosomes.dev_pointer();
          auto min         = adt_protein.get_min();
          auto max         = adt_protein.get_max();
          auto center      = adt_protein.get_center();
          auto maps        = (*dev).get_tex_dev_pointer();
          auto map_indexes = map_texture_index.dev_pointer();
          auto rand        = random_states.dev_pointer();
          auto l_scores    = ligand_scores.dev_pointer();
          auto best_chromo = best_chromosomes.dev_pointer();
          sycl::local_accessor<fp_type> shm_acc(sycl::range<1>(shared_mem), h);

          h.parallel_for(sycl::nd_range<1>{batch_ligands * subgroup_size, subgroup_size},
                         [=, this](sycl::nd_item<1> it) {
                           evaluate_fitness<n_atoms>(configuration.num_generations,
                                                     configuration.tournament_length,
                                                     configuration.mutation_prob,
                                                     configuration.population_number,
                                                     population_stride,
                                                     batch_atoms,
                                                     max_rotamers_per_ligand,
                                                     size_x,
                                                     size_xy,
                                                     size_xyz,
                                                     original_x,
                                                     original_y,
                                                     original_z,
                                                     scratch_x,
                                                     scratch_y,
                                                     scratch_z,
                                                     vol,
                                                     solpar,
                                                     charge,
                                                     nonbonds,
                                                     a1,
                                                     a2,
                                                     ca,
                                                     cb,
                                                     xb,
                                                     atoms,
                                                     rotamers,
                                                     fragments,
                                                     start,
                                                     stop,
                                                     chromo,
                                                     min,
                                                     max,
                                                     center,
                                                     maps,
                                                     map_indexes,
                                                     rand,
                                                     shm_acc,
                                                     l_scores,
                                                     best_chromo,
                                                     it);
                         });
        });
    });

    // Copy back chromosomes and scores
    best_chromosomes.copy_device2host();
    ligand_scores.copy_device2host();

    queue.wait();

    // update the ligand position with the best one that we found
    for (std::size_t index{0}; index < batch_ligands; ++index) {
      auto &adt_ligand        = vector_adt_ligands[index];
      auto &ligand            = incoming_batch.molecules[index];
      const int num_atoms     = adt_ligand.get_num_atoms();
      const auto num_rotamers = adt_ligand.get_num_rotatable_bonds();

      // for (auto &ligand: std::span(incoming_batch.molecules.data(), incoming_batch.num_ligands)) {
      // Reset the random number generator to improve consistency
      apply<cpu_vectorization::AUTO>(adt_ligand.get_ligand_x_p(),
                                     adt_ligand.get_ligand_y_p(),
                                     adt_ligand.get_ligand_z_p(),
                                     *(best_chromosomes.host_pointer() + index),
                                     num_atoms,
                                     num_rotamers,
                                     adt_ligand.get_fragments_masks(),
                                     adt_ligand.get_fragmets_starts(),
                                     adt_ligand.get_fragments_stops());
      ligand->properties.assign(property_type::SCORE, std::to_string(ligand_scores.host_pointer()[index]));
    }
  }
} // namespace mudock
