#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/chem/autodock_ligand.hpp>
#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/cpp_implementation/center_of_mass.hpp>
#include <mudock/cpp_implementation/geometric_transformations.hpp>
#include <mudock/cpp_implementation/mutate.hpp>
#include <mudock/cpp_implementation/weed_bonds.hpp>
#include <mudock/sycl_implementation/evaluate_fitness.hpp>
#include <mudock/sycl_implementation/virtual_screen.hpp>

#define BUCKET_MULTIPLIER 3

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
        ligand_fragments_start(queue),
        frag_start_atom_indices(queue),
        frag_stop_atom_indices(queue),
        frag_indices_start(queue),
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

  template<int MAX_ATOMS>
  struct evaluate_fitness_kernel_tag {}; // just a name, no members

  void virtual_screen_sycl::operator()(batch &incoming_batch) {
    const std::size_t batch_atoms   = incoming_batch.batch_max_atoms;
    const std::size_t batch_ligands = incoming_batch.num_ligands;
    // Resize data structures
    const std::size_t tot_atoms_in_batch          = batch_ligands * batch_atoms;
    const std::size_t batch_rotamers              = batch_atoms - 3;
    const std::size_t batch_non_bonds             = batch_ligands * batch_atoms * batch_atoms;
    const std::size_t tot_rotamers_atoms_in_batch = tot_atoms_in_batch * batch_rotamers;
    const std::size_t tot_rotamers_in_batch       = batch_ligands * batch_rotamers;
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
    ligand_fragments_start.alloc(batch_ligands + 1);
    ligand_fragments_start()[0] = 0;
    frag_start_atom_indices.alloc(tot_rotamers_in_batch);
    frag_stop_atom_indices.alloc(tot_rotamers_in_batch);
    frag_indices_start.alloc(batch_ligands + 1);
    frag_indices_start()[0] = 0;
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
    index_nonbonds()[0] = 0;
    nonbond_a1.alloc(batch_non_bonds);
    nonbond_a2.alloc(batch_non_bonds);
    nonbond_cA.alloc(batch_non_bonds);
    nonbond_cB.alloc(batch_non_bonds);
    nonbond_xB.alloc(batch_non_bonds);
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
      const int num_atoms       = adt_ligand.get_num_atoms();
      ligand_num_atoms()[index] = num_atoms;
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

      std::memcpy((void *) (original_ligand_x() + stride_atoms), x_p, num_atoms * sizeof(fp_type));
      std::memcpy((void *) (original_ligand_y() + stride_atoms), y_p, num_atoms * sizeof(fp_type));
      std::memcpy((void *) (original_ligand_z() + stride_atoms), z_p, num_atoms * sizeof(fp_type));

      // Randomly initialize the population
      const auto num_rotamers      = adt_ligand.get_num_rotatable_bonds();
      ligand_num_rotamers()[index] = num_rotamers;
      assert(batch_rotamers > ligand.get()->num_rotamers());

      std::memcpy((void *) (ligand_fragments() + ligand_fragments_start()[index]),
                  adt_ligand.get_fragments_masks(),
                  num_atoms * num_rotamers * sizeof(fp_type));
      ligand_fragments_start()[index + 1] = ligand_fragments_start()[index] + (num_atoms * num_rotamers);
      std::memcpy((void *) (frag_start_atom_indices() + frag_indices_start()[index]),
                  adt_ligand.get_fragmets_starts(),
                  num_rotamers * sizeof(fp_type));
      std::memcpy((void *) (frag_stop_atom_indices() + frag_indices_start()[index]),
                  adt_ligand.get_fragments_stops(),
                  num_rotamers * sizeof(fp_type));
      frag_indices_start()[index + 1] = frag_indices_start()[index] + num_rotamers;

      const auto non_bond_size = adt_ligand.get_non_bond_size();
      std::memcpy((void *) (nonbond_a1() + index_nonbonds()[index]),
                  adt_ligand.get_non_bond_A(),
                  non_bond_size * sizeof(int));
      std::memcpy((void *) (nonbond_a2() + index_nonbonds()[index]),
                  adt_ligand.get_non_bond_B(),
                  non_bond_size * sizeof(int));
      std::memcpy((void *) (nonbond_cA() + index_nonbonds()[index]),
                  adt_ligand.get_non_bond_cA(),
                  non_bond_size * sizeof(fp_type));
      std::memcpy((void *) (nonbond_cB() + index_nonbonds()[index]),
                  adt_ligand.get_non_bond_cB(),
                  non_bond_size * sizeof(fp_type));
      std::memcpy((void *) (nonbond_xB() + index_nonbonds()[index]),
                  adt_ligand.get_non_bond_xB(),
                  non_bond_size * sizeof(int));
      index_nonbonds()[index + 1] = index_nonbonds()[index] + non_bond_size;

      // Autodock typing
      std::memcpy((void *) (ligand_vol() + stride_atoms),
                  adt_ligand.get_ligand_vol(),
                  num_atoms * sizeof(fp_type));
      std::memcpy((void *) (ligand_solpar() + stride_atoms),
                  adt_ligand.get_ligand_solpar(),
                  num_atoms * sizeof(fp_type));
      std::memcpy((void *) (ligand_charge() + stride_atoms),
                  adt_ligand.get_ligand_charge(),
                  num_atoms * sizeof(fp_type));

      std::memcpy((void *) (map_texture_index() + stride_atoms),
                  adt_ligand.get_atom_map_offsets(),
                  num_atoms * sizeof(int));
    }

    // Copy in
    const auto copy_frag_size         = ligand_fragments_start()[batch_ligands];
    const auto copy_frag_indices_size = frag_indices_start()[batch_ligands];
    const auto copy_nonbond_size      = index_nonbonds()[batch_ligands];

    ligand_num_atoms.copy_host2device();
    original_ligand_x.copy_host2device();
    original_ligand_y.copy_host2device();
    original_ligand_z.copy_host2device();
    ligand_fragments.copy_host2device(copy_frag_size);
    ligand_fragments_start.copy_host2device();
    ligand_num_rotamers.copy_host2device();
    frag_start_atom_indices.copy_host2device(copy_frag_indices_size);
    frag_stop_atom_indices.copy_host2device(copy_frag_indices_size);
    frag_indices_start.copy_host2device();
    ligand_vol.copy_host2device();
    ligand_solpar.copy_host2device();
    ligand_charge.copy_host2device();
    map_texture_index.copy_host2device();
    index_nonbonds.copy_host2device();
    nonbond_a1.copy_host2device(copy_nonbond_size);
    nonbond_a2.copy_host2device(copy_nonbond_size);
    nonbond_cA.copy_host2device(copy_nonbond_size);
    nonbond_cB.copy_host2device(copy_nonbond_size);
    nonbond_xB.copy_host2device(copy_nonbond_size);

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
          const auto num_gen        = configuration.num_generations;
          const auto tournament_len = configuration.tournament_length;
          const auto mutation_prob  = configuration.mutation_prob;
          const auto pop_num        = configuration.population_number;

          const auto size_x          = adt_protein.get_size_x();
          const auto size_xy         = adt_protein.get_size_xy();
          const auto size_xyz        = adt_protein.get_size_xyz();
          const auto original_x      = original_ligand_x.dev_pointer();
          const auto original_y      = original_ligand_y.dev_pointer();
          const auto original_z      = original_ligand_z.dev_pointer();
          auto scratch_x             = scratch_ligand_x.dev_pointer();
          auto scratch_y             = scratch_ligand_y.dev_pointer();
          auto scratch_z             = scratch_ligand_z.dev_pointer();
          const auto vol             = ligand_vol.dev_pointer();
          const auto solpar          = ligand_solpar.dev_pointer();
          const auto charge          = ligand_charge.dev_pointer();
          const auto nonbonds        = index_nonbonds.dev_pointer();
          const auto a1              = nonbond_a1.dev_pointer();
          const auto a2              = nonbond_a2.dev_pointer();
          const auto ca              = nonbond_cA.dev_pointer();
          const auto cb              = nonbond_cB.dev_pointer();
          const auto xb              = nonbond_xB.dev_pointer();
          const auto atoms           = ligand_num_atoms.dev_pointer();
          const auto rotamers        = ligand_num_rotamers.dev_pointer();
          const auto fragments       = ligand_fragments.dev_pointer();
          const auto fragments_start = ligand_fragments_start.dev_pointer();
          const auto start           = frag_start_atom_indices.dev_pointer();
          const auto stop            = frag_stop_atom_indices.dev_pointer();
          const auto indices_start   = frag_indices_start.dev_pointer();
          auto chromo                = chromosomes.dev_pointer();
          const auto min             = adt_protein.get_min();
          const auto max             = adt_protein.get_max();
          const auto center          = adt_protein.get_center();
          const auto maps            = (*dev).get_tex_dev_pointer();
          const auto map_indexes     = map_texture_index.dev_pointer();
          auto rand                  = random_states.dev_pointer();
          auto l_scores              = ligand_scores.dev_pointer();
          auto best_chromo           = best_chromosomes.dev_pointer();
          sycl::local_accessor<fp_type> shm_acc(sycl::range<1>(shared_mem), h);

          h.parallel_for<evaluate_fitness_kernel_tag<n_atoms>>(
              sycl::nd_range<1>{batch_ligands * subgroup_size, subgroup_size},
              [=, this](sycl::nd_item<1> it) {
                evaluate_fitness<n_atoms>(num_gen,
                                          tournament_len,
                                          mutation_prob,
                                          pop_num,
                                          population_stride,
                                          batch_atoms,
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
                                          fragments_start,
                                          start,
                                          stop,
                                          indices_start,
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
                                     *(best_chromosomes() + index),
                                     num_atoms,
                                     num_rotamers,
                                     adt_ligand.get_fragments_masks(),
                                     adt_ligand.get_fragmets_starts(),
                                     adt_ligand.get_fragments_stops());
      ligand->properties.assign(property_type::SCORE, std::to_string(ligand_scores()[index]));
    }
  }

  template<int MAX_ATOMS>
  int get_evaluate_fitness_batch(const sycl::device &dev) {
    const sycl::context ctx{dev};
    const int compute_units = dev.get_info<sycl::info::device::max_compute_units>();
    const int subgroup_size = dev.get_info<sycl::info::device::sub_group_sizes>()[0];

    // Fetch kernel device-specific limits for this device
    const auto kid         = sycl::get_kernel_id<evaluate_fitness_kernel_tag<MAX_ATOMS>>();
    const auto kb          = sycl::get_kernel_bundle<sycl::bundle_state::executable>(ctx, {dev}, {kid});
    const sycl::kernel krn = kb.get_kernel(kid);

    // Maximum work-group size the device allows for this kernel
    const std::size_t max_wg_size = krn.get_info<sycl::info::kernel_device_specific::work_group_size>(dev);

    // Preferred multiple (warp/wavefront alignment); useful when picking BLOCK_SIZE
    const std::size_t preferred_multiple =
        krn.get_info<sycl::info::kernel_device_specific::preferred_work_group_size_multiple>(dev);

    // Upper bound on concurrently resident work-groups per CU from wg-size alone
    const std::size_t wg_per_cu_cap = std::max<std::size_t>(1, max_wg_size / subgroup_size);

    // Portable heuristic: try to keep 2–4 work-groups per CU if possible.
    // You can tune this number per backend/workload.
    const std::size_t target_wg_per_cu = std::min<std::size_t>(wg_per_cu_cap, 4);

    return static_cast<int>(target_wg_per_cu * compute_units);
  }

  int compute_batch_size_vs(const sycl::device &d, const int num_atoms) {
    // populate the bucket dimension
    int bucket_size{0};
    constexpr_for<0, reorder_buffer::atoms_clusters.size(), 1>([&](const auto atoms_index) {
      const auto n_atoms = reorder_buffer::atoms_clusters[atoms_index];
      if (num_atoms == n_atoms)
        bucket_size = get_evaluate_fitness_batch<n_atoms>(d);
    });
    // TODO check if it can be made a compile error
    if (bucket_size == 0)
      throw std::runtime_error(
          "Compilation error: there is a bucket of atoms number which it is not handled.");
    return bucket_size * BUCKET_MULTIPLIER;
  }
} // namespace mudock
