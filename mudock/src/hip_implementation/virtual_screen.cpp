#include <cstddef>
#include <cstring>
#include <hip/hip_runtime.h>
#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/chem/autodock_ligand.hpp>
#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/cpp_implementation/center_of_mass.hpp>
#include <mudock/cpp_implementation/geometric_transformations.hpp>
#include <mudock/cpp_implementation/mutate.hpp>
#include <mudock/grid.hpp>
#include <mudock/hip_implementation/evaluate_fitness.hpp>
#include <mudock/hip_implementation/virtual_screen.hpp>
#include <mudock/utils.hpp>

namespace mudock {
  virtual_screen_hip::virtual_screen_hip(const knobs k, const std::shared_ptr<const device> dev)
      : configuration(k),
        dev(dev),
        stream(dev->get_stream()),
        original_ligand_x(stream()),
        original_ligand_y(stream()),
        original_ligand_z(stream()),
        scratch_ligand_x(stream()),
        scratch_ligand_y(stream()),
        scratch_ligand_z(stream()),
        ligand_vol(stream()),
        ligand_solpar(stream()),
        ligand_charge(stream()),
        ligand_num_hbond(stream()),
        ligand_num_atoms(stream()),
        ligand_num_rotamers(stream()),
        ligand_fragments(stream()),
        ligand_fragments_start(stream()),
        frag_start_atom_indices(stream()),
        frag_stop_atom_indices(stream()),
        frag_indices_start(stream()),
        index_nonbonds(stream()),
        nonbond_a1(stream()),
        nonbond_a2(stream()),
        nonbond_cA(stream()),
        nonbond_cB(stream()),
        nonbond_xB(stream()),
        map_texture_index(stream()),
        ligand_scores(stream()),
        chromosomes(stream()),
        best_chromosomes(stream()),
        hiprand_states(stream()) {}

  void virtual_screen_hip::operator()(batch<autodock_ligand> &incoming_batch) {
    const auto wavefront_size = dev->get_wavefront();

    const std::size_t batch_atoms   = incoming_batch.batch_max_atoms;
    const std::size_t batch_ligands = incoming_batch.num_ligands;
    // Resize data structures
    const std::size_t tot_atoms_in_batch = batch_ligands * batch_atoms;
    const std::size_t batch_rotamers     = batch_atoms - 3;
    const std::size_t batch_non_bonds    = batch_ligands * batch_atoms * batch_atoms;
    // const std::size_t tot_atoms_in_population     = tot_atoms_in_batch * configuration.population_number;
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
    ligand_num_hbond.alloc(tot_atoms_in_batch);
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
    for (std::size_t index{0}; index < batch_ligands; ++index) {
      auto &ligand = *incoming_batch.molecules[index];
      ligand.update_offsets(adt_protein);
      const int stride_atoms = index * batch_atoms;
      // Atoms and bonds
      const int num_atoms       = ligand.num_atoms();
      ligand_num_atoms()[index] = num_atoms;
      // TODO bonds
      // Place the molecule to the center of the target protein
      const auto x = ligand.x(), y = ligand.y(), z = ligand.z();

      const auto ligand_center_of_mass = compute_center_of_mass(x, y, z, num_atoms);
      const auto offset                = adt_protein.get_center() - ligand_center_of_mass;
      translate_molecule<cpu_vectorization::AUTO>(x, y, z, num_atoms, offset.x(), offset.y(), offset.z());

      std::memcpy((void *) (original_ligand_x() + stride_atoms), x, num_atoms * sizeof(fp_type));
      std::memcpy((void *) (original_ligand_y() + stride_atoms), y, num_atoms * sizeof(fp_type));
      std::memcpy((void *) (original_ligand_z() + stride_atoms), z, num_atoms * sizeof(fp_type));

      // Randomly initialize the population
      const auto num_rotamers      = ligand.num_rotamers();
      ligand_num_rotamers()[index] = num_rotamers;
      assert(batch_rotamers > ligand.get()->num_rotamers());

      std::memcpy((void *) (ligand_fragments() + ligand_fragments_start()[index]),
                  ligand.fragments_masks(),
                  num_atoms * num_rotamers * sizeof(fp_type));
      ligand_fragments_start()[index + 1] = ligand_fragments_start()[index] + (num_atoms * num_rotamers);
      std::memcpy((void *) (frag_start_atom_indices() + frag_indices_start()[index]),
                  ligand.fragmets_starts(),
                  num_rotamers * sizeof(fp_type));
      std::memcpy((void *) (frag_stop_atom_indices() + frag_indices_start()[index]),
                  ligand.fragments_stops(),
                  num_rotamers * sizeof(fp_type));
      frag_indices_start()[index + 1] = frag_indices_start()[index] + num_rotamers;

      const auto non_bond_size = ligand.non_bond_size();
      std::memcpy((void *) (nonbond_a1() + index_nonbonds()[index]),
                  ligand.non_bond_A(),
                  non_bond_size * sizeof(int));
      std::memcpy((void *) (nonbond_a2() + index_nonbonds()[index]),
                  ligand.non_bond_B(),
                  non_bond_size * sizeof(int));
      std::memcpy((void *) (nonbond_cA() + index_nonbonds()[index]),
                  ligand.non_bond_cA(),
                  non_bond_size * sizeof(fp_type));
      std::memcpy((void *) (nonbond_cB() + index_nonbonds()[index]),
                  ligand.non_bond_cB(),
                  non_bond_size * sizeof(fp_type));
      std::memcpy((void *) (nonbond_xB() + index_nonbonds()[index]),
                  ligand.non_bond_xB(),
                  non_bond_size * sizeof(int));
      index_nonbonds()[index + 1] = index_nonbonds()[index] + non_bond_size;

      // Autodock typing
      std::memcpy((void *) (ligand_vol() + stride_atoms), ligand.vol(), num_atoms * sizeof(fp_type));
      std::memcpy((void *) (ligand_solpar() + stride_atoms), ligand.solpar(), num_atoms * sizeof(fp_type));
      std::memcpy((void *) (ligand_charge() + stride_atoms), ligand.charge(), num_atoms * sizeof(fp_type));

      std::memcpy((void *) (map_texture_index() + stride_atoms),
                  ligand.atom_map_offsets(),
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
    ligand_num_hbond.copy_host2device();
    map_texture_index.copy_host2device();
    index_nonbonds.copy_host2device();
    nonbond_a1.copy_host2device(copy_nonbond_size);
    nonbond_a2.copy_host2device(copy_nonbond_size);
    nonbond_cA.copy_host2device(copy_nonbond_size);
    nonbond_cB.copy_host2device(copy_nonbond_size);
    nonbond_xB.copy_host2device(copy_nonbond_size);

    // Setup hip random
    // A state for each thread
    const int num_threads = batch_ligands * wavefront_size;
    if (configuration.seed.has_value())
      hiprand_states.alloc(num_threads, configuration.seed.value());
    else
      hiprand_states.alloc(num_threads);

    // Simulate the population evolution for the given amount of time
    const auto num_generations = configuration.num_generations;
    // The shared memory contains:
    // - each chromosome's score at last population evaluation
    // Max due to the reduction at the end to find the highest scores per each chromosome
    const std::size_t min_energy_reduction_s_mem =
        std::max(configuration.population_number, static_cast<std::size_t>(wavefront_size)) * sizeof(fp_type);
    const std::size_t shared_mem = min_energy_reduction_s_mem;

    constexpr_for<0, reorder_buffer<autodock_ligand>::atoms_clusters.size(), 1>([&](const auto atoms_index) {
      const auto n_atoms = reorder_buffer<autodock_ligand>::atoms_clusters[atoms_index];
      if (batch_atoms == n_atoms)
        evaluate_fitness<n_atoms>
            <<<batch_ligands, wavefront_size, shared_mem, stream()>>>(num_generations,
                                                                      configuration.tournament_length,
                                                                      configuration.mutation_prob,
                                                                      configuration.population_number,
                                                                      population_stride,
                                                                      batch_atoms,
                                                                      adt_protein.get_size_x(),
                                                                      adt_protein.get_size_xy(),
                                                                      adt_protein.get_size_xyz(),
                                                                      original_ligand_x.dev_pointer(),
                                                                      original_ligand_y.dev_pointer(),
                                                                      original_ligand_z.dev_pointer(),
                                                                      scratch_ligand_x.dev_pointer(),
                                                                      scratch_ligand_y.dev_pointer(),
                                                                      scratch_ligand_z.dev_pointer(),
                                                                      ligand_vol.dev_pointer(),
                                                                      ligand_solpar.dev_pointer(),
                                                                      ligand_charge.dev_pointer(),
                                                                      index_nonbonds.dev_pointer(),
                                                                      nonbond_a1.dev_pointer(),
                                                                      nonbond_a2.dev_pointer(),
                                                                      nonbond_cA.dev_pointer(),
                                                                      nonbond_cB.dev_pointer(),
                                                                      nonbond_xB.dev_pointer(),
                                                                      ligand_num_atoms.dev_pointer(),
                                                                      ligand_num_rotamers.dev_pointer(),
                                                                      ligand_fragments.dev_pointer(),
                                                                      ligand_fragments_start.dev_pointer(),
                                                                      frag_start_atom_indices.dev_pointer(),
                                                                      frag_stop_atom_indices.dev_pointer(),
                                                                      frag_indices_start.dev_pointer(),
                                                                      chromosomes.dev_pointer(),
                                                                      (*dev).get_tex_dev_pointer(),
                                                                      map_texture_index.dev_pointer(),
                                                                      hiprand_states.dev_pointer(),
                                                                      ligand_scores.dev_pointer(),
                                                                      best_chromosomes.dev_pointer());
    });
    // Copy back chromosomes and scores
    best_chromosomes.copy_device2host();
    ligand_scores.copy_device2host();

    MUDOCK_CHECK_KERNELCALL();
    MUDOCK_CHECK(hipStreamSynchronize(stream()));

    // update the ligand position with the best one that we found
    for (std::size_t index{0}; index < batch_ligands; ++index) {
      auto &ligand            = *incoming_batch.molecules[index];
      const int num_atoms     = ligand.num_atoms();
      const auto num_rotamers = ligand.num_rotamers();

      // for (auto &ligand: std::span(incoming_batch.molecules.data(), incoming_batch.num_ligands)) {
      // Reset the random number generator to improve consistency
      apply<cpu_vectorization::AUTO>(ligand.x(),
                                     ligand.y(),
                                     ligand.z(),
                                     *(best_chromosomes() + index),
                                     num_atoms,
                                     num_rotamers,
                                     ligand.fragments_masks(),
                                     ligand.fragmets_starts(),
                                     ligand.fragments_stops());
      ligand.properties.assign(property_type::SCORE, std::to_string(ligand_scores()[index]));
    }
  }
} // namespace mudock
