#include <cstddef>
#include <cstring>
#include <hip/hip_runtime.h>
#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/chem/autodock_ligand.hpp>
#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/cpp_implementation/center_of_mass.hpp>
#include <mudock/cpp_implementation/geometric_transformations.hpp>
#include <mudock/cpp_implementation/mutate.hpp>
#include <mudock/cpp_implementation/weed_bonds.hpp>
#include <mudock/grid.hpp>
#include <mudock/hip_implementation/evaluate_fitness.hpp>
#include <mudock/hip_implementation/virtual_screen.hpp>
#include <mudock/utils.hpp>

namespace mudock {
  // TODO add to bucketizer
  static constexpr std::size_t max_non_bonds{1 << 26};

  virtual_screen_hip::virtual_screen_hip(const knobs k, const std::shared_ptr<const device> dev)
      : configuration(k),
        dev(dev),
        stream(dev->get_stream()),
        original_ligand_x(stream),
        original_ligand_y(stream),
        original_ligand_z(stream),
        scratch_ligand_x(stream),
        scratch_ligand_y(stream),
        scratch_ligand_z(stream),
        ligand_vol(stream),
        ligand_solpar(stream),
        ligand_charge(stream),
        ligand_num_hbond(stream),
        ligand_num_atoms(stream),
        ligand_num_rotamers(stream),
        ligand_fragments(stream),
        frag_start_atom_indices(stream),
        frag_stop_atom_indices(stream),
        index_nonbonds(stream),
        nonbond_a1(stream),
        nonbond_a2(stream),
        nonbond_cA(stream),
        nonbond_cB(stream),
        nonbond_xB(stream),
        map_texture_index(stream),
        ligand_scores(stream),
        chromosomes(stream),
        best_chromosomes(stream),
        hiprand_states(stream) {}

  void virtual_screen_hip::operator()(batch &incoming_batch) {
    const auto wavefront_size = dev->get_wavefront();

    const std::size_t batch_atoms    = incoming_batch.batch_max_atoms;
    const std::size_t batch_rotamers = incoming_batch.batch_max_rotamers;
    const std::size_t batch_ligands  = incoming_batch.num_ligands;
    // Resize data structures
    const std::size_t tot_atoms_in_batch = batch_ligands * batch_atoms;
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
    frag_start_atom_indices.alloc(tot_rotamers_in_batch);
    frag_stop_atom_indices.alloc(tot_rotamers_in_batch);
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
      const int stride_masks                    = index * batch_rotamers * batch_atoms;
      const int stride_rotamers                 = index * batch_rotamers;
      assert(batch_rotamers > ligand.get()->num_rotamers());

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
                  adt_ligand.get_atom_map_index(),
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
    ligand_num_hbond.copy_host2device();
    map_texture_index.copy_host2device();
    index_nonbonds.copy_host2device();
    nonbond_a1.copy_host2device();
    nonbond_a2.copy_host2device();
    nonbond_cA.copy_host2device();
    nonbond_cB.copy_host2device();
    nonbond_xB.copy_host2device();

    // Setup hip random
    // A state for each thread
    const int num_threads = batch_ligands * wavefront_size;
    hiprand_states.alloc(num_threads);

    // Simulate the population evolution for the given amount of time
    const auto num_generations = configuration.num_generations;
    // The shared memory contains:
    // - each chromosome's score at last population evaluation
    // Max due to the reduction at the end to find the highest scores per each chromosome
    const std::size_t min_energy_reduction_s_mem =
        std::max(configuration.population_number, static_cast<std::size_t>(wavefront_size)) * sizeof(fp_type);
    const std::size_t shared_mem = min_energy_reduction_s_mem;

    constexpr_for<0, reorder_buffer::atoms_clusters.size(), 1>([&](const auto atoms_index) {
      const auto n_atoms = reorder_buffer::atoms_clusters[atoms_index];
      constexpr_for<0, reorder_buffer::rotamer_clusters.size(), 1>([&](const auto rotamers_index) {
        const auto n_rotamers = reorder_buffer::rotamer_clusters[rotamers_index];
        if (batch_atoms == n_atoms && batch_rotamers == n_rotamers)
          evaluate_fitness<n_atoms, n_rotamers>
              <<<batch_ligands, wavefront_size, shared_mem, stream>>>(num_generations,
                                                                      configuration.tournament_length,
                                                                      configuration.mutation_prob,
                                                                      configuration.population_number,
                                                                      population_stride,
                                                                      batch_atoms,
                                                                      batch_rotamers,
                                                                      max_non_bonds,
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
                                                                      frag_start_atom_indices.dev_pointer(),
                                                                      frag_stop_atom_indices.dev_pointer(),
                                                                      chromosomes.dev_pointer(),
                                                                      dev.get()->atom_tex.dev_pointer(),
                                                                      map_texture_index.dev_pointer(),
                                                                      hiprand_states.dev_pointer(),
                                                                      ligand_scores.dev_pointer(),
                                                                      best_chromosomes.dev_pointer());
      });
    });

    MUDOCK_CHECK_KERNELCALL();
    MUDOCK_CHECK(hipStreamSynchronize(stream));

    // Copy back chromosomes and scores
    // FIXME move before the synchronize
    best_chromosomes.copy_device2host();
    ligand_scores.copy_device2host();

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
