#include "mudock/chem/precomputed_layer.hpp"
#include <mudock/likwid_utils.hpp>
#include "mudock/molecule.hpp"
#include <mudock/chem/autodock_grid_types.hpp>
#include <cmath>    
#include <iostream> 

#include <chrono>  
#include <iostream>

namespace mudock {

void precomputed_protein::prepare_fused_maps(const fp_type* grid_maps, const static_molecule& ligand) {
 
    const auto size_a = ligand.num_atoms(); 

    //dimensione singola griglia
    const std::size_t map_flat_size = sx * sy * sz;

    std::vector<fp_type> lig_charges(size_a);
    std::vector<fp_type> lig_abs_charges(size_a);
    std::vector<int> lig_atom_types(size_a);
    MUDOCK_CPP_MARKER_START("Fase_Setup");
    auto t_start = std::chrono::high_resolution_clock::now();
    for(std::size_t a = 0; a < size_a; ++a) {
        lig_charges[a] = ligand.charge(a);
        lig_abs_charges[a] = std::abs(lig_charges[a]);
        lig_atom_types[a] = static_cast<int>(ligand.autodock_type(a));
    }

    const int ELEC_IDX  = static_cast<int>(autodock_grid_type::ELEC); 
    const int DSOLV_IDX = static_cast<int>(autodock_grid_type::DESOLV);

    //puntatori alle mappe di inizio elec e desolv
    const fp_type* elec_grid  = grid_maps + (ELEC_IDX * map_flat_size);
    const fp_type* dsolv_grid = grid_maps + (DSOLV_IDX * map_flat_size);

    for(std::size_t index_a = 0; index_a < size_a; ++index_a) {
      
      const auto q = lig_charges[index_a];
      const auto abs_q = lig_abs_charges[index_a];
      const int atom_type_idx = lig_atom_types[index_a]; 

      // Puntatore di partenza per la mappa del tipo di atomo corrente
      const fp_type* atom_grid = grid_maps + (atom_type_idx * map_flat_size);
      std::size_t voxel_idx = 0;

      for (std::size_t index_z = 0; index_z < sz; ++index_z) {
        for (std::size_t index_y = 0; index_y < sy; ++index_y) {
          for (std::size_t index_x = 0; index_x < sx; ++index_x) {
             
             const auto grid_elec   = elec_grid[voxel_idx];
             const auto grid_desolv = dsolv_grid[voxel_idx];
             const auto atom_contr  = atom_grid[voxel_idx];
             
             auto electrostatic_contr = grid_elec * q;
             auto desolvation_contr   = grid_desolv * abs_q;
             
             fused_data.get(index_x, index_y, index_z, index_a) = electrostatic_contr + desolvation_contr + atom_contr;
             voxel_idx++;
          }
        }
      }
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> t_diff = t_end - t_start;
    std::cout << "[PROFILAZIONE CUSTOM] Tempo fusione mappa ligando: " << t_diff.count() << " secondi" << std::endl;

    MUDOCK_CPP_MARKER_STOP("Fase_Setup");     
}

} // namespace mudock