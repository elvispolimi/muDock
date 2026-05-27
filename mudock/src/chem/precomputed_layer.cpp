#include "mudock/chem/precomputed_layer.hpp"
#include <mudock/likwid_utils.hpp>
#include "mudock/molecule.hpp"
#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/chem/autodock_ligand.hpp>
#include <cmath>    
#include <iostream> 
 

namespace mudock {

void precomputed_protein::prepare_fused_maps(const fp_type* grid_maps, 
                                            const autodock_ligand& adt_ligand, 
                                            const static_molecule& ligand) {
    
    const auto size_a = adt_ligand.num_atoms(); 
    const std::size_t map_flat_size = sx * sy * sz;
    
    // Puntatore grezzo per la scrittura diretta
    fp_type* raw_fused_ptr = const_cast<fp_type*>(fused_data.data());

    MUDOCK_CPP_MARKER_START("Fase_Setup");
    // Indici per le mappe globali ELEC e DESOLV
    const int ELEC_IDX  = static_cast<int>(autodock_grid_type::ELEC); 
    const int DSOLV_IDX = static_cast<int>(autodock_grid_type::DESOLV);

    const fp_type* elec_grid  = grid_maps + (ELEC_IDX * map_flat_size);
    const fp_type* dsolv_grid = grid_maps + (DSOLV_IDX * map_flat_size);

   

    for(std::size_t index_a = 0; index_a < size_a; ++index_a) {
        
        const fp_type q = ligand.charge(index_a);
        const fp_type abs_q = std::abs(q);
        
        const auto ff_type = ligand.autodock_type(index_a);
        const autodock_grid_type grid_type = autodock_grid_from_ff(ff_type);
        const std::size_t grid_idx = static_cast<std::size_t>(grid_type);
        const fp_type* atom_grid = grid_maps + (grid_idx * map_flat_size);

        
        for (std::size_t index_z = 0; index_z < sz; ++index_z) {
            for (std::size_t index_y = 0; index_y < sy; ++index_y) {
                for (std::size_t index_x = 0; index_x < sx; ++index_x) {
                    
                    // Indice 1D della griglia originale [Z][Y][X]
                    std::size_t voxel_idx = (index_z * sy * sx) + (index_y * sx) + index_x;

                    const auto grid_elec   = elec_grid[voxel_idx];
                    const auto grid_desolv = dsolv_grid[voxel_idx];
                    const auto atom_contr  = atom_grid[voxel_idx];
                    
                    auto total_val = atom_contr + (grid_elec * q) + (grid_desolv * abs_q);
                    
                    
                    // Scrittura nell'ordine [Atomo][Z][Y][X]
                    std::size_t final_idx = (index_a * map_flat_size) + voxel_idx;    
                    raw_fused_ptr[final_idx] = total_val;
                }
            }
        }
    }

    MUDOCK_CPP_MARKER_STOP("Fase_Setup");     
}

} // namespace mudock