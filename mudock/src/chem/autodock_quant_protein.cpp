#include "autodock_quant_protein.hpp"
#include <mudock/likwid_utils.hpp>
#include "mudock/molecule.hpp"
#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/chem/autodock_ligand.hpp>
#include <cmath>
#include <iostream>
#include <chrono>

namespace mudock {

void autodock_quant_protein::prepare_fused_maps(const autodock_protein* base_protein) {
    
    
    const std::size_t sx = base_protein->get_size_x();
    const std::size_t size_xy = base_protein->get_size_xy();
    const std::size_t map_flat_size = base_protein->get_map_flat_size();
    const std::size_t sy = size_xy / sx;
    const std::size_t sz = map_flat_size / size_xy;
    
   
    const fp_type* grid_maps = base_protein->get_maps_pointer();

    const int num_bins = get_thresholds.size() + 1;

    
    fp_type* raw_fused_ptr = const_cast<fp_type*>(quantized_fused_maps.data());

    MUDOCK_CPP_MARKER_START("Fase_Setup_Quant");
    auto t_start = std::chrono::high_resolution_clock::now();

    const int ELEC_IDX  = static_cast<int>(autodock_grid_type::ELEC); 
    const int DSOLV_IDX = static_cast<int>(autodock_grid_type::DESOLV);

    const fp_type* elec_grid  = grid_maps + (ELEC_IDX * map_flat_size);
    const fp_type* dsolv_grid = grid_maps + (DSOLV_IDX * map_flat_size);

    for(int bin_idx = 0; bin_idx < num_bins; ++bin_idx) {
        
        // Calcoliamo la carica rappresentativa di questo bin
        const fp_type q = calculate_bin_center(bin_idx);
        const fp_type abs_q = std::abs(q);
        
        for (std::size_t index_z = 0; index_z < sz; ++index_z) {
            for (std::size_t index_y = 0; index_y < sy; ++index_y) {
                for (std::size_t index_x = 0; index_x < sx; ++index_x) {
                    
                    std::size_t voxel_idx = (index_z * sy * sx) + (index_y * sx) + index_x;

                    const auto grid_elec   = elec_grid[voxel_idx];
                    const auto grid_desolv = dsolv_grid[voxel_idx];
                    
                    auto total_val = (grid_elec * q) + (grid_desolv * abs_q);
                    
                    std::size_t final_idx = (bin_idx * map_flat_size) + voxel_idx;    
                    raw_fused_ptr[final_idx] = total_val;
                }
            }
        }
    }

    auto t_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> t_diff = t_end - t_start;
    std::cout << "[PROFILAZIONE] Tempo fusione Ramo 2 (Mappe " << num_bins << "): " << t_diff.count() << "s" << std::endl;

    MUDOCK_CPP_MARKER_STOP("Fase_Setup_Quant");     
}

// Implementazione della funzione di supporto
fp_type autodock_quant_protein::calculate_bin_center(int bin_index) const {
    const auto& thresh = get_thresholds();
    if (bin_index == 0) return thresh[0] - 0.5; 
    if (bin_index >= thresh.size()) return thresh.back() + 0.5; 
    return (thresh[bin_index] + thresh[bin_index - 1]) / 2.0;
}

} // namespace mudock