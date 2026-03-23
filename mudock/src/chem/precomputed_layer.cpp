#include "mudock/chem/precomputed_layer.hpp"
#include "mudock/molecule.hpp"
#include <cmath>    
#include <iostream> 

namespace mudock {

void precomputed_protein::prepare_fused_maps(const static_molecule& ligand) {

    const auto size_x = get_size_x();
    const auto size_y = get_size_xy() / size_x;
    const auto size_z = get_size_xyz() / get_size_xy();
    
    const auto size_a = ligand.num_atoms();  

    for (std::size_t index_z = 0; index_z < size_z; ++index_z) {
      for (std::size_t index_y = 0; index_y < size_y; ++index_y) {
        for (std::size_t index_x = 0; index_x < size_x; ++index_x) {
           
           const auto grid_elec   = get_eletrostatic().get(index_x, index_y, index_z);
           const auto grid_desolv = get_desolvation().get(index_x, index_y, index_z);

           for(std::size_t index_a = 0; index_a < size_a; ++index_a) {
            
            auto electrostatic_contr = grid_elec * ligand.charge(index_a);
            
            auto desolvation_contr = grid_desolv * std::abs(ligand.charge(index_a));
            
            auto atom_type = ligand.autodock_type(index_a);

            auto atom_contr = get_atom_map(static_cast<autodock_grid_type>(atom_type)).get(index_x, index_y, index_z);
            
            fused_data.get(index_x, index_y, index_z, index_a) = electrostatic_contr + desolvation_contr + atom_contr;
            
          }
        }
      }
    }     
}

} // namespace mudock