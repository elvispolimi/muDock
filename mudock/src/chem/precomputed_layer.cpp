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
    // in caso fosse corretto scambiare memoria per efficienza temporale si potrebbe fare cosi

    // std::vector<fp_type> lig_charges(size_a);
    // std::vector<fp_type> lig_abs_charges(size_a);
    // std::vector<autodock_grid_type> lig_atom_types(size_a);

    // for(std::size_t a = 0; a < size_a; ++a) {
    //     lig_charges[a] = ligand.charge(a);
    //     lig_abs_charges[a] = std::abs(lig_charges[a]);
    //     lig_atom_types[a] = static_cast<autodock_grid_type>(ligand.autodock_type(a));
    // }

    for (std::size_t index_z = 0; index_z < size_z; ++index_z) {
      for (std::size_t index_y = 0; index_y < size_y; ++index_y) {
        for (std::size_t index_x = 0; index_x < size_x; ++index_x) {
           
           const auto grid_elec   = get_eletrostatic().get(index_x, index_y, index_z);
           const auto grid_desolv = get_desolvation().get(index_x, index_y, index_z);

           for(std::size_t index_a = 0; index_a < size_a; ++index_a) {
            //probabilmente conviene salvarsi i valori dell'atomo del ligando a index_a piuttosto che chiederle ripetutamente
            //la scrittura della forse ottimizzazione è scritta sotto

            auto electrostatic_contr = grid_elec * ligand.charge(index_a);
            //auto electrostatic_contr = grid_elec * lig_charges[index_a];

            auto desolvation_contr = grid_desolv * std::abs(ligand.charge(index_a));
            //auto desolvation_contr = grid_desolv * lig_abs_charges[index_a];
            
            //questo si leverebbe
            auto atom_type = ligand.autodock_type(index_a);

            auto atom_contr = get_atom_map(static_cast<autodock_grid_type>(atom_type)).get(index_x, index_y, index_z);
            //auto atom_contr = get_atom_map(lig_atom_types[index_a]).get(index_x, index_y, index_z);
            
            fused_data.get(index_x, index_y, index_z, index_a) = electrostatic_contr + desolvation_contr + atom_contr;
            
          }
        }
      }
    }     
}

} // namespace mudock