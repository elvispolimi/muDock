#pragma once

#include <mudock/chem/autodock_protein.hpp>
#include <mudock/molecule.hpp>
#include <mudock/grid/mdspan.hpp>
#include <mudock/grid/space_grid.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {

  struct precomputed_protein : public autodock_protein {
  
  private:
    md_vector<fp_type, 4> fused_data;
    
    // La funzione che farà il calcolo
    void prepare_fused_maps(const static_molecule& ligand);

  public:
    precomputed_protein(const point3D min, 
                        const point3D max, 
                        const fp_type resolution, 
                        dynamic_molecule& _molecule, 
                        const static_molecule& _ligand)
        : autodock_protein(min, max, resolution, _molecule) 
    {
       
        const auto sx = get_size_x();
        const auto sy = get_size_xy() / sx;
        const auto sz = get_size_xyz() / get_size_xy();

       
        fused_data = md_vector<fp_type, 4>(sx, sy, sz, _ligand.num_atoms());
        
        prepare_fused_maps(_ligand);
    }

   
    [[nodiscard]] inline auto get_fused_map(const int atom_index) {
      
      const auto sx = get_size_x();
      const auto sy = get_size_xy() / sx;
      const auto sz = get_size_xyz() / get_size_xy();

      
      return space_grid_view<fp_type>{
          get_min(),
          get_max(),
          get_center(),
          get_eletrostatic()._inv_resolution,
          fused_data.get_slice(
              md_index<4>{sx, sy, sz, static_cast<std::size_t>(atom_index)},
              md_index<3>{sx, sy, sz}
          )
      };
    }
  };

} // namespace mudock