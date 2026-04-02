#pragma once

#include <mudock/molecule.hpp>
#include <mudock/grid/mdspan.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {

  struct precomputed_protein{
  
  private:
    md_vector<fp_type, 4> fused_data;

    std::size_t sx, sy, sz;
    
    // La funzione che farà il calcolo
    void prepare_fused_maps(const fp_type* grid_maps, const static_molecule& ligand);

  public:
    precomputed_protein(const fp_type* grid_maps, int size_x, int size_y, int size_z, const static_molecule& _ligand) 
    {
       
        sx = size_x;
        sy = size_y;
        sz = size_z;

        fused_data = md_vector<fp_type, 4>(sx, sy, sz, _ligand.num_atoms());
        
        prepare_fused_maps(grid_maps, _ligand);
    }
    //per il precomputed_adt_score mappa piatta da caricare diretta
    [[nodiscard]] inline const fp_type* get_raw_data() const { 
    return fused_data.data(); 
}
   
    // [[nodiscard]] inline auto get_fused_map(const int atom_index) {
      
    //   return space_grid_view<fp_type>{
    //       p_min,
    //       p_max,
    //       p_center,
    //       p_inv_res,
    //       fused_data.get_slice(
    //           md_index<4>{sx, sy, sz, static_cast<std::size_t>(atom_index)},
    //           md_index<3>{sx, sy, sz}
    //       )
    //   };
    // }
  };

} // namespace mudock