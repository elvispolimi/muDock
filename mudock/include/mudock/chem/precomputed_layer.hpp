#pragma once
#include <mudock/chem/autodock_ligand.hpp>
#include <mudock/molecule.hpp>
#include <mudock/grid/mdspan.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {

  struct precomputed_protein{
  
  private:
    md_vector<fp_type, 4> fused_data;

    std::size_t sx, sy, sz;
    
    // La funzione che farà il calcolo
    void prepare_fused_maps(const fp_type* grid_maps, const autodock_ligand& adt_ligand, const static_molecule& ligand);

  public:
    precomputed_protein(const fp_type* grid_maps, int size_x, int size_y, int size_z, const autodock_ligand& adt_ligand, const static_molecule& ligand) 
    {
       
        sx = size_x;
        sy = size_y;
        sz = size_z;

        fused_data = md_vector<fp_type, 4>(adt_ligand.num_atoms(), sz, sy, sx);
        
        prepare_fused_maps(grid_maps, adt_ligand, ligand);
    }
    
    //per il precomputed_adt_score mappa piatta da caricare diretta
    [[nodiscard]] inline const fp_type* get_raw_data() const { 
    return fused_data.data(); 
}
   
  };

} // namespace mudock