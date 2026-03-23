#pragma once

#include <mudock/chem/autodock_protein.hpp>
#include <mudock/molecule.hpp>
#include <mudock/grid/mdspan.hpp>
#include <mudock/grid/space_grid.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {

  struct precomputed_protein : public autodock_protein {
  
  private:
    md_container<std::vector<fp_type>, 4> fused_data;
    // La funzione che farà il calcolo
    void prepare_fused_maps(const static_molecule& ligand);

  public:
    precomputed_protein(const point3D min, const point3D max, const fp_type resolution, 
                        dynamic_molecule& _molecule, const static_molecule& _ligand)
        : autodock_protein(min, max, resolution, _molecule) {
        
        prepare_fused_maps(_ligand);
    }

    // L'unica differenza con il get precedente è che secondo me
    //  invece di 'autodock_grid_type' è necessario passare l'indice dell'atomo dello specifico ligando
    inline auto get_fused_map(const int atom_index) const {
      return space_grid_view<const fp_type>{
          get_min(),
          get_max(),
          get_center(),
          get_eletrostatic()._inv_resolution, //non essendoci un get l'ho presa da una mappa a caso
          fused_data.get_slice(
              md_index<4>{get_size_x(), get_size_y(), get_size_z(), atom_index},
              md_index<3>{get_size_x(), get_size_y(), get_size_z()}
          )
      };
    }
  };

} // namespace mudock