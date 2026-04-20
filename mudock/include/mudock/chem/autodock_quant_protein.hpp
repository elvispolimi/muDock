#pragma once
#include <mudock/chem/autodock_ligand.hpp>
#include <mudock/molecule.hpp>
#include <mudock/grid/mdspan.hpp>
#include <mudock/type_alias.hpp>
#include <mudock/chem/autodock_protein.hpp> 

#include <vector>
#include <algorithm> 

namespace mudock {

  struct autodock_quant_protein {
  
  private:
    // La nostra mappa fusa. Le 4 dimensioni saranno: [num_bins][sz][sy][sx]
    md_vector<fp_type, 4> quantized_fused_maps;

    // L'array dei confini 
    static const std::vector<fp_type>& get_thresholds(){
        static const std::vector<fp_type> thresh = {
            -1.00000, 
            -0.73010, -0.62000, -0.55120, -0.28820, -0.15000, -0.14350, 
             0.00000, 
             0.08250,  0.14350,  0.15000,  0.16000,  0.28000,  0.37000,  
             0.40000,  0.45000,  0.54380,
             1.00000  
        };
        return thresh;
    }
    
    // Metodo di supporto per calcolare la carica da moltiplicare
    fp_type calculate_bin_center(int bin_index) const;

    
    void prepare_fused_maps(const autodock_protein* base_protein);

  public:
         autodock_quant_protein(const autodock_protein* base_protein) 
    {
    
        auto sx = base_protein->get_size_x();
        auto sy = base_protein->get_size_xy() / sx;
        auto sz = base_protein->get_size_xyz() / base_protein->get_size_xy();

        int num_bins = get_thresholds.size() + 1;

        quantized_fused_maps = md_vector<fp_type, 4>(num_bins, sz, sy, sx);
        prepare_fused_maps(base_protein);
    }
    

    
    [[nodiscard]] inline const fp_type* get_raw_data() const { 
        return quantized_fused_maps.data(); 
    }
  };

  // Da richiamare nel file adt_score.hpp quando si setta il batch dei ligandi
  // associa ad ogni atom_index del ligando il corrispettivo bin (l'indice) in modo da accedervi in O(1) nello score
  inline std::vector<int> build_atom_to_bin_map(const static_molecule& ligand, const std::vector<fp_type>& thresh) {
      std::vector<int> atom_bins(ligand.num_atoms());
      for (int i = 0; i < ligand.num_atoms(); ++i) {
          fp_type q = ligand.charge(i);
          // Ricerca binaria
          auto it = std::upper_bound(thresh.begin(), thresh.end(), q);
          atom_bins[i] = std::distance(thresh.begin(), it);
      }
      return atom_bins;
  }

} // namespace mudock