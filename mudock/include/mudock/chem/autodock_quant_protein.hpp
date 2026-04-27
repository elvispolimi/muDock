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
    // Metodo di supporto per calcolare la carica da moltiplicare
    fp_type calculate_bin_center(int bin_index) const;

    
    void prepare_fused_maps(const autodock_protein* base_protein);


  public:
            // 16 valori -> 17 mappe con k-Means accuratezza <0.5 tempo_score = 3.301 PICCO_RAM = 417.7MB
            /*
            inline static const std::vector<fp_type> thresholds = {
                -0.81145, -0.69867, -0.60693, -0.44784, -0.24571, -0.07543, 0.03971, 0.11567, 0.21422, 0.32101, 0.38444, 0.42745, 0.47576, 0.53418, 0.62383, 0.77917
            };
            */
            // 15 valori -> 16 mappe con K-Means accuratezza <0.5 tempo_score = 3.2133 PICCO_RAM = 400MB
            /*
            static const std::vector<fp_type> thresholds = {
                -0.81145, -0.69867, -0.60693, -0.44784, -0.24571, -0.07543, 0.03971, 0.11567, 0.21422, 0.32101, 0.38444, 0.42855, 0.49415, 0.60294, 0.77467
            };
            */
            // 14 valori -> 15 mappe con K-Means accuratezza <0.5 tempo_score = 3.2003 PICCO_RAM = 394.7MB
            /*
            static const std::vector<fp_type> thresholds = {
                -0.81145, -0.69867, -0.60693, -0.44784, -0.24571, -0.07543, 0.03971, 0.11567, 0.21468, 0.33126, 0.41632, 0.49156, 0.60294, 0.77467
            };
            */
            //13 valori -> 14 mappe con K-means accuratezza <0.5 tempo_score = 3.1989 PICCO_RAM = 383MB
            /*
            static const std::vector<fp_type> thresholds = {
                -0.81145, -0.69867, -0.60693, -0.44784, -0.24571, -0.07543, 0.03971, 0.11567, 0.21468, 0.33295, 0.43259, 0.55169, 0.75161
            };
            */
            // 12 valori -> 13 mappe con K-means accuratezza <0.5 tempo_score = 3.2002 PICCO_RAM = 371.2MB
            /*
            static const std::vector<fp_type> thresholds = {
                -0.78353, -0.62541, -0.45223, -0.24571, -0.07543, 0.03971, 0.11567, 0.21468, 0.33295, 0.43259, 0.55169, 0.75161
            };
            */
            // 11 valori -> 12 mappe con K-Means accuratezza <0.8 tempo_score = 3.1798 PICCO_RAM = 360MB STOP
            
            inline static const std::vector<fp_type> thresholds = {
                -0.78353, -0.62541, -0.45223, -0.24569, -0.07442, 0.07171, 0.20985, 0.33288, 0.43259, 0.55169, 0.75161
            };

         autodock_quant_protein(const autodock_protein* base_protein) 
    {
        auto sx = base_protein->get_size_x();
        auto sy = base_protein->get_size_xy() / sx;
        auto sz = base_protein->get_size_xyz() / base_protein->get_size_xy();
        
        int num_bins = thresholds.size() + 1;
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