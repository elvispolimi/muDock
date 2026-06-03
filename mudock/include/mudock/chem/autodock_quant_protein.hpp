#pragma once
#include <mudock/chem/autodock_ligand.hpp>
#include <mudock/molecule.hpp>
#include <mudock/grid/mdspan.hpp>
#include <mudock/type_alias.hpp>
#include <mudock/chem/autodock_protein.hpp> 
#include <mudock/likwid_utils.hpp>

#include <vector>
#include <algorithm> 

namespace mudock {

  struct autodock_quant_protein {
  
  private:
    // fused quantized maps: [bin_index][z][y][x]
    md_vector<fp_type, 4> quantized_fused_maps;
    // helper function to calculate the center of a bin given its index
    fp_type calculate_bin_center(int bin_index) const;

    
    void prepare_fused_maps(const autodock_protein* base_protein);


  public:
            // 16 values -> 17 maps 
            /*
            inline static constexpr std::array<fp_type, 16> thresholds = {
                -0.81145f, -0.69867f, -0.60693f, -0.44784f, -0.24571f, -0.07543f, 0.03971f, 0.11567f, 0.21422f, 0.32101f, 0.38444f, 0.42745f, 0.47576f, 0.53418f, 0.62383f, 0.77917f
            };
            */
            // 15 values -> 16 maps
            /*
            inline static constexpr std::array<fp_type, 15> thresholds = {
                -0.81145f, -0.69867f, -0.60693f, -0.44784f, -0.24571f, -0.07543f, 0.03971f, 0.11567f, 0.21422f, 0.32101f, 0.38444f, 0.42855f, 0.49415f, 0.60294f, 0.77467f
            };
            */
            // 14 values -> 15 maps 
            /*
            inline static constexpr std::array<fp_type, 14> thresholds = {
                -0.81145f, -0.69867f, -0.60693f, -0.44784f, -0.24571f, -0.07543f, 0.03971f, 0.11567f, 0.21468f, 0.33126f, 0.41632f, 0.49156f, 0.60294f, 0.77467f
            };
            */
            //13 values -> 14 maps
            
            inline static constexpr std::array<fp_type, 13> thresholds = {
                -0.81145f, -0.69867f, -0.60693f, -0.44784f, -0.24571f, -0.07543f, 0.03971f, 0.11567f, 0.21468f, 0.33295f, 0.43259f, 0.55169f, 0.75161f
            };
            
            // 12 values -> 13 maps 
            /*
            inline static constexpr std::array<fp_type, 12> thresholds = {
                -0.78353f, -0.62541f, -0.45223f, -0.24571f, -0.07543f, 0.03971f, 0.11567f, 0.21468f, 0.33295f, 0.43259f, 0.55169f, 0.75161f
            };
            */
            // 11 values -> 12 maps
            /*
            inline static constexpr std::array<fp_type, 11> thresholds = {
                -0.78353f, -0.62541f, -0.45223f, -0.24569f, -0.07442f, 0.07171f, 0.20985f, 0.33288f, 0.43259f, 0.55169f, 0.75161f
            };
            */
            // 7 values -> 8 maps
            /*
            inline static constexpr std::array<fp_type, 7> thresholds = {
                -0.77754f, -0.59979f, -0.34481f, -0.08214f, 0.07318f, 0.26056f, 0.49312f
            };*/
         autodock_quant_protein(const autodock_protein* base_protein) 
    {
        auto sx = base_protein->get_size_x();
        auto sy = base_protein->get_size_xy() / sx;
        auto sz = base_protein->get_size_xyz() / base_protein->get_size_xy();
        
        int num_bins = static_cast<int>(thresholds.size()) + 1;
        quantized_fused_maps = md_vector<fp_type, 4>(num_bins, sz, sy, sx);
        prepare_fused_maps(base_protein);
        
    }
    
    [[nodiscard]] inline const fp_type* get_raw_data() const { 
        return quantized_fused_maps.data(); 
    }
  };

  /*This function will create a mapping from each ligand atom index to the corresponding bin index based on the charge thresholds.
   This allows for O(1) access during scoring, as we can directly retrieve the bin index for each atom without performing a search.*/
  inline std::vector<int> build_atom_to_bin_map(const static_molecule& ligand, const std::array<fp_type, 13>& thresh) {
      std::vector<int> atom_bins(ligand.num_atoms());
      for (int i = 0; i < ligand.num_atoms(); ++i) {
          fp_type q = ligand.charge(i);
          // Binary search to find the appropriate bin for the charge q
          auto it = std::upper_bound(thresh.begin(), thresh.end(), q);
          atom_bins[i] = static_cast<int>(std::distance(thresh.begin(), it));
      }
      return atom_bins;
  }

} // namespace mudock