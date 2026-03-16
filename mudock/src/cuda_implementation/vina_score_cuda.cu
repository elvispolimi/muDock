#pragma once

#include <cstdio>
#include <mudock/compute/devices_memory.hpp>
#include <mudock/compute/vina_score_kernel.hpp>
#include <mudock/compute/devices_memory.hpp>
#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/cuda_implementation/vina_score_cuda.cuh>
#include <mudock/cuda_implementation/cuda_texture.cuh>
#include <mudock/cuda_implementation/cuda_utils.cuh>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>
#include <mudock/utils.hpp>
#include <stdexcept>

#include <mudock/type_alias.hpp>

/// -10.219824: risultato finale variazione solo di +0.000009       <- with parallelisation
/// Score inter -15.313681, Score intra -1.478089, Score -10.219815 <- with no parallelisation
/// Score inter -15.313679, Score intra -1.478089, Score -10.219813 <- real

#define BUCKET_MULTIPLIER 3

namespace mudock {

    __device__ static constexpr fp_type GAUSS1_COEFF_CUDA{- 0.035579};
    __device__ static constexpr fp_type GAUSS2_COEFF_CUDA{- 0.005156};
    __device__ static constexpr fp_type REPULSION_COEFF_CUDA{0.840245};
    __device__ static constexpr fp_type HYDROPHOBIC_COEFF_CUDA{- 0.035069};
    __device__ static constexpr fp_type H_BOND_COEFF_CUDA{- 0.587439};
    __device__ static constexpr fp_type NROT_COEFF_CUDA{0.05846};

    __device__ inline fp_type distance(fp_type x, fp_type y, fp_type z) {
        return sqrt( x*x + y*y + z*z );
    }

    __device__ inline fp_type gauss1(const size_t idx, const fp_type* __restrict__ dst_mtx) {
        fp_type gauss1 = 0;
        if(dst_mtx[idx] != 0) gauss1 = exp(- pow(dst_mtx[idx] / 0.5, 2));
        return gauss1;
    }

    __device__ inline fp_type gauss2(const size_t idx, const fp_type* __restrict__ dst_mtx) {
        fp_type gauss2 = 0;
        if(dst_mtx[idx] != 0) gauss2 = exp(- pow((dst_mtx[idx] - 3) / 2, 2));
        return gauss2;
    }

    __device__ inline fp_type repulsion(const size_t idx, const fp_type* __restrict__ dst_mtx) {
        return pow((dst_mtx[idx] < 0) * dst_mtx[idx], 2);
    }

    __device__ inline fp_type hydrophobic(const size_t idx, const fp_type* __restrict__ dst_mtx, const int* __restrict__ rec_lig_is_hydrophobic) {
        if(rec_lig_is_hydrophobic[idx] < 0) return 0; // Sentinel value, no interaction
        bool hydro_1 = rec_lig_is_hydrophobic[idx] && (dst_mtx[idx] <= 0.5);
        bool hydro_2_cond = rec_lig_is_hydrophobic[idx] && (dst_mtx[idx] > 0.5) && (dst_mtx[idx] < 1.5);
        fp_type hydro_2 = 1.5 * hydro_2_cond - hydro_2_cond * dst_mtx[idx];
        return hydro_1 + hydro_2;
    }

    __device__ inline fp_type hbonding(const size_t idx, const fp_type* __restrict__ dst_mtx, const int* __restrict__ rec_lig_is_hb) {
        if(rec_lig_is_hb[idx] < 0) return 0; // Sentinel value, no interaction
        bool h_bond_1 = rec_lig_is_hb[idx] && (dst_mtx[idx] <= -0.7);
        bool h_bond_2_cond = rec_lig_is_hb[idx] && (dst_mtx[idx] < 0) && (dst_mtx[idx] > -0.7);
        fp_type h_bond_2 = h_bond_2_cond * (- dst_mtx[idx]) / 0.7;
        return h_bond_1 + h_bond_2;
    }


    __device__ inline fp_type score_function(
        const fp_type* __restrict__ dst_mtx, 
        const int* __restrict__ ij_is_hydrophobic,
        const int* __restrict__ ij_is_hbond,
        const size_t size
    ) {

        fp_type g1 = 0; 
        fp_type g2 = 0; 
        fp_type rep = 0;
        fp_type hydro = 0;
        fp_type hbond = 0;

        for(size_t i = threadIdx.x; i < size; i += blockDim.x) {
            if(ij_is_hydrophobic[i] < 0) continue; // Sentinel value, no interaction
            g1  += gauss1(i, dst_mtx);
            g2  += gauss2(i, dst_mtx);
            rep += repulsion(i, dst_mtx);
            hydro += hydrophobic(i, dst_mtx, ij_is_hydrophobic);
            hbond += hbonding(i, dst_mtx, ij_is_hbond);
        }

        return GAUSS1_COEFF_CUDA * g1 + GAUSS2_COEFF_CUDA * g2 + REPULSION_COEFF_CUDA * rep + HYDROPHOBIC_COEFF_CUDA * hydro + H_BOND_COEFF_CUDA * hbond;
    }

    __device__ inline void parse_data(
        /// Protein data
        const size_t num_atoms_protein,
        const fp_type* __restrict__ protein_x,
        const fp_type* __restrict__ protein_y,
        const fp_type* __restrict__ protein_z,
        const int* __restrict__ p_is_hbond_acceptor,
        const int* __restrict__ p_is_hbond_donor,
        const int* __restrict__ p_is_hydrophobic,
        const fp_type* __restrict__ p_vdw_radius,

        ///Ligand data
        const size_t num_atoms_ligand,
        const fp_type* __restrict__ ligand_x,
        const fp_type* __restrict__ ligand_y,
        const fp_type* __restrict__ ligand_z,
        const int* __restrict__ l_is_hbond_acceptor,
        const int* __restrict__ l_is_hbond_donor,
        const int* __restrict__ l_is_hydrophobic,
        const fp_type* __restrict__ l_vdw_radius,

        /// Output buffers
        fp_type* __restrict__ dst_mtx,
        int* __restrict__ is_hbond,
        int* __restrict__ is_hydrophobic
    ){

        fp_type vdw_sum;
        for (size_t proteinIdx = 0; proteinIdx < num_atoms_protein; proteinIdx++) {
            for(size_t ligandIdx = threadIdx.x; ligandIdx < num_atoms_ligand; ligandIdx += blockDim.x) {
                fp_type dst = distance(
                    (protein_x[proteinIdx] - ligand_x[ligandIdx]),
                    (protein_y[proteinIdx] - ligand_y[ligandIdx]),
                    (protein_z[proteinIdx] - ligand_z[ligandIdx])
                );

                size_t idx = ligandIdx + (proteinIdx * num_atoms_ligand);
                
                if(dst > 8){
                    is_hydrophobic[idx] = -1; /// Sentinel value. No interaction.
                }
                else{
                    vdw_sum = p_vdw_radius[proteinIdx] + l_vdw_radius[ligandIdx];
                    dst_mtx[idx] = dst - vdw_sum;
                    is_hbond[idx] = (p_is_hbond_acceptor[proteinIdx] && l_is_hbond_donor[ligandIdx]) || (l_is_hbond_acceptor[ligandIdx] && p_is_hbond_donor[proteinIdx]);
                    is_hydrophobic[idx] = p_is_hydrophobic[proteinIdx] && l_is_hydrophobic[ligandIdx];
                }

            }
        }
    }

    __device__ inline void parse_intra_data(
        const fp_type* __restrict__ ligand_x,
        const fp_type* __restrict__ ligand_y,
        const fp_type* __restrict__ ligand_z,
        const int* __restrict__ l_is_hbond_acceptor,
        const int* __restrict__ l_is_hbond_donor,
        const int* __restrict__ l_is_hydrophobic,
        const fp_type* __restrict__ l_vdw_radius,
        const int* __restrict__ interacting_pairs_first,
        const int* __restrict__ interacting_pairs_second,
        const size_t num_interacting_pairs,

        /// Output buffers
        fp_type* __restrict__ intra_dst_mtx,
        int* __restrict__ intra_is_hbond,
        int* __restrict__ intra_is_hydrophobic
    ){

        fp_type vdw_sum;
        for(size_t i = threadIdx.x; i < num_interacting_pairs; i += blockDim.x) {
            int atom_1 = interacting_pairs_first[i];
            int atom_2 = interacting_pairs_second[i];

            fp_type dst = distance(
                (ligand_x[atom_1] - ligand_x[atom_2]),
                (ligand_y[atom_1] - ligand_y[atom_2]),
                (ligand_z[atom_1] - ligand_z[atom_2])
            );
            
            if(dst > 8){
                intra_is_hydrophobic[i] = -1; /// Sentinel value. No interaction.
            }
            else {
                vdw_sum = l_vdw_radius[atom_1] + l_vdw_radius[atom_2];
                intra_dst_mtx[i] = dst - vdw_sum;
                intra_is_hbond[i] = (l_is_hbond_acceptor[atom_1] && l_is_hbond_donor[atom_2]) || (l_is_hbond_acceptor[atom_2] && l_is_hbond_donor[atom_1]);
                intra_is_hydrophobic[i] = l_is_hydrophobic[atom_1] && l_is_hydrophobic[atom_2];
            }
        }
    }

    __device__ inline fp_type scoring_cuda(  
        /// Protein data
        const size_t num_atoms_protein,
        const fp_type* __restrict__ protein_x,
        const fp_type* __restrict__ protein_y,
        const fp_type* __restrict__ protein_z,
        const int* __restrict__ p_is_hbond_acceptor,
        const int* __restrict__ p_is_hbond_donor,
        const int* __restrict__ p_is_hydrophobic,
        const fp_type* __restrict__ p_vdw_radius,

        ///Ligand data
        const size_t num_atoms_ligand,
        const fp_type* __restrict__ ligand_x,
        const fp_type* __restrict__ ligand_y,
        const fp_type* __restrict__ ligand_z,
        const int* __restrict__ l_is_hbond_acceptor,
        const int* __restrict__ l_is_hbond_donor,
        const int* __restrict__ l_is_hydrophobic,
        const fp_type* __restrict__ l_vdw_radius,
        const size_t active_torsions,
        const int* __restrict__ interacting_pairs_first,
        const int* __restrict__ interacting_pairs_second,
        const size_t num_interacting_pairs,
        
        /// Buffers (worst dim: num_atoms_ligand * num_atoms_protein)
        fp_type* __restrict__ dst_mtx,
        int* __restrict__ is_hbond,
        int* __restrict__ is_hydrophobic

    ){
        /// TODO: essere sicuri che tutti gli elementi siano diversi dall'idrogeno
        size_t mtx_size = 0;

        parse_data(
            num_atoms_protein,
            protein_x,
            protein_y,
            protein_z,
            p_is_hbond_acceptor,
            p_is_hbond_donor,
            p_is_hydrophobic,
            p_vdw_radius,
            
            num_atoms_ligand,
            ligand_x,
            ligand_y,
            ligand_z,
            l_is_hbond_acceptor,
            l_is_hbond_donor,
            l_is_hydrophobic,
            l_vdw_radius,
            
            dst_mtx, 
            is_hbond, 
            is_hydrophobic
        );
        
        mtx_size = num_atoms_ligand * num_atoms_protein;
        fp_type inter_score = score_function(dst_mtx, is_hydrophobic, is_hbond, mtx_size);
        
        parse_intra_data(
            ligand_x,
            ligand_y,
            ligand_z,
            l_is_hbond_acceptor,
            l_is_hbond_donor,
            l_is_hydrophobic,
            l_vdw_radius,
            interacting_pairs_first,
            interacting_pairs_second,
            num_interacting_pairs,
            
            dst_mtx,
            is_hbond, 
            is_hydrophobic
        );

        mtx_size = num_interacting_pairs;
        fp_type intra_score = score_function(dst_mtx, is_hydrophobic, is_hbond, mtx_size);
        
        fp_type score = (inter_score + intra_score) / ( 1 + NROT_COEFF_CUDA * active_torsions);
        
        // printf("dst_mtx len %ld\n", dst_mtx.size());
        // printf("intra_dst_mtx len %ld\n", intra_dst_mtx.size());
        // printf("Score inter %f, Score intra %f, Score %f\n", inter_score, intra_score, score); 
        
        return score;
    }

  __global__ void calc_energy(const int atom_stride,
                              const int scores_per_ligand,
                              const int *__restrict__ num_atoms_b, 
                              const fp_type *__restrict__ x_scratch,
                              const fp_type *__restrict__ y_scratch,
                              const fp_type *__restrict__ z_scratch,
                              const int* __restrict__ ligand_is_hbond_acceptor,
                              const int* __restrict__ ligand_is_hbond_donor,
                              const int* __restrict__ ligand_is_hydrophobic,
                              const fp_type* __restrict__ ligand_vdw_radius,
                              const int *__restrict__ ligand_active_torsions,
                              const int* __restrict__ ligand_interacting_pairs_first, 
                              const int* __restrict__ ligand_interacting_pairs_second,
                              const int* __restrict__ ligand_num_interacting_pairs, 
                              const int* __restrict__ ligand_interacting_pairs_offset, 
                              const size_t num_atoms_protein,
                              const fp_type* __restrict__ protein_x,
                              const fp_type* __restrict__ protein_y,
                              const fp_type* __restrict__ protein_z,
                              const int* __restrict__ p_is_hbond_acceptor,
                              const int* __restrict__ p_is_hbond_donor,
                              const int* __restrict__ p_is_hydrophobic,
                              const fp_type* __restrict__ p_vdw_radius,
                              fp_type *__restrict__ scores,
                              
                              fp_type* __restrict__ dst_mtx_b,
                              int* __restrict__ is_hbond_b,
                              int* __restrict__ is_hydrophobic_b
                              ) {

    const int ligand_id       = blockIdx.x;
    const int local_thread_id = threadIdx.x;

    const int num_atoms_ligand    = num_atoms_b[ligand_id];
    const int active_torsions = ligand_active_torsions[ligand_id];
    const int num_interacting_pairs = ligand_num_interacting_pairs[ligand_id];
    const int interacting_pairs_offset = ligand_interacting_pairs_offset[ligand_id];
    const int stride       = ligand_id * atom_stride;
    
    fp_type* dst_mtx = dst_mtx_b + num_atoms_ligand * num_atoms_protein * ligand_id;
    int* is_hbond = is_hbond_b + num_atoms_ligand * num_atoms_protein * ligand_id;
    int* is_hydrophobic = is_hydrophobic_b + num_atoms_ligand * num_atoms_protein * ligand_id;
    
    const fp_type* l_scratch_x = x_scratch + stride * scores_per_ligand;
    const fp_type* l_scratch_y = y_scratch + stride * scores_per_ligand;
    const fp_type* l_scratch_z = z_scratch + stride * scores_per_ligand;
    const int* l_is_hbond_acceptor = ligand_is_hbond_acceptor + stride * scores_per_ligand;
    const int* l_is_hbond_donor = ligand_is_hbond_donor + stride * scores_per_ligand;
    const int* l_is_hydrophobic = ligand_is_hydrophobic + stride * scores_per_ligand;
    const fp_type* l_vdw_radius = ligand_vdw_radius + stride * scores_per_ligand;
    /// TODO: check this cus i am not really sure about his behaviour in case of multiple scores per ligand
    const int* interacting_pairs_first = ligand_interacting_pairs_first + interacting_pairs_offset * scores_per_ligand;
    const int* interacting_pairs_second = ligand_interacting_pairs_second + interacting_pairs_offset * scores_per_ligand;

    fp_type* scores_l = scores + ligand_id * scores_per_ligand;

    for (int scores_index = 0; scores_index < scores_per_ligand; ++scores_index) {
      // Copy original coordinates
      const fp_type* ligand_x = l_scratch_x + scores_index * atom_stride;
      const fp_type* ligand_y = l_scratch_y + scores_index * atom_stride;
      const fp_type* ligand_z = l_scratch_z + scores_index * atom_stride;

      // Calculate energy 
      fp_type result = scoring_cuda(num_atoms_protein, 
          protein_x, 
          protein_y,
          protein_z, 
          p_is_hbond_acceptor, 
          p_is_hbond_donor, 
          p_is_hydrophobic, 
          p_vdw_radius, 
          num_atoms_ligand, 
          ligand_x, 
          ligand_y, 
          ligand_z, 
          l_is_hbond_acceptor,
          l_is_hbond_donor,
          l_is_hydrophobic,
          l_vdw_radius,
          active_torsions,
          interacting_pairs_first, 
          interacting_pairs_second, 
          num_interacting_pairs,
          dst_mtx,
          is_hbond,
          is_hydrophobic
          );   


#ifdef MUDOCK_TEST
#else
#endif

#pragma unroll
      for (int offset = BLOCK_SIZE / 2; offset > 0; offset /= 2) {
        result += __shfl_down_sync(0xffffffff, result, offset);
      }

      if (local_thread_id == 0) {
        scores_l[scores_index]         = result;
      }
    }
  }

  template<>
  void vina_score_kernel<queue_cuda>::operator()() {
    const int dev_id = q->get_id();
    

    void* args[] = {(void*) &batch_atoms,
                    (void*) &scores_per_ligand,
                    (void*) &num_ligand_atoms_b,
                    (void*) &x_scratch_b,
                    (void*) &y_scratch_b,
                    (void*) &z_scratch_b,
                    (void*) &l_is_hbond_acceptor_b,
                    (void*) &l_is_hbond_donor_b,
                    (void*) &l_is_hydrophobic_b,
                    (void*) &l_vdw_radius_b,
                    (void*) &active_torsions_b,
                    (void*) &interacting_pairs_first_b,
                    (void*) &interacting_pairs_second_b,
                    (void*) &num_interacting_pairs_b,
                    (void*) &interacting_pairs_offset_b,
                    (void*) &num_atoms_protein,
                    (void*) &protein_x,
                    (void*) &protein_y,
                    (void*) &protein_z,
                    (void*) &p_is_hbond_acceptor,
                    (void*) &p_is_hbond_donor,
                    (void*) &p_is_hydrophobic,
                    (void*) &p_vdw_radius,
                    (void*) &scores_b,
                    (void*) &dst_mtx_b,
                    (void*) &is_hbond_b,
                    (void*) &is_hydrophobic_b
    };
    constexpr_switch_bucket<0, reorder_buffer<static_molecule>::get_num_atom_clusters(), 1>(
        [&](const auto atom_index) {
          const auto max_atoms = reorder_buffer<static_molecule>::atoms_clusters[atom_index];
          q->launch_kernel((void*) calc_energy, args, batch_ligands, BLOCK_SIZE);
        },
        batch_atoms,
        reorder_buffer<static_molecule>::atoms_clusters.data());

  }; // namespace mudock

  template<int MAX_ATOMS>
  int get_evaluate_fitness_batch(const int device_id) {
    MUDOCK_CHECK(cudaSetDevice(device_id));
    cudaDeviceProp props;
    MUDOCK_CHECK(cudaGetDeviceProperties(&props, device_id));
    int num_block_per_SM = 0;
    // TODO
    MUDOCK_CHECK(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&num_block_per_SM,
                                                               calc_energy,
                                                               BLOCK_SIZE,
                                                               0));
    // TODO check the return value
    return num_block_per_SM * props.multiProcessorCount;
  }

  template<>
  int get_vina_score_batch<queue_cuda>(const int atoms, std::shared_ptr<queue_cuda> q_b) {
    // populate the bucket dimension
    int bucket_size{0};
    const int device_id = q_b->get_id();
    constexpr_for<0, reorder_buffer<static_molecule>::get_num_atom_clusters(), 1>([&](const auto atom_index) {
      const auto n_atoms = reorder_buffer<static_molecule>::atoms_clusters[atom_index];
      if (atoms == n_atoms)
        bucket_size = get_evaluate_fitness_batch<n_atoms>(device_id);
    });
    if (bucket_size == 0)
      throw std::runtime_error(
          "Compilation error: there is a bucket of atoms number which it is not handled.");

    mudock::info("CUDA Bucket size for ", atoms, " atoms ", bucket_size * BUCKET_MULTIPLIER, " ligands.");
    return bucket_size * BUCKET_MULTIPLIER;
  };
}

