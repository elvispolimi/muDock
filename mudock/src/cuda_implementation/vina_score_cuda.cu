#include "mudock/cuda_implementation/queue_cuda.cuh"
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

#define VINA_BLOCK_SIZE 256
#define BUCKET_MULTIPLIER 22

#ifndef EPSILON
#define EPSILON 1e-6f
#endif

#define IS_DIFF_FROM_ZERO(dst) (fabsf(dst) > EPSILON)

namespace mudock {
  
    __device__ static constexpr fp_type GAUSS1_COEFF_CUDA{- 0.035579f};
    __device__ static constexpr fp_type GAUSS2_COEFF_CUDA{- 0.005156f};
    __device__ static constexpr fp_type REPULSION_COEFF_CUDA{0.840245f};
    __device__ static constexpr fp_type HYDROPHOBIC_COEFF_CUDA{- 0.035069f};
    __device__ static constexpr fp_type H_BOND_COEFF_CUDA{- 0.587439f};
    __device__ static constexpr fp_type NROT_COEFF_CUDA{0.05846f};

#define PARAM_HA    x
#define PARAM_HD    y
#define PARAM_HYDRO z
#define PARAM_VDW   w

#define get_ligand_par(l_par, par) ((l_par).PARAM_##par)

__device__ inline float warp_reduce_sum(float val)
{
    constexpr unsigned int FULL_MASK{0xffffffff};
#pragma unroll
    for (size_t offset{16}; offset > 0; offset /= 2)
    {
        val += __shfl_down_sync(FULL_MASK, val, offset);
    }
    // Only the first thread in the warp will return the correct result.
    return val;
}

template <size_t NUM_THREADS, size_t NUM_WARPS = NUM_THREADS / 32>
__device__ inline fp_type block_reduce_sum_v2(fp_type val, fp_type shared_data[NUM_WARPS])
{
    val = warp_reduce_sum(val);

    if (threadIdx.x % 32 == 0)
    {
        shared_data[threadIdx.x / 32] = val;
    }

    __syncthreads();

    // First warp reduce 
    fp_type block_sum{0.0f};
    if (threadIdx.x < 32)
    {
        fp_type warp_val = (threadIdx.x < NUM_WARPS) ? shared_data[threadIdx.x] : 0.0f;
        block_sum = warp_reduce_sum(warp_val);
    }

    return block_sum;
}

    __device__ inline fp_type gauss1(const fp_type dst) {
      fp_type x = dst * 2.0f;
      return (fp_type)(IS_DIFF_FROM_ZERO(dst)) * __expf(-(x*x));
    }

    __device__ inline fp_type gauss2(const fp_type dst) {
      fp_type x = (dst - 3.0f) * 0.5f;
      return (fp_type)(IS_DIFF_FROM_ZERO(dst)) * __expf(-(x*x));
    }

    __device__ inline fp_type repulsion(const fp_type dst) {
      return (fp_type)(dst < 0.0f) * (dst * dst);
    }

    __device__ inline fp_type hydrophobic(const fp_type dst, const fp_type is_hydro) {
      fp_type hydro_1 = (fp_type)(dst <= 0.5f); 
      fp_type val = 1.5f - dst;
      fp_type hydro_2 = fmaxf(0.0f, fminf(1.0f, val)) * (fp_type)(dst > 0.5f);
      return is_hydro * (hydro_1 + hydro_2); 
    }

    __device__ inline fp_type hbonding(const fp_type dst, const fp_type is_hb) {
      fp_type h_bond_1 = (fp_type)(dst <= -0.7f);
      fp_type val = -dst * 1.42857f;
      fp_type h_bond_2 = fmaxf(0.0f, fminf(1.0f, val)) * (fp_type)(dst > -0.7f && dst < 0.0f);
      return is_hb * (h_bond_1 + h_bond_2);
    }

    __device__ inline fp_type compute_pair_energy(
        fp_type x1, fp_type y1, fp_type z1,
        fp_type x2, fp_type y2, fp_type z2,
        fp_type vdw1, fp_type vdw2,
        fp_type is_hba1, fp_type is_hbd1, fp_type is_hba2, fp_type is_hbd2,
        fp_type is_hydro1, fp_type is_hydro2
        ) {

      const fp_type dx = x1 - x2;
      const fp_type dy = y1 - y2;
      const fp_type dz = z1 - z2;
      fp_type d2 = (dx * dx) + (dy * dy) + (dz * dz);

      if (d2 > 64.0f) return 0.0f;

      fp_type dst = sqrtf(d2);

      dst -= (vdw1 + vdw2);

      const fp_type is_h = fminf(1.0f, (is_hba1 * is_hbd2) + (is_hba2 * is_hbd1));
      const fp_type is_hydro = is_hydro1 * is_hydro2;

      fp_type res = GAUSS1_COEFF_CUDA * gauss1(dst) +
        GAUSS2_COEFF_CUDA * gauss2(dst) +
        REPULSION_COEFF_CUDA * repulsion(dst) +
        HYDROPHOBIC_COEFF_CUDA * hydrophobic(dst, is_hydro) +
        H_BOND_COEFF_CUDA * hbonding(dst, is_h);

      return res;
    }


 template<int MAX_ATOMS>
    __device__ inline fp_type score_inter(
        /// Protein data
        const int num_atoms_protein,
        const fp_type* __restrict__ protein_x,
        const fp_type* __restrict__ protein_y,
        const fp_type* __restrict__ protein_z,
        const int* __restrict__ p_is_hbond_acceptor,
        const int* __restrict__ p_is_hbond_donor,
        const int* __restrict__ p_is_hydrophobic,
        const fp_type* __restrict__ p_vdw_radius,

        ///Ligand data
        const int num_atoms_ligand,
        const float4* ligand_coords, 
        const float4* ligand_par
        ) {
      fp_type total = 0;
      for (int pIdx = threadIdx.x; pIdx < num_atoms_protein; pIdx += blockDim.x) {

        const fp_type px = protein_x[pIdx];
        const fp_type py = protein_y[pIdx];
        const fp_type pz = protein_z[pIdx];
        const fp_type pvdw = p_vdw_radius[pIdx];
        // Explicit conversion to fp_type
        const fp_type phba = (fp_type) p_is_hbond_acceptor[pIdx];
        const fp_type phbd = (fp_type) p_is_hbond_donor[pIdx];
        const fp_type phf  = (fp_type) p_is_hydrophobic[pIdx];

#pragma unroll 8 // unsing unroll on MAX_ATOMS cause register spilling. Trying I found 8 as a good balance.
        for (int lIdx = 0; lIdx < num_atoms_ligand; lIdx++) {

          const float4 coords = ligand_coords[lIdx];
          const float4 par    = ligand_par[lIdx];

          total += compute_pair_energy(
              px, coords.x, 
              py, coords.y,
              pz, coords.z,
              pvdw, get_ligand_par(par, VDW),
              phba, phbd,
              get_ligand_par(par, HA), get_ligand_par(par, HD),
              phf, get_ligand_par(par, HYDRO)
              );   
        }
      }
      return total;       
    }


    __device__ inline fp_type score_intra(
        const float4* ligand_coords, 
        const float4* ligand_par,
        const int* __restrict__ interacting_pairs_first,
        const int* __restrict__ interacting_pairs_second,
        const int num_interacting_pairs
        ) {
      fp_type total = 0;
      for (int i = threadIdx.x; i < num_interacting_pairs; i += blockDim.x) {
        const int a1 = interacting_pairs_first[i];
        const int a2 = interacting_pairs_second[i];
 
        const float4 a1c = ligand_coords[a1];
        const float4 a2c = ligand_coords[a2];
        
        const float4 a1p = ligand_par[a1];
        const float4 a2p = ligand_par[a2];

        total += compute_pair_energy(
            a1c.x, a2c.x,
            a1c.y, a2c.y,
            a1c.z, a2c.z,
            get_ligand_par(a1p, VDW),   get_ligand_par(a2p, VDW),
            get_ligand_par(a1p, HA),    get_ligand_par(a1p, HD),
            get_ligand_par(a2p, HA),    get_ligand_par(a2p, HD),
            get_ligand_par(a1p, HYDRO), get_ligand_par(a2p, HYDRO)
            );     
      }
      return total;
    }

 template<int MAX_ATOMS>
    __device__ inline fp_type scoring_cuda(  
        /// Protein data
        const int num_atoms_protein,
        const fp_type* __restrict__ protein_x,
        const fp_type* __restrict__ protein_y,
        const fp_type* __restrict__ protein_z,
        const int* __restrict__ p_is_hbond_acceptor,
        const int* __restrict__ p_is_hbond_donor,
        const int* __restrict__ p_is_hydrophobic,
        const fp_type* __restrict__ p_vdw_radius,

        ///Ligand data
        const int num_atoms_ligand,
        const float4* ligand_coords, 
        const float4* l_par,
        const int active_torsions,
        const int* __restrict__ interacting_pairs_first,
        const int* __restrict__ interacting_pairs_second,
        const int num_interacting_pairs
    ){

        fp_type inter_score = score_inter<MAX_ATOMS>(
            num_atoms_protein,
            protein_x,
            protein_y,
            protein_z,
            p_is_hbond_acceptor,
            p_is_hbond_donor,
            p_is_hydrophobic,
            p_vdw_radius,
            num_atoms_ligand,
            ligand_coords, 
            l_par
        );
        
        fp_type intra_score = score_intra(
            ligand_coords, 
            l_par,
            interacting_pairs_first,
            interacting_pairs_second,
            num_interacting_pairs
        );

        fp_type score = (inter_score + intra_score) / ( 1 + NROT_COEFF_CUDA * active_torsions);

        // if(threadIdx.x == 0) printf("Score inter %f, Score intra %f, Score %f\n", inter_score, intra_score, score); 
        return score;
    }

template<int MAX_ATOMS>
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
                              const int num_atoms_protein,
                              const fp_type* __restrict__ protein_x,
                              const fp_type* __restrict__ protein_y,
                              const fp_type* __restrict__ protein_z,
                              const int* __restrict__ p_is_hbond_acceptor,
                              const int* __restrict__ p_is_hbond_donor,
                              const int* __restrict__ p_is_hydrophobic,
                              const fp_type* __restrict__ p_vdw_radius,
                              fp_type *__restrict__ scores                              
                              ) {

    const int ligand_id       = blockIdx.x;
    const int local_thread_id = threadIdx.x;

    //if(local_thread_id == 0) printf("MAX_ATOMS: %d\n", MAX_ATOMS);
  
    const int num_atoms_ligand    = num_atoms_b[ligand_id];
    const int active_torsions = ligand_active_torsions[ligand_id];
    const int num_interacting_pairs = ligand_num_interacting_pairs[ligand_id];
    const int interacting_pairs_offset = ligand_interacting_pairs_offset[ligand_id];
    const int stride       = ligand_id * atom_stride;

    const fp_type* l_scratch_x = x_scratch + stride * scores_per_ligand;
    const fp_type* l_scratch_y = y_scratch + stride * scores_per_ligand;
    const fp_type* l_scratch_z = z_scratch + stride * scores_per_ligand;
    const int* l_is_hbond_acceptor = ligand_is_hbond_acceptor + stride;
    const int* l_is_hbond_donor = ligand_is_hbond_donor + stride;
    const int* l_is_hydrophobic = ligand_is_hydrophobic + stride;
    const fp_type* l_vdw_radius = ligand_vdw_radius + stride;

    const int* interacting_pairs_first = ligand_interacting_pairs_first + interacting_pairs_offset;
    const int* interacting_pairs_second = ligand_interacting_pairs_second + interacting_pairs_offset;

    fp_type* scores_l = scores + ligand_id * scores_per_ligand;

    __shared__ float4 ligand_coords[MAX_ATOMS];
    __shared__ float4 ligand_par[MAX_ATOMS];
    __shared__ fp_type warp_shared_buffer[VINA_BLOCK_SIZE/32];

    for (int i = local_thread_id; i < num_atoms_ligand; i += blockDim.x) {
      ligand_par[i] = make_float4(
          (float) l_is_hbond_acceptor[i], 
          (float) l_is_hbond_donor[i], 
          (float) l_is_hydrophobic[i], 
          (float) l_vdw_radius[i]
          );  
    }

    __syncthreads();

    for (int scores_index = 0; scores_index < scores_per_ligand; ++scores_index) {

      const fp_type* ligand_x = l_scratch_x + scores_index * atom_stride;
      const fp_type* ligand_y = l_scratch_y + scores_index * atom_stride;
      const fp_type* ligand_z = l_scratch_z + scores_index * atom_stride;

      for (int i = local_thread_id; i < num_atoms_ligand; i += blockDim.x) {
        ligand_coords[i] = make_float4(ligand_x[i], ligand_y[i], ligand_z[i], .0f);
      }
      __syncthreads();

      // Calculate energy 
      fp_type result = scoring_cuda<MAX_ATOMS>(num_atoms_protein, 
          protein_x, 
          protein_y,
          protein_z, 
          p_is_hbond_acceptor, 
          p_is_hbond_donor, 
          p_is_hydrophobic, 
          p_vdw_radius, 
          num_atoms_ligand, 
          ligand_coords, 
          ligand_par,
          active_torsions,
          interacting_pairs_first, 
          interacting_pairs_second, 
          num_interacting_pairs
          );   

#ifdef MUDOCK_TEST
#else
#endif

      fp_type final_result = block_reduce_sum_v2<VINA_BLOCK_SIZE>(result, warp_shared_buffer);

      if (local_thread_id == 0) {
        scores_l[scores_index]         = final_result;
      }

      __syncthreads(); 
    }
  }

  template<>
  void vina_score_kernel<queue_cuda>::operator()() {
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
                    (void*) &scores_b
    };
    constexpr_switch_bucket<0, reorder_buffer<static_molecule>::get_num_atom_clusters(), 1>(
        [&](const auto atom_index) {
          const auto max_atoms = reorder_buffer<static_molecule>::atoms_clusters[atom_index];
          q->launch_kernel((void*) calc_energy<max_atoms>, args, batch_ligands, VINA_BLOCK_SIZE);
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
                                                               calc_energy<MAX_ATOMS>,
                                                               VINA_BLOCK_SIZE,
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

