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
#define MAX_LIGAND_ATOMS 256

    /**
     *  TODO: add the boolean problem in the relation https://gemini.google.com/share/eb3cf48bb2eb
     */


namespace mudock {
  
    typedef struct {
      fp_type x, y, z;
    } coords_t;

    __device__ static constexpr fp_type GAUSS1_COEFF_CUDA{- 0.035579f};
    __device__ static constexpr fp_type GAUSS2_COEFF_CUDA{- 0.005156f};
    __device__ static constexpr fp_type REPULSION_COEFF_CUDA{0.840245f};
    __device__ static constexpr fp_type HYDROPHOBIC_COEFF_CUDA{- 0.035069f};
    __device__ static constexpr fp_type H_BOND_COEFF_CUDA{- 0.587439f};
    __device__ static constexpr fp_type NROT_COEFF_CUDA{0.05846f};

    __device__ inline fp_type distance(fp_type x, fp_type y, fp_type z) {
        return sqrtf( x*x + y*y + z*z );
    }

    __device__ inline fp_type gauss1(const fp_type dst) {
      fp_type x = dst * 2.0f;
      return (dst != 0.0f) ? __expf(-(x*x)) : 0.0f;
    }

    __device__ inline fp_type gauss2(const fp_type dst) {
      fp_type x = (dst - 3) * 0.5f;
      return (dst != 0.0f) ? __expf(-(x*x)) : 0.0f;
    }

    __device__ inline fp_type repulsion(const fp_type dst) {
      return (dst < 0.0f) ? dst*dst : 0.0f;
    }

    __device__ inline fp_type hydrophobic(const fp_type dst, const int rec_lig_is_hydrophobic) {
      fp_type hydro_1 = (dst <= 0.5f) ? 1.0f : 0.0f;
      fp_type hydro_2 = (dst > 0.5f & dst < 1.5f) ? (1.5f - dst) : 0.0f;
      return (rec_lig_is_hydrophobic) ? hydro_1 + hydro_2 : 0.0f;
    }

    __device__ inline fp_type hbonding(const fp_type dst, const int rec_lig_is_hb) {
      fp_type h_bond_1 = (dst <= -0.7f) ? 1.0f : 0.0f;
      fp_type h_bond_2 = (dst < 0.0f & dst > -0.7f) ? (-dst * 1.42857f ) : 0.0f;  // 1.4285714285714286 = 1/0.7f
      return (rec_lig_is_hb) ? (h_bond_1 + h_bond_2) : 0.0f;
    }

    __device__ inline fp_type compute_pair_energy(
        fp_type dx, fp_type dy, fp_type dz,
        fp_type vdw1, fp_type vdw2,
        int is_hba1, int is_hbd1, int is_hba2, int is_hbd2,
        int is_hydro1, int is_hydro2
        ) {
      fp_type dst = distance(dx, dy, dz);

      dst -= (vdw1 + vdw2);
      const int is_h = (is_hba1 & is_hbd2) | (is_hba2 & is_hbd1);
      const int is_hydro = is_hydro1 & is_hydro2;

      fp_type res = GAUSS1_COEFF_CUDA * gauss1(dst) +
        GAUSS2_COEFF_CUDA * gauss2(dst) +
        REPULSION_COEFF_CUDA * repulsion(dst) +
        HYDROPHOBIC_COEFF_CUDA * hydrophobic(dst, is_hydro) +
        H_BOND_COEFF_CUDA * hbonding(dst, is_h);

      return (dst <= 8.0f) ? res : 0.0f;
    }

    __device__ inline fp_type score_inter(
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
        const coords_t* ligand_coords, 
        const int* __restrict__ l_is_hbond_acceptor,
        const int* __restrict__ l_is_hbond_donor,
        const int* __restrict__ l_is_hydrophobic,
        const fp_type* __restrict__ l_vdw_radius
        ) {
      fp_type total = 0;
      for (size_t pIdx = threadIdx.x; pIdx < num_atoms_protein; pIdx += blockDim.x) {

        fp_type px = protein_x[pIdx];
        fp_type py = protein_y[pIdx];
        fp_type pz = protein_z[pIdx];
        fp_type prv = p_vdw_radius[pIdx];
        int phba = p_is_hbond_acceptor[pIdx];
        int phbd = p_is_hbond_donor[pIdx];
        int phf = p_is_hydrophobic[pIdx];

        for (size_t lIdx = 0; lIdx < num_atoms_ligand; lIdx++) {
          total += compute_pair_energy(
              px - ligand_coords[lIdx].x,
              py - ligand_coords[lIdx].y,
              pz - ligand_coords[lIdx].z,
              prv, l_vdw_radius[lIdx],
              phba, phbd,
              l_is_hbond_acceptor[lIdx], l_is_hbond_donor[lIdx],
              phf, l_is_hydrophobic[lIdx]
              );
        }
      }
      return total;       
    }


    __device__ inline fp_type score_intra(
        const coords_t* ligand_coords, 
        const int* __restrict__ l_is_hbond_acceptor,
        const int* __restrict__ l_is_hbond_donor,
        const int* __restrict__ l_is_hydrophobic,
        const fp_type* __restrict__ l_vdw_radius,
        const int* __restrict__ interacting_pairs_first,
        const int* __restrict__ interacting_pairs_second,
        const size_t num_interacting_pairs
        ) {
      fp_type total = 0;
      for (size_t i = threadIdx.x; i < num_interacting_pairs; i += blockDim.x) {
        int a1 = interacting_pairs_first[i];
        int a2 = interacting_pairs_second[i];
        total += compute_pair_energy(
            ligand_coords[a1].x - ligand_coords[a2].x,
            ligand_coords[a1].y - ligand_coords[a2].y,
            ligand_coords[a1].z - ligand_coords[a2].z,
            l_vdw_radius[a1], l_vdw_radius[a2],
            l_is_hbond_acceptor[a1], l_is_hbond_donor[a1],
            l_is_hbond_acceptor[a2], l_is_hbond_donor[a2],
            l_is_hydrophobic[a1], l_is_hydrophobic[a2]
            );
      }
      return total;
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
        const coords_t* ligand_coords, 
        const int* __restrict__ l_is_hbond_acceptor,
        const int* __restrict__ l_is_hbond_donor,
        const int* __restrict__ l_is_hydrophobic,
        const fp_type* __restrict__ l_vdw_radius,
        const size_t active_torsions,
        const int* __restrict__ interacting_pairs_first,
        const int* __restrict__ interacting_pairs_second,
        const size_t num_interacting_pairs
    ){

        fp_type inter_score = score_inter(
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
            l_is_hbond_acceptor,
            l_is_hbond_donor,
            l_is_hydrophobic,
            l_vdw_radius
        );
        
        fp_type intra_score = score_intra(
            ligand_coords,
            l_is_hbond_acceptor,
            l_is_hbond_donor,
            l_is_hydrophobic,
            l_vdw_radius,
            interacting_pairs_first,
            interacting_pairs_second,
            num_interacting_pairs
        );

        fp_type score = (inter_score + intra_score) / ( 1 + NROT_COEFF_CUDA * active_torsions);
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
                              fp_type *__restrict__ scores                              
                              ) {

    const int ligand_id       = blockIdx.x;
    const int local_thread_id = threadIdx.x;

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

    /// TODO: check this cus i am not really sure about his behaviour in case of multiple scores per ligand
    const int* interacting_pairs_first = ligand_interacting_pairs_first + interacting_pairs_offset;
    const int* interacting_pairs_second = ligand_interacting_pairs_second + interacting_pairs_offset;

    fp_type* scores_l = scores + ligand_id * scores_per_ligand;

    __shared__ coords_t ligand_coords[MAX_LIGAND_ATOMS];

    for (int scores_index = 0; scores_index < scores_per_ligand; ++scores_index) {

      // Copy original coordinates
      for (int i = local_thread_id; i < num_atoms_ligand; i += blockDim.x) {
        ligand_coords[i].x = l_scratch_x[scores_index * atom_stride + i];
        ligand_coords[i].y = l_scratch_y[scores_index * atom_stride + i];
        ligand_coords[i].z = l_scratch_z[scores_index * atom_stride + i];
      }

      __syncwarp();

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
          ligand_coords,
          l_is_hbond_acceptor,
          l_is_hbond_donor,
          l_is_hydrophobic,
          l_vdw_radius,
          active_torsions,
          interacting_pairs_first, 
          interacting_pairs_second, 
          num_interacting_pairs
          );   

      __syncwarp();
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
                    (void*) &scores_b
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

