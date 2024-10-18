#include <alloca.h>
#include <cstddef>
#include <cstring>
#include <cuda_runtime.h>
#include <mudock/chem/autodock_parameters.hpp>
#include <mudock/chem/grid_const.hpp>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/cuda_implementation/cuda_check_error_macro.cuh>
#include <polygeist/cuda_random.cuh>
#include <mudock/grid/mdindex.hpp>
#include <mudock/grid/point3D.hpp>
#include <mudock/type_alias.hpp>
#include <mudock/utils.hpp>

#define FLATTENED_3D(x, y, z, index_n_x, index_n_xy) (index_n_xy * z + y * index_n_x + x)
#define MAX(a, b)                                    ((a) > (b) ? (a) : (b))
#define MIN(a, b)                                    ((a) < (b) ? (a) : (b))
// Keep it to 32 to enable warp optimizations
#define BLOCK_SIZE 32

// TODO implement better way to define population
// NOTE: the population size can be at most 100
#define POPULATION_NUMBER 100
#define SHARED_MEM        (BLOCK_SIZE * 2 + 100)
namespace mudock {
  // calc_energy.cu
  __device__ fp_type calc_ddd_Mehler_Solmajer_cuda(fp_type distance) {
    const fp_type lambda{0.003627};
    const fp_type epsilon0{78.4};
    const fp_type A{-8.5525};
    const fp_type B = epsilon0 - A;
    const fp_type rk{7.7839};
    const fp_type lambda_B = -lambda * B;

    fp_type epsilon = A + B / (fp_type{1} + rk * expf(lambda_B * distance));

    if (epsilon < std::numeric_limits<fp_type>::epsilon()) {
      epsilon = fp_type{1.0};
    }
    return epsilon;
  }

  // TODO check math functions if they are from CUDA
  __device__ fp_type calc_intra_energy(const fp_type* ligand_x,
                                       const fp_type* ligand_y,
                                       const fp_type* ligand_z,
                                       const fp_type* ligand_vol,
                                       const fp_type* ligand_solpar,
                                       const fp_type* ligand_charge,
                                       const int* ligand_num_hbond,
                                       const fp_type* ligand_Rij_hb,
                                       const fp_type* ligand_Rii,
                                       const fp_type* ligand_epsij_hb,
                                       const fp_type* ligand_epsii,
                                       const int ligand_num_nonbonds,
                                       const int* __restrict__ ligand_nonbond_a1,
                                       const int* __restrict__ ligand_nonbond_a2) {
    fp_type elect_total_eintcal{0}, emap_total_eintcal{0}, dmap_total_eintcal{0};
    // TODO OPT: precompute value from branch rework_cpp
    for (int nonbond_list = threadIdx.x; nonbond_list < ligand_num_nonbonds; nonbond_list += blockDim.x) {
      const int& a1 = ligand_nonbond_a1[nonbond_list];
      const int& a2 = ligand_nonbond_a2[nonbond_list];

      const fp_type distance_two = powf(fabs(ligand_x[a1] - ligand_x[a2]), fp_type{2}) +
                                   powf(fabs(ligand_y[a1] - ligand_y[a2]), fp_type{2}) +
                                   powf(fabs(ligand_z[a1] - ligand_z[a2]), fp_type{2});
      const fp_type distance_two_clamp = MAX(distance_two, RMIN_ELEC * RMIN_ELEC);
      const fp_type distance           = std::sqrt(distance_two_clamp);

      //  Calculate  Electrostatic  Energy
      const fp_type r_dielectric = fp_type{1} / (distance * calc_ddd_Mehler_Solmajer_cuda(distance));
      const fp_type e_elec =
          ligand_charge[a1] * ligand_charge[a2] * ELECSCALE * autodock_parameters::coeff_estat * r_dielectric;
      elect_total_eintcal += e_elec;

      // Calcuate desolv
      const fp_type nb_desolv = (ligand_vol[a2] * (ligand_solpar[a1] + qsolpar * fabsf(ligand_charge[a1])) +
                                 ligand_vol[a1] * (ligand_solpar[a2] + qsolpar * fabsf(ligand_charge[a2])));

      const fp_type e_desolv = autodock_parameters::coeff_desolv *
                               expf(fp_type{-0.5} / (sigma * sigma) * distance_two_clamp) * nb_desolv;
      dmap_total_eintcal += e_desolv;

      fp_type e_vdW_Hb{0};
      if (distance_two_clamp < nbc2) {
        const int& hbond_i        = ligand_num_hbond[a1];
        const int& hbond_j        = ligand_num_hbond[a2];
        const fp_type& Rij_hb_i   = ligand_Rij_hb[a1];
        const fp_type& Rij_hb_j   = ligand_Rij_hb[a2];
        const fp_type& Rii_i      = ligand_Rii[a1];
        const fp_type& Rii_j      = ligand_Rii[a2];
        const fp_type& epsij_hb_i = ligand_epsij_hb[a1];
        const fp_type& epsij_hb_j = ligand_epsij_hb[a2];
        const fp_type& epsii_i    = ligand_epsii[a1];
        const fp_type& epsii_j    = ligand_epsii[a2];

        int xA = 12;
        int xB = 6;

        fp_type Rij{0}, epsij{0};
        if ((hbond_i == 1 || hbond_i == 2) && hbond_j > 2) {
          Rij   = Rij_hb_j;
          epsij = epsij_hb_j;
          xB    = 10;
        } else if ((hbond_i > 2) && (hbond_j == 1 || hbond_j == 2)) {
          Rij   = Rij_hb_i;
          epsij = epsij_hb_i;
          xB    = 10;
        } else {
          Rij   = (Rii_i + Rii_j) / fp_type{2};
          epsij = sqrt(epsii_i * epsii_j);
        }
        if (xA != xB) {
          const fp_type tmp = epsij / (xA - xB);
          const fp_type cA  = tmp * powf(Rij, xA) * xB;
          const fp_type cB  = tmp * powf(Rij, xB) * xA;

          const fp_type rA = powf(distance, static_cast<fp_type>(xA));
          const fp_type rB = powf(distance, static_cast<fp_type>(xB));

          // TODO EINTCLAMP problem polygeist with constexpr?
          e_vdW_Hb = MIN(100000, (cA / rA - cB / rB));
        }
      }
      emap_total_eintcal += e_vdW_Hb;
    }
    return elect_total_eintcal + emap_total_eintcal + dmap_total_eintcal;
  }

  // geometric_transformations.cu
  __device__ void translate_molecule_cuda(fp_type* __restrict__ x,
                                          fp_type* __restrict__ y,
                                          fp_type* __restrict__ z,
                                          const fp_type* offset_x,
                                          const fp_type* offset_y,
                                          const fp_type* offset_z,
                                          const int num_atoms) {
    for (int i = threadIdx.x; i < num_atoms; i += blockDim.x) {
      x[i] += *offset_x;
      y[i] += *offset_y;
      z[i] += *offset_z;
    }
  }

  __device__ void rotate_molecule_cuda(fp_type* __restrict__ x,
                                       fp_type* __restrict__ y,
                                       fp_type* __restrict__ z,
                                       const fp_type* angle_x,
                                       const fp_type* angle_y,
                                       const fp_type* angle_z,
                                       const int num_atoms) {
    // compute the angles sine and cosine
    const auto rad_x = deg_to_rad(*angle_x), rad_y = deg_to_rad(*angle_y), rad_z = deg_to_rad(*angle_z);
    const auto cx = std::cos(rad_x), sx = std::sin(rad_x);
    const auto cy = std::cos(rad_y), sy = std::sin(rad_y);
    const auto cz = std::cos(rad_z), sz = std::sin(rad_z);

    // compute the rotation matrix defined as Rz*Ry*Rx
    const auto m00 = cy * cz;
    const auto m01 = sx * sy * cz - cx * sz;
    const auto m02 = cx * sy * cz + sx * sz;
    const auto m10 = cy * sz;
    const auto m11 = sx * sy * sz + cx * cz;
    const auto m12 = cx * sy * sz - sx * cz;
    const auto m20 = -sy;
    const auto m21 = sx * cy;
    const auto m22 = cx * cy;

    // apply the rotation matrix
    for (int i = threadIdx.x; i < num_atoms; i += blockDim.x) {
      const auto prev_x = x[i], prev_y = y[i], prev_z = z[i];
      x[i] = prev_x * m00 + prev_y * m01 + prev_z * m02;
      y[i] = prev_x * m10 + prev_y * m11 + prev_z * m12;
      z[i] = prev_x * m20 + prev_y * m21 + prev_z * m22;
    }
  }

  __device__ void rotate_fragment_cuda(fp_type* __restrict__ x,
                                       fp_type* __restrict__ y,
                                       fp_type* __restrict__ z,
                                       const int* bitmask,
                                       const int start_index,
                                       const int stop_index,
                                       const fp_type* angle,
                                       const int num_atoms) {
    // compute the axis vector (and some properties)
    const auto origx = x[start_index], origy = y[start_index], origz = z[start_index];
    const auto destx = x[stop_index], desty = y[stop_index], destz = z[stop_index];
    const auto u = destx - origx;
    const auto v = desty - origy;
    const auto w = destz - origz;

    const auto u2 = u * u, v2 = v * v, w2 = w * w;
    const auto l2 = u * u + v * v + w * w;
    // Check if origin and dest coincide
    // No need to continue the intramolecular energy will be very high
    if (isinf(l2) || l2 == fp_type{0} || isnan(l2))
      // TODO print error?
      return;
    const auto l = std::sqrt(l2);

    // compute the angle sine and cosine
    const auto rad = deg_to_rad(*angle);
    const auto s = std::sin(rad), c = std::cos(rad);
    const auto one_minus_c = fp_type{1} - c;
    const auto ls          = l * s;

    // compute the rotation matrix (rodrigues' rotation formula)
    const auto m00 = (u2 + (v2 + w2) * c) / l2;
    const auto m01 = (u * v * one_minus_c - w * l * s) / l2;
    const auto m02 = (u * w * one_minus_c + v * l * s) / l2;
    const auto m03 =
        ((origx * (v2 + w2) - u * (origy * v + origz * w)) * one_minus_c + (origy * w - origz * v) * ls) / l2;

    const auto m10 = (u * v * one_minus_c + w * ls) / l2;
    const auto m11 = (v2 + (u2 + w2) * c) / l2;
    const auto m12 = (v * w * one_minus_c - u * ls) / l2;
    const auto m13 =
        ((origy * (u2 + w2) - v * (origx * u + origz * w)) * one_minus_c + (origz * u - origx * w) * ls) / l2;

    const auto m20 = (u * w * one_minus_c - v * ls) / l2;
    const auto m21 = (v * w * one_minus_c + u * ls) / l2;
    const auto m22 = (w2 + (u2 + v2) * c) / l2;
    const auto m23 =
        ((origz * (u2 + v2) - w * (origx * u + origy * v)) * one_minus_c + (origx * v - origy * u) * ls) / l2;

    // apply the rotation matrix
    for (int i = threadIdx.x; i < num_atoms; i += blockDim.x) {
      if (bitmask[i] == 1) {
        const auto prev_x = x[i], prev_y = y[i], prev_z = z[i];
        x[i] = prev_x * m00 + prev_y * m01 + prev_z * m02 + m03;
        y[i] = prev_x * m10 + prev_y * m11 + prev_z * m12 + m13;
        z[i] = prev_x * m20 + prev_y * m21 + prev_z * m22 + m23;
      }
    }
  }

  // mutate.cu
  __device__ void apply_cuda(fp_type* __restrict__ x,
                             fp_type* __restrict__ y,
                             fp_type* __restrict__ z,
                             const chromosome& chromosome,
                             const int* __restrict__ fragments,
                             const int* __restrict__ fragments_start_index,
                             const int* __restrict__ fragments_stop_index,
                             const int num_rotamers,
                             const int stride_atoms,
                             const int num_atoms) {
    // apply rigid transformations
    translate_molecule_cuda(x, y, z, &chromosome[0], &chromosome[1], &chromosome[2], num_atoms);
    rotate_molecule_cuda(x, y, z, &chromosome[3], &chromosome[4], &chromosome[5], num_atoms);

    // change the molecule shape
    for (int i = 0; i < num_rotamers; ++i) {
      const int* bitmask = fragments + i * stride_atoms;
      rotate_fragment_cuda(x,
                           y,
                           z,
                           bitmask,
                           fragments_start_index[i],
                           fragments_stop_index[i],
                           &chromosome[6 + i],
                           num_atoms);
    }
  }

  // evaluate_fitness.cu
  static constexpr fp_type coordinate_step{0.2};
  static constexpr fp_type angle_step{4};

  __device__ fp_type trilinear_interpolation_cuda(const fp_type coord[],
                                                  const fp_type* tex,
                                                  const int index_n_x,
                                                  const int index_n_xy) {
    // Interpolation CUDA
    const int u0      = coord[0];
    const fp_type p0u = coord[0] - static_cast<fp_type>(u0);
    const fp_type p1u = fp_type{1} - p0u;

    const int v0      = coord[1];
    const fp_type p0v = coord[1] - static_cast<fp_type>(v0);
    const fp_type p1v = fp_type{1} - p0v;

    const int w0      = coord[2];
    const fp_type p0w = coord[2] - static_cast<fp_type>(w0);
    const fp_type p1w = fp_type{1} - p0w;

    const fp_type pu[2] = {p1u, p0u};
    const fp_type pv[2] = {p1v, p0v};
    const fp_type pw[2] = {p1w, p0w};
    fp_type value{0};
#pragma unroll
    for (int i = 0; i <= 1; i++)
#pragma unroll
      for (int t = 0; t <= 1; t++)
#pragma unroll
        for (int n = 0; n <= 1; n++) {
          const fp_type tmp = tex[FLATTENED_3D(u0 + n, v0 + t, w0 + i, index_n_x, index_n_xy)];
          value += pu[n] * pv[t] * pw[i] * tmp;
        }
    return value;
  }

    // Generate the next random number
  __device__ fp_type gen_random(XORWOWState& st) {
    /* Algorithm "xorwow" from p. 5 of Marsaglia, "Xorshift RNGs" */
    unsigned int t = st.state[4];

    const unsigned int s = st.state[0]; /* Perform a contrived 32-bit rotate. */
    st.state[4]             = st.state[3];
    st.state[3]             = st.state[2];
    st.state[2]             = st.state[1];
    st.state[1]             = s;

    t ^= t >> 2;
    t ^= t << 1;
    t ^= s ^ (s << 4);
    st.state[0] = t;
    st.index += 362437;
    return static_cast<fp_type>(t + st.index) / static_cast<fp_type>(std::numeric_limits<unsigned int>::max());
  }

  template<typename T>
  __device__ const T random_gen_cuda(XORWOWState& state, const T min, const T max) {
    fp_type value{0};
    if constexpr (is_debug()) {
      // TODO value here for debug
      value = fp_type{0.4};
    } else {
      value = gen_random(state);
    }
    return static_cast<T>((value * static_cast<fp_type>(max - min)) + min);
  }

  __device__ int get_selection_distribution(XORWOWState& state, const int* population_number) {
    return random_gen_cuda<int>(state, 0, *population_number - 1);
  };

  __device__ fp_type get_init_change_distribution(XORWOWState& state) {
    return random_gen_cuda<fp_type>(state, -45, 45);
  }
  __device__ fp_type get_mutation_change_distribution(XORWOWState& state) {
    return random_gen_cuda<fp_type>(state, -10, 10);
  };
  __device__ fp_type get_mutation_coin_distribution(XORWOWState& state) {
    return random_gen_cuda<fp_type>(state, 0, 1);
  };
  __device__ int get_crossover_distribution(XORWOWState& state, const int* num_rotamers) {
    return random_gen_cuda<int>(state, 0, 6 + *num_rotamers);
  };

  __device__ int tournament_selection_cuda(XORWOWState& state,
                                           const int tournament_length,
                                           const int chromosome_number,
                                           const fp_type* __restrict__ scores) {
    const int num_iterations = tournament_length;
    int best_individual      = get_selection_distribution(state, &chromosome_number);
    for (int i = 0; i < num_iterations; ++i) {
      auto contended = get_selection_distribution(state, &chromosome_number);
      if (scores[contended] < scores[best_individual]) {
        best_individual = contended;
      }
    }
    return best_individual;
  }

  // TODO check the syncthreads
  //TODO OPT: template parameter based on number of atoms, rotamers, chromosomes and population
  // Interesting the usage of the bucketizer
  __global__ void evaluate_fitness(const int num_generations,
                                   const int tournament_length,
                                   const fp_type mutation_prob,
                                   const int chromosome_number,
                                   const int chromosome_stride,
                                   const int atom_stride,
                                   const int rotamers_stride,
                                   const int nonbond_stride,
                                   const fp_type* __restrict__ original_ligand_x,
                                   const fp_type* __restrict__ original_ligand_y,
                                   const fp_type* __restrict__ original_ligand_z,
                                   fp_type* __restrict__ scratch_ligand_x,
                                   fp_type* __restrict__ scratch_ligand_y,
                                   fp_type* __restrict__ scratch_ligand_z,
                                   const fp_type* __restrict__ ligand_vol,
                                   const fp_type* __restrict__ ligand_solpar,
                                   const fp_type* __restrict__ ligand_charge,
                                   const int* __restrict__ ligand_num_hbond,
                                   const fp_type* __restrict__ ligand_Rij_hb,
                                   const fp_type* __restrict__ ligand_Rii,
                                   const fp_type* __restrict__ ligand_epsij_hb,
                                   const fp_type* __restrict__ ligand_epsii,
                                   const int* __restrict__ ligand_num_nonbonds,
                                   const int* __restrict__ ligand_nonbond_a1,
                                   const int* __restrict__ ligand_nonbond_a2,
                                   const int* __restrict__ ligand_num_atoms,
                                   const int* __restrict__ ligand_num_rotamers,
                                   const int* __restrict__ ligand_fragments,
                                   const int* __restrict__ frag_start_atom_index,
                                   const int* __restrict__ frag_stop_atom_index,
                                   chromosome* __restrict__ chromosomes,
                                   const fp_type minimum_x,
                                   const fp_type minimum_y,
                                   const fp_type minimum_z,
                                   const fp_type maximum_x,
                                   const fp_type maximum_y,
                                   const fp_type maximum_z,
                                   const fp_type center_x,
                                   const fp_type center_y,
                                   const fp_type center_z,
                                   const int index_n_x,
                                   const int index_n_xy,
                                   const fp_type* const __restrict__* const __restrict__ atom_textures,
                                   const int* __restrict__ atom_tex_indexes,
                                   const fp_type* __restrict__ electro_texture,
                                   const fp_type* __restrict__ desolv_texture,
                                   XORWOWState* __restrict__ state,
                                   fp_type* __restrict__ ligand_scores,
                                   chromosome* __restrict__ best_chromosomes) {
    const int ligand_id        = blockIdx.x;
    const int local_thread_id  = threadIdx.x;
    const int thread_per_block = blockDim.x;
    const int global_thread_id = local_thread_id + thread_per_block * ligand_id;

    const int num_atoms    = ligand_num_atoms[ligand_id];
    const int num_nonbonds = ligand_num_nonbonds[ligand_id];
    const int num_rotamers = ligand_num_rotamers[ligand_id];

    const fp_type* l_original_ligand_x = original_ligand_x + ligand_id * atom_stride;
    const fp_type* l_original_ligand_y = original_ligand_y + ligand_id * atom_stride;
    const fp_type* l_original_ligand_z = original_ligand_z + ligand_id * atom_stride;
    fp_type* l_scratch_ligand_x        = scratch_ligand_x + ligand_id * atom_stride;
    fp_type* l_scratch_ligand_y        = scratch_ligand_y + ligand_id * atom_stride;
    fp_type* l_scratch_ligand_z        = scratch_ligand_z + ligand_id * atom_stride;
    const fp_type* l_ligand_vol        = ligand_vol + ligand_id * atom_stride;
    const fp_type* l_ligand_solpar     = ligand_solpar + ligand_id * atom_stride;
    const fp_type* l_ligand_charge     = ligand_charge + ligand_id * atom_stride;
    const int* l_ligand_num_hbond      = ligand_num_hbond + ligand_id * atom_stride;
    const fp_type* l_ligand_Rij_hb     = ligand_Rij_hb + ligand_id * atom_stride;
    const fp_type* l_ligand_Rii        = ligand_Rii + ligand_id * atom_stride;
    const fp_type* l_ligand_epsij_hb   = ligand_epsij_hb + ligand_id * atom_stride;
    const fp_type* l_ligand_epsii      = ligand_epsii + ligand_id * atom_stride;
    chromosome* l_chromosomes          = chromosomes + ligand_id * chromosome_stride;
    // Point to the next population buffer
    chromosome* l_next_chromosomes = chromosomes + ligand_id * chromosome_stride + chromosome_number;
    // Output chromosome with the best one
    const auto* l_fragments             = ligand_fragments + ligand_id * atom_stride * rotamers_stride;
    const auto* l_frag_start_atom_index = frag_start_atom_index + ligand_id * rotamers_stride;
    const auto* l_frag_stop_atom_index  = frag_stop_atom_index + ligand_id * rotamers_stride;
    const auto* l_atom_tex_indexes      = atom_tex_indexes + ligand_id * atom_stride;
    const int* l_ligand_nonbond_a1      = ligand_nonbond_a1 + ligand_id * nonbond_stride;
    const int* l_ligand_nonbond_a2      = ligand_nonbond_a2 + ligand_id * nonbond_stride;
    XORWOWState& l_state                = (state[global_thread_id]);

    // Shared memory
    __shared__ fp_type shared_data[SHARED_MEM];
    fp_type* s_reduction_first   = shared_data;
    fp_type* s_reduction_second  = shared_data + blockDim.x;
    fp_type* s_chromosome_scores = shared_data + 2 * blockDim.x;
    // Initialize shared mem
    for (int chromosome_index = local_thread_id; chromosome_index < MAX(chromosome_number, thread_per_block);
         chromosome_index += thread_per_block)
      s_chromosome_scores[chromosome_index] = 0x7ff0000000000000; // Set initial score value, infinity value

    // Generate initial population
    for (int chromosome_index = local_thread_id; chromosome_index < chromosome_number;
         chromosome_index += thread_per_block) {
      chromosome& chromo = *(l_chromosomes + chromosome_index);
#pragma unroll
      for (int i{0}; i < 3; ++i) { // initialize the rigid translation
        chromo[i] = get_init_change_distribution(l_state) * coordinate_step;
      }
#pragma unroll
      for (int i{3}; i < 6 + num_rotamers; ++i) { // initialize the rotations
        chromo[i] = get_init_change_distribution(l_state) * angle_step;
      }
    }
    __syncthreads();
    // TODO maybe template parameter?
    for (int generation = 0; generation < num_generations; ++generation) {
      for (int chromosome_index = 0; chromosome_index < chromosome_number; ++chromosome_index) {
        // Copy original coordinates
        // TODO OPT: shared memory for coordinate ?
        for (int atom_index = local_thread_id; atom_index < num_atoms; atom_index += thread_per_block) {
          l_scratch_ligand_x[atom_index] = l_original_ligand_x[atom_index];
          l_scratch_ligand_y[atom_index] = l_original_ligand_y[atom_index];
          l_scratch_ligand_z[atom_index] = l_original_ligand_z[atom_index];
        }
        // Modify coordinates
        apply_cuda(l_scratch_ligand_x,
                   l_scratch_ligand_y,
                   l_scratch_ligand_z,
                   *(l_chromosomes + chromosome_index),
                   l_fragments,
                   l_frag_start_atom_index,
                   l_frag_stop_atom_index,
                   num_rotamers,
                   atom_stride,
                   num_atoms);

        // Calculate energy
        fp_type elect_total_trilinear = 0;
        fp_type emap_total_trilinear  = 0;
        fp_type dmap_total_trilinear  = 0;
        for (int atom_index = local_thread_id; atom_index < num_atoms; atom_index += thread_per_block) {
          fp_type coord_tex[3]{l_scratch_ligand_x[atom_index],
                               l_scratch_ligand_y[atom_index],
                               l_scratch_ligand_z[atom_index]};

          if (coord_tex[0] < minimum_x || coord_tex[0] > maximum_x || coord_tex[1] < minimum_y ||
              coord_tex[1] > maximum_y || coord_tex[2] < minimum_z || coord_tex[2] > maximum_z) {
            // Is outside
            const fp_type distance_two = powf(fabs(coord_tex[0] - center_x), fp_type{2}) +
                                         powf(fabs(coord_tex[1] - center_y), fp_type{2}) +
                                         powf(fabs(coord_tex[2] - center_z), fp_type{2});

            const fp_type epenalty = distance_two * ENERGYPENALTY;
            elect_total_trilinear += epenalty;
            emap_total_trilinear += epenalty;
          } else {
            // Is inside
            // Center atom coordinates on the grid center
            coord_tex[0] = (l_scratch_ligand_x[atom_index] - minimum_x) * inv_spacing,
            coord_tex[1] = (l_scratch_ligand_y[atom_index] - minimum_y) * inv_spacing;
            coord_tex[2] = (l_scratch_ligand_z[atom_index] - minimum_z) * inv_spacing;

            elect_total_trilinear +=
                trilinear_interpolation_cuda(coord_tex, electro_texture, index_n_x, index_n_xy) *
                l_ligand_charge[atom_index];
            dmap_total_trilinear +=
                trilinear_interpolation_cuda(coord_tex, desolv_texture, index_n_x, index_n_xy) *
                fabsf(l_ligand_charge[atom_index]);
            emap_total_trilinear +=
                trilinear_interpolation_cuda(coord_tex,
                                             atom_textures[l_atom_tex_indexes[atom_index]],
                                             index_n_x,
                                             index_n_xy);
          }
        }
        fp_type total_trilinear = elect_total_trilinear + dmap_total_trilinear + emap_total_trilinear;

        fp_type total_eintcal{0};
        if (num_rotamers > 0)
          total_eintcal += calc_intra_energy(l_scratch_ligand_x,
                                             l_scratch_ligand_y,
                                             l_scratch_ligand_z,
                                             l_ligand_vol,
                                             l_ligand_solpar,
                                             l_ligand_charge,
                                             l_ligand_num_hbond,
                                             l_ligand_Rij_hb,
                                             l_ligand_Rii,
                                             l_ligand_epsij_hb,
                                             l_ligand_epsii,
                                             num_nonbonds,
                                             l_ligand_nonbond_a1,
                                             l_ligand_nonbond_a2);

        // Reduction
        s_reduction_first[local_thread_id]  = total_trilinear;
        s_reduction_second[local_thread_id] = total_eintcal;
        for (int offset = thread_per_block / 2; offset > 0; offset /= 2) {
          if (local_thread_id < offset) {
            s_reduction_first[local_thread_id] += s_reduction_first[local_thread_id + offset];
            s_reduction_second[local_thread_id] += s_reduction_second[local_thread_id + offset];
          }
          __syncthreads(); // Synchronize threads after each reduction step
        }

        if (local_thread_id == 0) {
          total_trilinear                       = s_reduction_first[local_thread_id];
          total_eintcal                         = s_reduction_second[local_thread_id];
          const fp_type tors_free_energy        = num_rotamers * autodock_parameters::coeff_tors;
          s_chromosome_scores[chromosome_index] = total_trilinear + total_eintcal + tors_free_energy;
        }
      }

      // Generate the new population
      for (int chromosome_index = local_thread_id; chromosome_index < chromosome_number;
           chromosome_index += thread_per_block) {
        chromosome& next_chromosome = *(l_next_chromosomes + chromosome_index);

        const int best_individual_1 =
            tournament_selection_cuda(l_state, tournament_length, chromosome_number, s_chromosome_scores);
        const int best_individual_2 =
            tournament_selection_cuda(l_state, tournament_length, chromosome_number, s_chromosome_scores);

        // generate the offspring
        const int split_index = get_crossover_distribution(l_state, &num_rotamers);
        memcpy(next_chromosome.data(), &(l_chromosomes[best_individual_1][0]), split_index * sizeof(fp_type));
        const int parent2_copy_size = 6 + num_rotamers - split_index;
        if (parent2_copy_size > 0)
          memcpy(next_chromosome.data() + split_index,
                 &(l_chromosomes[best_individual_2][split_index]),
                 parent2_copy_size * sizeof(fp_type));

// mutate the offspring
#pragma unroll
        for (int i{0}; i < 3; ++i) {
          if (get_mutation_coin_distribution(l_state) < mutation_prob)
            next_chromosome[i] += get_mutation_change_distribution(l_state) * coordinate_step;
        }
#pragma unroll
        for (int i{3}; i < 6 + num_rotamers; ++i) {
          if (get_mutation_coin_distribution(l_state) < mutation_prob) {
            next_chromosome[i] += get_mutation_change_distribution(l_state) * angle_step;
          }
        }
      }

      // Swap the actual population with the next one
      // TODO check const
      chromosome* const tmp_chromosomes = l_chromosomes;
      l_chromosomes                     = l_next_chromosomes;
      l_next_chromosomes                = tmp_chromosomes;
      __syncthreads();
    }

    // Output

    // Compute the maximum value within the warp
    int min_index     = local_thread_id;
    fp_type min_score = s_chromosome_scores[min_index];
    for (int chromosome_index = local_thread_id + thread_per_block; chromosome_index < chromosome_number;
         chromosome_index += thread_per_block) {
      if (min_score > s_chromosome_scores[chromosome_index]) {
        min_index = chromosome_index;
        min_score = s_chromosome_scores[chromosome_index];
      }
    }

    // Reduction
    s_reduction_first[local_thread_id]  = min_score;
    s_reduction_second[local_thread_id] = min_index;
    for (int offset = thread_per_block / 2; offset > 0; offset /= 2) {
      if (local_thread_id < offset) {
        const fp_type other_min_score = s_reduction_first[local_thread_id + offset];
        const int other_min_index     = s_reduction_second[local_thread_id + offset];
        if (other_min_score < s_reduction_first[local_thread_id]) {
          s_reduction_first[local_thread_id]  = other_min_score;
          s_reduction_second[local_thread_id] = other_min_index;
        }
      }
      __syncthreads(); // Synchronize threads after each reduction step
    }
    min_index = s_reduction_second[0];
    if (local_thread_id == 0) {
      ligand_scores[ligand_id] = s_reduction_first[0];
      // memcpy(out_chromosome.data(), &(in_best_chromosome[0]), sizeof(fp_type) * (6 + num_rotamers));
    }
    chromosome& in_best_chromosome = l_chromosomes[min_index];
    chromosome& out_chromosome     = best_chromosomes[ligand_id];
    for (int gene_index = local_thread_id; gene_index < (6 + num_rotamers); gene_index += thread_per_block)
      // printf("%f\n", in_best_chromosome[gene_index]);
      out_chromosome[gene_index] = in_best_chromosome[gene_index];
  }

  void call_kernel(const int batch_ligands,
                   const int num_generations,
                   const int tournament_length,
                   const fp_type mutation_prob,
                   const int chromosome_number,
                   const int chromosome_stride,
                   const int atom_stride,
                   const int rotamers_stride,
                   const int nonbond_stride,
                   const fp_type* __restrict__ original_ligand_x,
                   const fp_type* __restrict__ original_ligand_y,
                   const fp_type* __restrict__ original_ligand_z,
                   fp_type* __restrict__ scratch_ligand_x,
                   fp_type* __restrict__ scratch_ligand_y,
                   fp_type* __restrict__ scratch_ligand_z,
                   const fp_type* __restrict__ ligand_vol,
                   const fp_type* __restrict__ ligand_solpar,
                   const fp_type* __restrict__ ligand_charge,
                   const int* __restrict__ ligand_num_hbond,
                   const fp_type* __restrict__ ligand_Rij_hb,
                   const fp_type* __restrict__ ligand_Rii,
                   const fp_type* __restrict__ ligand_epsij_hb,
                   const fp_type* __restrict__ ligand_epsii,
                   const int* __restrict__ ligand_num_nonbonds,
                   const int* __restrict__ ligand_nonbond_a1,
                   const int* __restrict__ ligand_nonbond_a2,
                   const int* __restrict__ ligand_num_atoms,
                   const int* __restrict__ ligand_num_rotamers,
                   const int* __restrict__ ligand_fragments,
                   const int* __restrict__ frag_start_atom_index,
                   const int* __restrict__ frag_stop_atom_index,
                   chromosome* __restrict__ chromosomes,
                   const fp_type minimum_x,
                   const fp_type minimum_y,
                   const fp_type minimum_z,
                   const fp_type maximum_x,
                   const fp_type maximum_y,
                   const fp_type maximum_z,
                   const fp_type center_x,
                   const fp_type center_y,
                   const fp_type center_z,
                   const int index_n_x,
                   const int index_n_xy,
                   const fp_type* const __restrict__* const __restrict__ atom_textures,
                   const int* __restrict__ atom_tex_indexes,
                   const fp_type* __restrict__ electro_texture,
                   const fp_type* __restrict__ desolv_texture,
                   XORWOWState* __restrict__ state,
                   fp_type* __restrict__ ligand_scores,
                   chromosome* __restrict__ best_chromosomes) {
    evaluate_fitness<<<batch_ligands, BLOCK_SIZE>>>(num_generations,
                                                    tournament_length,
                                                    mutation_prob,
                                                    chromosome_number,
                                                    chromosome_stride,
                                                    atom_stride,
                                                    rotamers_stride,
                                                    nonbond_stride,
                                                    original_ligand_x,
                                                    original_ligand_y,
                                                    original_ligand_z,
                                                    scratch_ligand_x,
                                                    scratch_ligand_y,
                                                    scratch_ligand_z,
                                                    ligand_vol,
                                                    ligand_solpar,
                                                    ligand_charge,
                                                    ligand_num_hbond,
                                                    ligand_Rij_hb,
                                                    ligand_Rii,
                                                    ligand_epsij_hb,
                                                    ligand_epsii,
                                                    ligand_num_nonbonds,
                                                    ligand_nonbond_a1,
                                                    ligand_nonbond_a2,
                                                    ligand_num_atoms,
                                                    ligand_num_rotamers,
                                                    ligand_fragments,
                                                    frag_start_atom_index,
                                                    frag_stop_atom_index,
                                                    chromosomes,
                                                    minimum_x,
                                                    minimum_y,
                                                    minimum_z,
                                                    maximum_x,
                                                    maximum_y,
                                                    maximum_z,
                                                    center_x,
                                                    center_y,
                                                    center_z,
                                                    index_n_x,
                                                    index_n_xy,
                                                    atom_textures,
                                                    atom_tex_indexes,
                                                    electro_texture,
                                                    desolv_texture,
                                                    state,
                                                    ligand_scores,
                                                    best_chromosomes);
    MUDOCK_CHECK_KERNELCALL();
    MUDOCK_CHECK(cudaDeviceSynchronize());
  }
} // namespace mudock
