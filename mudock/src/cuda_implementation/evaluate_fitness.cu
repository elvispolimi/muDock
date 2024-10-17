#include <cstdio>
#include <curand_kernel.h>
#include <mudock/cuda_implementation/calc_energy.cuh>
#include <mudock/cuda_implementation/cuda_check_error_macro.cuh>
#include <mudock/cuda_implementation/evaluate_fitness.cuh>
#include <mudock/cuda_implementation/mutate.cuh>
#include <mudock/grid.hpp>
#include <mudock/utils.hpp>

namespace mudock {
  __device__ __constant__ fp_type map_min_const[3];
  __device__ __constant__ fp_type map_max_const[3];
  __device__ __constant__ fp_type map_center_const[3];

  void setup_constant_memory(const point3D& minimum_coord,
                             const point3D& maximum_coord,
                             const point3D& center) {
    const fp_type l_map_min[3]{minimum_coord.x, minimum_coord.y, minimum_coord.z};
    const fp_type l_map_max[3]{maximum_coord.x, maximum_coord.y, maximum_coord.z};
    const fp_type l_map_center[3]{center.x, center.y, center.z};

    MUDOCK_CHECK(
        cudaMemcpyToSymbol(map_min_const, &l_map_min, 3 * sizeof(fp_type), 0, cudaMemcpyHostToDevice));
    MUDOCK_CHECK(
        cudaMemcpyToSymbol(map_max_const, &l_map_max, 3 * sizeof(fp_type), 0, cudaMemcpyHostToDevice));
    MUDOCK_CHECK(
        cudaMemcpyToSymbol(map_center_const, &l_map_center, 3 * sizeof(fp_type), 0, cudaMemcpyHostToDevice));
  }

  __device__ fp_type trilinear_interpolation_cuda(const fp_type coord[], const cudaTextureObject_t& tex) {
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
          const fp_type tmp = tex3D<fp_type>(tex, u0 + n, v0 + t, w0 + i);
          value += pu[n] * pv[t] * pw[i] * tmp;
        }
    return value;
  }

  template<typename T>
  __device__ const T random_gen_cuda(curandState& state, const T min, const T max) {
    fp_type value;
    if constexpr (is_debug()) {
      // TODO value here for debug
      value = fp_type{0.4};
    } else {
      value = curand_uniform(&state);
    }
    return static_cast<T>((value * static_cast<fp_type>(max - min)) + min);
  }

  __device__ int get_selection_distribution(curandState& state, const int* population_number) {
    return random_gen_cuda<int>(state, 0, *population_number - 1);
  };

  __device__ fp_type get_init_change_distribution(curandState& state) {
    return random_gen_cuda<fp_type>(state, -45, 45);
  }
  __device__ fp_type get_mutation_change_distribution(curandState& state) {
    return random_gen_cuda<fp_type>(state, -10, 10);
  };
  __device__ fp_type get_mutation_coin_distribution(curandState& state) {
    return random_gen_cuda<fp_type>(state, 0, 1);
  };
  __device__ int get_crossover_distribution(curandState& state, const int* num_rotamers) {
    return random_gen_cuda<int>(state, 0, 6 + *num_rotamers);
  };

  __device__ int tournament_selection_cuda(curandState& state,
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
} // namespace mudock
