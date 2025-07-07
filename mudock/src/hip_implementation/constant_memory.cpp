#include <mudock/grid.hpp>
#include <mudock/hip_implementation/hip_check_error_macro.hpp>
#include <mudock/utils.hpp>

namespace mudock {
  __device__ __constant__ fp_type map_min_const[3];
  __device__ __constant__ fp_type map_max_const[3];
  __device__ __constant__ fp_type map_center_const[3];

  void setup_constant_memory(const point3D& minimum_coord,
                             const point3D& maximum_coord,
                             const point3D& center) {
    const auto l_map_min    = minimum_coord.get_component_p();
    const auto l_map_max    = maximum_coord.get_component_p();
    const auto l_map_center = center.get_component_p();

    MUDOCK_CHECK(hipMemcpyToSymbol(map_min_const, l_map_min, 3 * sizeof(fp_type), 0, hipMemcpyHostToDevice));
    MUDOCK_CHECK(hipMemcpyToSymbol(map_max_const, l_map_max, 3 * sizeof(fp_type), 0, hipMemcpyHostToDevice));
    MUDOCK_CHECK(
        hipMemcpyToSymbol(map_center_const, l_map_center, 3 * sizeof(fp_type), 0, hipMemcpyHostToDevice));
    MUDOCK_CHECK(hipDeviceSynchronize());
  }
} // namespace mudock
