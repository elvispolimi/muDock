#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/chem/autodock_protein.hpp>
#include <mudock/grid/mdindex.hpp>
#include <mudock/grid/space_grid.hpp>
#include <mudock/hip_implementation/constant_memory.hpp>
#include <mudock/hip_implementation/device.hpp>
#include <mudock/hip_implementation/hip_check_error_macro.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {
  device::device(const std::size_t gpu_id, const autodock_protein& adt_protein)
      : id(gpu_id),
        center_maps(adt_protein.get_center()),
        stream(get_stream()),
        adt_protein(adt_protein),
        atom_tex(stream(), adt_protein) {
    // Allocate grid maps
    const point<fp_type, 3>&minimum(adt_protein.get_min()), maximum(adt_protein.get_max()),
        center(adt_protein.get_center());

    // Grid spacing fixed to 0.5 Angstrom
    setup_constant_memory(minimum, maximum, center);
    MUDOCK_CHECK(hipStreamSynchronize(stream()));

    MUDOCK_CHECK(hipDeviceGetAttribute(&wavefront_size, hipDeviceAttributeWarpSize, id));
  }

  hipStream_wrapper device::get_stream() const { return {static_cast<int>(id)}; }
} // namespace mudock
