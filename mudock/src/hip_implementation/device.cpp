#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/chem/autodock_protein.hpp>
#include <mudock/grid/mdindex.hpp>
#include <mudock/grid/space_grid.hpp>
#include <mudock/hip_implementation/constant_memory.hpp>
#include <mudock/hip_implementation/device.hpp>
#include <mudock/hip_implementation/hip_check_error_macro.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {
  // TODO pack together maps using float4 data structures
  void init_texture_memory(const space_grid_view<const fp_type> grid_map,
                           hipTexture_wrapper& tex_obj,
                           const hipStream_t& stream) {
    // Create 3D CUDA array for the texture
    hipArray_t d_array;
    hipChannelFormatDesc channel_desc = hipCreateChannelDesc<fp_type>();

    const md_index<3> map_index = grid_map.get_space_index();
    hipExtent extent            = make_hipExtent(map_index.size_x(), map_index.size_y(), map_index.size_z());
    MUDOCK_CHECK(hipMalloc3DArray(&d_array, &channel_desc, extent, hipArrayDefault));

    // Copy data from host to the 3D hip array
    hipMemcpy3DParms copyParams = {0};
    copyParams.srcPtr           = make_hipPitchedPtr((void*) grid_map.data(),
                                           extent.width * sizeof(fp_type),
                                           extent.width,
                                           extent.height);
    copyParams.srcArray         = nullptr;
    copyParams.srcPos           = make_hipPos(0, 0, 0);
    copyParams.dstArray         = d_array;
    copyParams.dstPtr           = hipPitchedPtr{};
    copyParams.dstPos           = make_hipPos(0, 0, 0);
    copyParams.extent           = extent;
    copyParams.kind             = hipMemcpyHostToDevice;
    MUDOCK_CHECK(hipMemcpy3DAsync(&copyParams, stream));

    // Create texture object
    hipResourceDesc res_desc;
    memset(&res_desc, 0, sizeof(hipResourceDesc));
    res_desc.resType         = hipResourceTypeArray;
    res_desc.res.array.array = d_array;

    hipTextureDesc tex_desc;
    memset(&tex_desc, 0, sizeof(hipTextureDesc));
    tex_desc.addressMode[0] = hipAddressModeClamp;
    tex_desc.addressMode[1] = hipAddressModeClamp;
    tex_desc.addressMode[2] = hipAddressModeClamp;
#ifdef MUDOCK_TEST
    tex_desc.filterMode = hipFilterModePoint;
#else
    tex_desc.filterMode = hipFilterModeLinear;
#endif
    tex_desc.readMode         = hipReadModeElementType;
    tex_desc.normalizedCoords = false;

    MUDOCK_CHECK(hipCreateTextureObject(&tex_obj(), &res_desc, &tex_desc, NULL));
  }

  device::device(const std::size_t gpu_id, const autodock_protein& adt_protein)
      : id(gpu_id),
        center_maps(adt_protein.get_center()),
        stream(get_stream()),
        atom_tex(stream()),
        adt_protein(adt_protein) {
    // Allocate grid maps
    const point<fp_type, 3>&minimum(adt_protein.get_min()), maximum(adt_protein.get_max()),
        center(adt_protein.get_center());

    atom_tex.alloc(num_autodock_grids());
    for (size_t index{0}; index < num_autodock_grids(); index++) {
      auto& tex = atom_tex.host_pointer()[index];

      const auto& grid = adt_protein.get_atom_map(static_cast<autodock_grid_type>(index));
      init_texture_memory(grid, tex, stream());
      MUDOCK_CHECK(hipStreamSynchronize(stream()));
    }

    atom_tex.copy_host2device();

    // Grid spacing fixed to 0.5 Angstrom
    setup_constant_memory(minimum, maximum, center);
    MUDOCK_CHECK(hipStreamSynchronize(stream()));

    MUDOCK_CHECK(hipDeviceGetAttribute(&wavefront_size, hipDeviceAttributeWarpSize, id));
  }

  hipStream_t device::get_stream() const {
    hipStream_t cs;
    MUDOCK_CHECK(hipSetDevice(static_cast<int>(id)));
    MUDOCK_CHECK(hipStreamCreate(&cs););
    return cs;
  }
} // namespace mudock
