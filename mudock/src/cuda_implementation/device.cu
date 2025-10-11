#include "mudock/cuda_implementation/cuda_texture.cuh"

#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/chem/autodock_protein.hpp>
#include <mudock/cuda_implementation/cuda_check_error_macro.cuh>
#include <mudock/cuda_implementation/device.cuh>
#include <mudock/cuda_implementation/evaluate_fitness.cuh>
#include <mudock/grid/mdindex.hpp>
#include <mudock/grid/space_grid.hpp>
#include <mudock/type_alias.hpp>

namespace mudock {
  // TODO pack together maps using float4 data structures
  void init_texture_memory(const space_grid_view<const fp_type> grid_map,
                           cudaTexture_wrapper& tex_obj,
                           const cudaStream_t& stream) {
    // Create 3D CUDA array for the texture
    cudaArray* d_array;
    cudaChannelFormatDesc channel_desc = cudaCreateChannelDesc<fp_type>();

    const md_index<3> map_index = grid_map.get_space_index();
    cudaExtent extent           = make_cudaExtent(map_index.size_x(), map_index.size_y(), map_index.size_z());
    MUDOCK_CHECK(cudaMalloc3DArray(&d_array, &channel_desc, extent));

    // Copy data from host to the 3D CUDA array
    cudaMemcpy3DParms copyParams = {0};
    copyParams.srcPtr            = make_cudaPitchedPtr((void*) grid_map.data(),
                                            extent.width * sizeof(fp_type),
                                            extent.width,
                                            extent.height);
    copyParams.srcArray          = nullptr;
    copyParams.srcPos            = make_cudaPos(0, 0, 0);
    copyParams.dstArray          = d_array;
    copyParams.dstPtr            = cudaPitchedPtr{};
    copyParams.dstPos            = make_cudaPos(0, 0, 0);
    copyParams.extent            = extent;
    copyParams.kind              = cudaMemcpyHostToDevice;
    MUDOCK_CHECK(cudaMemcpy3DAsync(&copyParams, stream));

    // Create texture object
    cudaResourceDesc res_desc;
    memset(&res_desc, 0, sizeof(cudaResourceDesc));
    res_desc.resType         = cudaResourceTypeArray;
    res_desc.res.array.array = d_array;

    cudaTextureDesc tex_desc;
    memset(&tex_desc, 0, sizeof(cudaTextureDesc));
    tex_desc.addressMode[0] = cudaAddressModeClamp;
    tex_desc.addressMode[1] = cudaAddressModeClamp;
    tex_desc.addressMode[2] = cudaAddressModeClamp;
#ifdef MUDOCK_TEST
    tex_desc.filterMode = cudaFilterModePoint;
#else
    tex_desc.filterMode = cudaFilterModeLinear;
#endif
    tex_desc.readMode         = cudaReadModeElementType;
    tex_desc.normalizedCoords = false;

    MUDOCK_CHECK(cudaCreateTextureObject(&tex_obj(), &res_desc, &tex_desc, NULL));
  }

  device::device(const std::size_t gpu_id, const autodock_protein& adt_protein)
      : id(gpu_id),
        center_maps(adt_protein.get_center()),
        adt_protein(adt_protein),
        stream(get_stream()),
        atom_tex(stream()) {
    // Allocate grid maps
    const point<fp_type, 3>&minimum(adt_protein.get_min()), maximum(adt_protein.get_max()),
        center(adt_protein.get_center());

    atom_tex.alloc(num_autodock_grids());
    for (size_t index{0}; index < num_autodock_grids(); index++) {
      auto& tex = atom_tex.host_pointer()[index];

      const auto& grid = adt_protein.get_atom_map(static_cast<autodock_grid_type>(index));
      init_texture_memory(grid, tex, stream());
      MUDOCK_CHECK(cudaStreamSynchronize(stream()));
    }

    atom_tex.copy_host2device();

    // Grid spacing fixed to 0.5 Angstrom
    setup_constant_memory(minimum, maximum, center);
    MUDOCK_CHECK(cudaStreamSynchronize(stream()));
  }

  cudaStream_wrapper device::get_stream() const { return cudaStream_wrapper(id); }

} // namespace mudock
