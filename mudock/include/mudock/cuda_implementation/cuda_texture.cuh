#pragma once

#include <cuda_runtime.h>
#include <memory>
#include <mudock/cuda_implementation/cuda_utils.cuh>
#include <mudock/type_alias.hpp>

namespace mudock {
  struct cuda_texture_wrapper {
    cuda_texture_wrapper(const int dev_id_,
                         const int width,
                         const int height,
                         const int depth,
                         const fp_type* grid_map)
        : dev_id(dev_id_) {
      MUDOCK_CHECK(cudaSetDevice(dev_id));
      // Create 3D CUDA array for the texture
      cudaChannelFormatDesc channel_desc = cudaCreateChannelDesc<fp_type>();

      cudaExtent extent = make_cudaExtent(width, height, depth);
      MUDOCK_CHECK(cudaMalloc3DArray(&array_, &channel_desc, extent));

      // Copy data from host to the 3D CUDA array
      cudaMemcpy3DParms copyParams = {0};
      copyParams.srcPtr =
          make_cudaPitchedPtr((void*) grid_map, extent.width * sizeof(fp_type), extent.width, extent.height);
      copyParams.srcArray = nullptr;
      copyParams.srcPos   = make_cudaPos(0, 0, 0);
      copyParams.dstArray = array_;
      copyParams.dstPtr   = cudaPitchedPtr{};
      copyParams.dstPos   = make_cudaPos(0, 0, 0);
      copyParams.extent   = extent;
      copyParams.kind     = cudaMemcpyHostToDevice;
      MUDOCK_CHECK(cudaMemcpy3DAsync(&copyParams));

      // Create texture object
      cudaResourceDesc res_desc;
      memset(&res_desc, 0, sizeof(cudaResourceDesc));
      res_desc.resType         = cudaResourceTypeArray;
      res_desc.res.array.array = array_;

      cudaTextureDesc tex_desc;
      memset(&tex_desc, 0, sizeof(cudaTextureDesc));
      tex_desc.addressMode[0] = cudaAddressModeClamp;
      tex_desc.addressMode[1] = cudaAddressModeClamp;
      tex_desc.addressMode[2] = cudaAddressModeClamp;
      tex_desc.filterMode = cudaFilterModePoint;
      tex_desc.readMode         = cudaReadModeElementType;
      tex_desc.normalizedCoords = false;

      MUDOCK_CHECK(cudaCreateTextureObject(&tex_, &res_desc, &tex_desc, NULL));
    };

    // Keep the destructor not throwable because if it is a global variable it gets destroyed after the driver deallocated everything
    ~cuda_texture_wrapper() noexcept {
      cudaSetDevice(dev_id);
      // TODO fix the try catch
      if (tex_) {
        cudaDestroyTextureObject(tex_);
      }
      if (array_) {
        cudaFreeArray(array_);
      }
    }

    // forbid copying, allow moving if you want
    cuda_texture_wrapper(const cuda_texture_wrapper&)            = delete;
    cuda_texture_wrapper& operator=(const cuda_texture_wrapper&) = delete;

    cuda_texture_wrapper(cuda_texture_wrapper&& other)            = delete;
    cuda_texture_wrapper& operator=(cuda_texture_wrapper&& other) = delete;

    cudaTextureObject_t tex_ = 0;
    cudaArray_t array_       = nullptr;
    const int dev_id         = 0;
  };

  struct cuda_texture_devices {
    std::vector<std::unique_ptr<cuda_texture_wrapper>> textures;
    cudaTextureObject_t* textures_dev_p = nullptr;
    const int dev_id                    = 0;

    cuda_texture_devices(const int dev_id_,
                         const int map_index_x,
                         const int map_index_xy,
                         const int map_index_xyz,
                         const fp_type* grid_map,
                         const int num_maps)
        : dev_id(dev_id_) {
      MUDOCK_CHECK(cudaSetDevice(dev_id));
      std::vector<cudaTextureObject_t> tex_host_p;
      for (int index{0}; index < num_maps; index++) {
        textures.emplace_back(std::make_unique<cuda_texture_wrapper>(dev_id_,
                                                                     map_index_x,
                                                                     map_index_xy / map_index_x,
                                                                     map_index_xyz / map_index_xy,
                                                                     grid_map + map_index_xyz * index));
        tex_host_p.push_back((*textures.back()).tex_);
      }

      MUDOCK_CHECK(cudaMalloc(&textures_dev_p, sizeof(cudaTextureObject_t) * num_maps));
      MUDOCK_CHECK(cudaMemcpyAsync(textures_dev_p,
                                   tex_host_p.data(),
                                   num_maps * sizeof(cudaTextureObject_t),
                                   cudaMemcpyHostToDevice));
      MUDOCK_CHECK(cudaDeviceSynchronize());
    };

    // no copy
    cuda_texture_devices(const cuda_texture_devices&)            = delete;
    cuda_texture_devices& operator=(const cuda_texture_devices&) = delete;

    // no move
    cuda_texture_devices(cuda_texture_devices&& other)            = delete;
    cuda_texture_devices& operator=(cuda_texture_devices&& other) = delete;

    // Keep it noexcept see upper destructor
    ~cuda_texture_devices() noexcept {
      if (textures_dev_p) {
        cudaSetDevice(dev_id);
        cudaFree(textures_dev_p);
        textures_dev_p = nullptr;
      }
    };
  };
} // namespace mudock
