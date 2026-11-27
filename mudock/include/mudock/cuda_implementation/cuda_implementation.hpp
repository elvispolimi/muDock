#pragma once

#if defined(MUDOCK_USE_CUDA)
  #include <mudock/cuda_implementation/scoring_cuda.cuh>
  #include <mudock/cuda_implementation/geom_transform_cuda.cuh>
  #include <mudock/cuda_implementation/queue_cuda.cuh>
  #include <mudock/compute/genetic.hpp>
#endif
