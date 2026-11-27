#pragma once

// CUDA scorer implementations used by generic pipeline dispatch belong here.
// genetic_cuda.cu includes this aggregate and resolves the implementation
// through scoring_backend_t::device_kernel, so adding a scorer does not add a
// scorer-specific dependency to the genetic backend. The pipeline registry
// still controls which scoring functions are instantiated and dispatched.
#include <mudock/cuda_implementation/adt_score_cuda.cuh>
