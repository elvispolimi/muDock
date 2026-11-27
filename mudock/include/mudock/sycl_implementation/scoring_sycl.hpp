#pragma once

// SYCL scorer implementations used by generic pipeline dispatch belong here.
// genetic_sycl.cpp includes this aggregate and resolves the implementation
// through scoring_backend_t::device_kernel, so adding a scorer does not add a
// scorer-specific dependency to the genetic backend. The pipeline registry
// still controls which scoring functions are instantiated and dispatched.
#include <mudock/sycl_implementation/adt_score_sycl.hpp>
