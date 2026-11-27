#pragma once

#if defined(MUDOCK_USE_SYCL)
  #include <mudock/compute/genetic.hpp>
  #include <mudock/sycl_implementation/geom_transform_sycl.hpp>
  #include <mudock/sycl_implementation/queue_sycl.hpp>
  #include <mudock/sycl_implementation/scoring_sycl.hpp>
#endif
