#pragma once

#if defined(MUDOCK_USE_HIP)
  #include <mudock/compute/genetic.hpp>
  #include <mudock/hip_implementation/geom_transform_hip.hpp>
  #include <mudock/hip_implementation/queue_hip.hpp>
  #include <mudock/hip_implementation/scoring_hip.hpp>
#endif
