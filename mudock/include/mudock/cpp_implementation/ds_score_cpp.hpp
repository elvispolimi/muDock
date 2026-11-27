#pragma once

#include <mudock/compute/drugscore/ds_score.hpp>
#include <mudock/compute/drugscore/ds_score_kernel.hpp>
#include <mudock/cpp_implementation/queue_cpp.hpp>

namespace mudock {

  template<>
  void ds_score_kernel<queue_cpp>::operator()();

} // namespace mudock
