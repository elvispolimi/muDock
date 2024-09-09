#pragma once

#include "chromosome.hpp"

#include <span>

namespace mudock {
  void apply(std::span<fp_type> x,
             std::span<fp_type> y,
             std::span<fp_type> z,
             const chromosome& c,
             const fragments<static_containers>& fragments);
}
