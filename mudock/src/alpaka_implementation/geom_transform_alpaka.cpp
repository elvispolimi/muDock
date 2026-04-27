#include <mudock/alpaka_implementation/geom_transform_alpaka.hpp>

#include <stdexcept>

namespace mudock {
  template<>
  void geom_kernel<queue_alpaka>::operator()() {
    throw std::runtime_error("Alpaka geometric transform kernel is not implemented yet");
  }
} // namespace mudock
