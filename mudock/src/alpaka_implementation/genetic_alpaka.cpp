#include <mudock/alpaka_implementation/genetic_alpaka.hpp>

#include <stdexcept>

namespace mudock {
  template<>
  void genetic_kernel<queue_alpaka>::initialize() {
    throw std::runtime_error("Alpaka genetic initialize kernel is not implemented yet");
  }

  template<>
  void genetic_kernel<queue_alpaka>::operator()() {
    throw std::runtime_error("Alpaka genetic iterate kernel is not implemented yet");
  }

  template<>
  void genetic_kernel<queue_alpaka>::finalize() {
    throw std::runtime_error("Alpaka genetic finalize kernel is not implemented yet");
  }
} // namespace mudock
