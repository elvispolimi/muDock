#include <mudock/alpaka_implementation/object_alpaka.hpp>
#include <mudock/alpaka_implementation/queue_alpaka.hpp>
#include <mudock/implementation_types.hpp>

#include <cstdlib>
#include <iostream>
#include <memory>
#include <vector>

int main() {
#ifndef MUDOCK_USE_ALPAKA
  return EXIT_FAILURE;
#else
  const auto impl = mudock::get_impl_type("ALPAKA");
  if (impl != mudock::implementation_type::ALPAKA) {
    return EXIT_FAILURE;
  }

  auto queue = std::make_shared<mudock::queue_alpaka>(0, mudock::device_type::CPU);
  (*queue)();

  std::vector<int> host{1, 2, 3, 4, 5, 6, 7, 8};
  std::vector<int> roundtrip(host.size(), 0);

  mudock::object<int, mudock::queue_alpaka> device_object{queue};
  device_object.alloc(host.size());
  device_object.copy_host2device(host.data());
  queue->synchronize();
  device_object.copy_device2host(roundtrip.data());
  queue->synchronize();

  if (roundtrip != host) {
    std::cerr << "Alpaka object host/device roundtrip failed" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Alpaka queue/object smoke test passed" << std::endl;
  return EXIT_SUCCESS;
#endif
}
