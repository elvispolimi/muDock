#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <mudock/alpaka_implementation/queue_alpaka.hpp>
#include <mudock/compute/buffer.hpp>
#include <mudock/implementation_types.hpp>
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

  mudock::buffer_vector<int, mudock::queue_alpaka> device_buffer{queue};
  device_buffer.alloc(host.size());
  std::copy(host.begin(), host.end(), device_buffer.host_pointer());
  device_buffer.copy_host2device();
  queue->synchronize();
  device_buffer.copy_device2host();
  queue->synchronize();
  std::copy(device_buffer.host_pointer(),
            device_buffer.host_pointer() + device_buffer.num_elements(),
            roundtrip.begin());

  if (roundtrip != host) {
    std::cerr << "Alpaka buffer host/device roundtrip failed" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Alpaka queue/buffer smoke test passed" << std::endl;
  return EXIT_SUCCESS;
#endif
}
