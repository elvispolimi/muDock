#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mudock/alpaka_implementation/alpaka_dot_product.hpp>
#include <mudock/alpaka_implementation/queue_alpaka.hpp>
#include <mudock/compute/buffer.hpp>
#include <vector>

int main() {
  using clock                          = std::chrono::steady_clock;
  constexpr std::size_t vector_size    = 1u << 24;
  constexpr int dot_product_iterations = 96;
  constexpr float lhs_value            = 1.5F;
  constexpr float rhs_value            = 2.0F;

  std::vector<float> lhs(vector_size, lhs_value);
  std::vector<float> rhs(vector_size, rhs_value);
  {
    auto queue = std::make_shared<mudock::queue_alpaka>(0, mudock::device_type::CPU);
    mudock::buffer_vector<float, mudock::queue_alpaka> lhs_dev{queue};
    mudock::buffer_vector<float, mudock::queue_alpaka> rhs_dev{queue};
    mudock::buffer_vector<float, mudock::queue_alpaka> copy_dev{queue};
    std::vector<float> copied(lhs.size(), 0.0F);

    lhs_dev.alloc(lhs.size());
    rhs_dev.alloc(rhs.size());
    copy_dev.alloc(lhs.size());
    std::copy(lhs.begin(), lhs.end(), lhs_dev.host_pointer());
    std::copy(rhs.begin(), rhs.end(), rhs_dev.host_pointer());
    lhs_dev.copy_host2device();
    rhs_dev.copy_host2device();
    copy_dev.copy_device2device(lhs_dev);
    queue->synchronize();
    copy_dev.copy_device2host();
    queue->synchronize();
    std::copy(copy_dev.host_pointer(), copy_dev.host_pointer() + copy_dev.num_elements(), copied.begin());

    if (copied != lhs) {
      std::cerr << "Alpaka buffer device-to-device copy failed" << std::endl;
      return EXIT_FAILURE;
    }
  }

  double const expected_dot =
      static_cast<double>(vector_size) * static_cast<double>(lhs_value) * static_cast<double>(rhs_value);
  double result        = 0.0;
  auto const dot_begin = clock::now();
  for (int i = 0; i < dot_product_iterations; ++i) { result = mudock::alpaka_demo::dot_product(lhs, rhs); }
  auto const dot_end     = clock::now();
  auto const dot_seconds = std::chrono::duration<double>(dot_end - dot_begin).count();

  std::cout << "alpaka dot product result: " << result << std::endl;
  std::cout << "alpaka dot product runtime over " << dot_product_iterations << " launches: " << std::fixed
            << std::setprecision(3) << dot_seconds << " s" << std::endl;

  if (std::fabs(result - expected_dot) > 1.0e-9 * expected_dot) {
    std::cerr << "Expected dot product " << expected_dot << " but got " << result << std::endl;
    return EXIT_FAILURE;
  }

  double reduction_result    = 0.0;
  auto const reduction_begin = clock::now();
  for (int i = 0; i < dot_product_iterations; ++i) {
    reduction_result = mudock::alpaka_demo::dot_product_reduction(lhs, rhs);
  }
  auto const reduction_end     = clock::now();
  auto const reduction_seconds = std::chrono::duration<double>(reduction_end - reduction_begin).count();

  std::cout << "alpaka reduction dot product result: " << reduction_result << std::endl;
  std::cout << "alpaka reduction dot product runtime over " << dot_product_iterations
            << " launches: " << std::fixed << std::setprecision(3) << reduction_seconds << " s" << std::endl;

  if (std::fabs(reduction_result - expected_dot) > 1.0e-9 * expected_dot) {
    std::cerr << "Expected reduction dot product " << expected_dot << " but got " << reduction_result
              << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "alpaka dot product and dot-product reduction tests passed" << std::endl;
  return EXIT_SUCCESS;
}
