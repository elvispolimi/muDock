#include <mudock/alpaka_implementation/alpaka_dot_product.hpp>

#include <alpaka/alpaka.hpp>
#include <alpaka/atomic/Traits.hpp>
#include <alpaka/example/ExampleDefaultAcc.hpp>

#include <cstddef>
#include <stdexcept>

namespace mudock::alpaka_demo {
  namespace {
    struct dot_product_kernel {
      ALPAKA_NO_HOST_ACC_WARNING
      template<typename TAcc>
      ALPAKA_FN_ACC auto operator()(TAcc const& acc,
                                    float const* lhs,
                                    float const* rhs,
                                    double* result,
                                    std::size_t num_elements) const -> void {
        static_assert(alpaka::Dim<TAcc>::value == 1u,
                      "dot_product_kernel expects a 1D Alpaka accelerator");

        for(auto i : alpaka::uniformElements(acc, num_elements)) {
          alpaka::atomicAdd(
            acc,
            result,
            static_cast<double>(lhs[i]) * static_cast<double>(rhs[i]),
            alpaka::hierarchy::Grids{});
        }
      }
    };
  } // namespace

  double dot_product(std::span<const float> lhs, std::span<const float> rhs) {
    if(lhs.size() != rhs.size()) {
      throw std::invalid_argument("dot_product requires vectors with the same size");
    }
    if(lhs.empty()) {
      return 0.0F;
    }

    using Dim = alpaka::DimInt<1u>;
    using Idx = std::size_t;
    using Acc = alpaka::ExampleDefaultAcc<Dim, Idx>;
    using DevAcc = alpaka::Dev<Acc>;
    using QueueAcc = alpaka::Queue<Acc, alpaka::Blocking>;
    using DevHost = alpaka::DevCpu;
    using BufAcc = alpaka::Buf<DevAcc, double, Dim, Idx>;

    auto const host_platform = alpaka::PlatformCpu{};
    auto const dev_host = alpaka::getDevByIdx(host_platform, 0);

    auto const acc_platform = alpaka::Platform<Acc>{};
    auto const dev_acc = alpaka::getDevByIdx(acc_platform, 0);
    QueueAcc queue(dev_acc);

    auto const extent = alpaka::Vec<Dim, Idx>{static_cast<Idx>(lhs.size())};
    auto const scalar_extent = alpaka::Vec<Dim, Idx>{Idx{1}};

    auto lhs_host = alpaka::createView(dev_host, lhs.data(), extent);
    auto rhs_host = alpaka::createView(dev_host, rhs.data(), extent);
    double result_value = 0.0;
    auto result_host = alpaka::createView(dev_host, &result_value, scalar_extent);

    auto lhs_acc = alpaka::allocBuf<float, Idx>(dev_acc, extent);
    auto rhs_acc = alpaka::allocBuf<float, Idx>(dev_acc, extent);
    BufAcc result_acc(alpaka::allocBuf<double, Idx>(dev_acc, scalar_extent));

    alpaka::memcpy(queue, lhs_acc, lhs_host);
    alpaka::memcpy(queue, rhs_acc, rhs_host);
    alpaka::memcpy(queue, result_acc, result_host);

    dot_product_kernel kernel;
    alpaka::KernelCfg<Acc> const kernel_cfg = {extent, Idx{256}};
    auto const work_div = alpaka::getValidWorkDiv(
      kernel_cfg,
      dev_acc,
      kernel,
      alpaka::getPtrNative(lhs_acc),
      alpaka::getPtrNative(rhs_acc),
      alpaka::getPtrNative(result_acc),
      static_cast<Idx>(lhs.size()));

    auto const task = alpaka::createTaskKernel<Acc>(
      work_div,
      kernel,
      alpaka::getPtrNative(lhs_acc),
      alpaka::getPtrNative(rhs_acc),
      alpaka::getPtrNative(result_acc),
      static_cast<Idx>(lhs.size()));

    alpaka::enqueue(queue, task);
    alpaka::memcpy(queue, result_host, result_acc);
    alpaka::wait(queue);

    return result_value;
  }

} // namespace mudock::alpaka_demo
