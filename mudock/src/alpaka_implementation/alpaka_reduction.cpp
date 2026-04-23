#include <mudock/alpaka_implementation/alpaka_dot_product.hpp>

#include <alpaka/alpaka.hpp>
#include <alpaka/block/shared/dyn/Traits.hpp>
#include <alpaka/block/sync/Traits.hpp>
#include <alpaka/example/ExampleDefaultAcc.hpp>

#include <cstddef>
#include <stdexcept>

namespace mudock::alpaka_demo {
  struct dot_product_partial_reduction_kernel {
      ALPAKA_NO_HOST_ACC_WARNING
      template<typename TAcc>
      ALPAKA_FN_ACC auto operator()(TAcc const& acc,
                                    float const* lhs,
                                    float const* rhs,
                                    double* partial_sums,
                                    std::size_t num_elements) const -> void {
        static_assert(alpaka::Dim<TAcc>::value == 1u,
                      "dot_product_partial_reduction_kernel expects a 1D Alpaka accelerator");

        auto const block_thread_idx = alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u];
        auto const block_idx = alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u];
        auto const block_thread_count = alpaka::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0u];
        double* const shared = alpaka::getDynSharedMem<double>(acc);

        double thread_sum = 0.0;
        for(auto i : alpaka::uniformElements(acc, num_elements)) {
          thread_sum += static_cast<double>(lhs[i]) * static_cast<double>(rhs[i]);
        }

        shared[block_thread_idx] = thread_sum;
        alpaka::syncBlockThreads(acc);

        for(auto stride = block_thread_count / 2u; stride > 0u; stride /= 2u) {
          if(block_thread_idx < stride) {
            shared[block_thread_idx] += shared[block_thread_idx + stride];
          }
          alpaka::syncBlockThreads(acc);
        }

        if(block_thread_idx == 0u) {
          partial_sums[block_idx] = shared[0];
        }
      }
    };

    struct final_reduction_kernel {
      ALPAKA_NO_HOST_ACC_WARNING
      template<typename TAcc>
      ALPAKA_FN_ACC auto operator()(TAcc const& acc,
                                    double const* partial_sums,
                                    double* result,
                                    std::size_t num_partials) const -> void {
        static_assert(alpaka::Dim<TAcc>::value == 1u,
                      "final_reduction_kernel expects a 1D Alpaka accelerator");

        auto const block_thread_idx = alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u];
        auto const block_thread_count = alpaka::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0u];
        double* const shared = alpaka::getDynSharedMem<double>(acc);

        double thread_sum = 0.0;
        for(std::size_t i = block_thread_idx; i < num_partials; i += block_thread_count) {
          thread_sum += partial_sums[i];
        }

        shared[block_thread_idx] = thread_sum;
        alpaka::syncBlockThreads(acc);

        for(auto stride = block_thread_count / 2u; stride > 0u; stride /= 2u) {
          if(block_thread_idx < stride) {
            shared[block_thread_idx] += shared[block_thread_idx + stride];
          }
          alpaka::syncBlockThreads(acc);
        }

        if(block_thread_idx == 0u) {
          result[0] = shared[0];
        }
      }
    };

  double dot_product_reduction(std::span<const float> lhs, std::span<const float> rhs) {
    if(lhs.size() != rhs.size()) {
      throw std::invalid_argument("dot_product_reduction requires vectors with the same size");
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
    auto const acc_props = alpaka::getAccDevProps<Acc>(dev_acc);

    auto const extent = alpaka::Vec<Dim, Idx>{static_cast<Idx>(lhs.size())};
    auto const scalar_extent = alpaka::Vec<Dim, Idx>{Idx{1}};
    Idx const block_size = std::min<Idx>(static_cast<Idx>(256u), static_cast<Idx>(acc_props.m_blockThreadCountMax));

    auto lhs_host = alpaka::createView(dev_host, lhs.data(), extent);
    auto rhs_host = alpaka::createView(dev_host, rhs.data(), extent);
    double result_value = 0.0;
    auto result_host = alpaka::createView(dev_host, &result_value, scalar_extent);

    auto lhs_acc = alpaka::allocBuf<float, Idx>(dev_acc, extent);
    auto rhs_acc = alpaka::allocBuf<float, Idx>(dev_acc, extent);
    BufAcc result_acc(alpaka::allocBuf<double, Idx>(dev_acc, scalar_extent));
    auto const partial_count =
      static_cast<Idx>((lhs.size() + static_cast<std::size_t>(block_size) - 1u) /
                       static_cast<std::size_t>(block_size));
    auto const partial_extent = alpaka::Vec<Dim, Idx>{partial_count};
    BufAcc partial_acc(alpaka::allocBuf<double, Idx>(dev_acc, partial_extent));

    alpaka::memcpy(queue, lhs_acc, lhs_host);
    alpaka::memcpy(queue, rhs_acc, rhs_host);
    alpaka::memcpy(queue, result_acc, result_host);

    dot_product_partial_reduction_kernel partial_kernel;
    auto const partial_work_div =
      alpaka::WorkDivMembers<Dim, Idx>(alpaka::Vec<Dim, Idx>{partial_count},
                                       alpaka::Vec<Dim, Idx>{block_size},
                                       alpaka::Vec<Dim, Idx>{Idx{1}});
    auto const partial_task = alpaka::createTaskKernel<Acc>(
      partial_work_div,
      partial_kernel,
      alpaka::getPtrNative(lhs_acc),
      alpaka::getPtrNative(rhs_acc),
      alpaka::getPtrNative(partial_acc),
      static_cast<Idx>(lhs.size()));
    alpaka::enqueue(queue, partial_task);

    final_reduction_kernel final_kernel;
    auto const final_work_div =
      alpaka::WorkDivMembers<Dim, Idx>(alpaka::Vec<Dim, Idx>{Idx{1}},
                                       alpaka::Vec<Dim, Idx>{block_size},
                                       alpaka::Vec<Dim, Idx>{Idx{1}});
    auto const final_task = alpaka::createTaskKernel<Acc>(
      final_work_div,
      final_kernel,
      alpaka::getPtrNative(partial_acc),
      alpaka::getPtrNative(result_acc),
      static_cast<Idx>(partial_count));
    alpaka::enqueue(queue, final_task);

    alpaka::memcpy(queue, result_host, result_acc);
    alpaka::wait(queue);

    return result_value;
  }

} // namespace mudock::alpaka_demo

namespace alpaka::trait {
  template<typename TAcc>
  struct BlockSharedMemDynSizeBytes<mudock::alpaka_demo::dot_product_partial_reduction_kernel, TAcc> {
    template<typename TVec, typename... TArgs>
    ALPAKA_FN_HOST_ACC static auto getBlockSharedMemDynSizeBytes(
      mudock::alpaka_demo::dot_product_partial_reduction_kernel const&,
      TVec const& block_thread_extent,
      TVec const&,
      TArgs&&...) -> std::size_t {
      return static_cast<std::size_t>(block_thread_extent.prod()) * sizeof(double);
    }
  };

  template<typename TAcc>
  struct BlockSharedMemDynSizeBytes<mudock::alpaka_demo::final_reduction_kernel, TAcc> {
    template<typename TVec, typename... TArgs>
    ALPAKA_FN_HOST_ACC static auto getBlockSharedMemDynSizeBytes(
      mudock::alpaka_demo::final_reduction_kernel const&,
      TVec const& block_thread_extent,
      TVec const&,
      TArgs&&...) -> std::size_t {
      return static_cast<std::size_t>(block_thread_extent.prod()) * sizeof(double);
    }
  };
} // namespace alpaka::trait
