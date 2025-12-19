#include <cstdlib>
#include <memory>
#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/cuda_implementation/cuda_utils.cuh>
#include <mudock/cuda_implementation/queue_cuda.cuh>
#include <mudock/log.hpp>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>
#include <mudock/utils.hpp>

namespace mudock {
  struct queue_cuda::impl {
    cudaStream_t stream;
    impl(const int id) {
      MUDOCK_CHECK(cudaSetDevice(static_cast<int>(id)));
      MUDOCK_CHECK(cudaStreamCreate(&stream));
    };
    ~impl() noexcept(false) { MUDOCK_CHECK(cudaStreamDestroy(stream)); };
  };

  queue_cuda::queue_cuda(const int _id): queue(_id), impl_(std::make_unique<impl>(_id)) {};
  queue_cuda::~queue_cuda() = default; // unique_ptr will destroy Impl

  queue_cuda::queue_cuda(queue_cuda&&) noexcept            = default;
  queue_cuda& queue_cuda::operator=(queue_cuda&&) noexcept = default;

  void queue_cuda::launch_kernel(void* f, const int batch_ligands, void* args[]) {
    // const std::size_t shared_mem =
    //     std::max(configuration.population_number, static_cast<std::size_t>(BLOCK_SIZE)) * sizeof(fp_type);
    MUDOCK_CHECK(cudaLaunchKernel(f, batch_ligands, BLOCK_SIZE, args, 0, impl_->stream));
    MUDOCK_CHECK_KERNELCALL();
  };

  void queue_cuda::alloc(void** ptr, const size_t bytes) {
    // TOOD check if required
    if (*ptr != nullptr)
      free(ptr);
    MUDOCK_CHECK(cudaMallocAsync(ptr, bytes, impl_->stream));
  };
  void queue_cuda::free(void** ptr) {
    if (*ptr != nullptr) {
      MUDOCK_CHECK(cudaFreeAsync(*ptr, impl_->stream));
      ptr = nullptr;
    }
  };
  void queue_cuda::set_to_value(void* ptr, const size_t num_bytes, const char value) {
    MUDOCK_CHECK(cudaMemsetAsync(ptr, value, num_bytes, impl_->stream));
  };
  void queue_cuda::copy_host2device(const void* host, void* device, const size_t num_bytes) {
    MUDOCK_CHECK(cudaMemcpyAsync(device, host, num_bytes, cudaMemcpyHostToDevice, impl_->stream));
  };
  void queue_cuda::copy_device2host(const void* device, void* host, const size_t num_bytes) {
    MUDOCK_CHECK(cudaMemcpyAsync(host, device, num_bytes, cudaMemcpyDeviceToHost, impl_->stream));
  };
  void queue_cuda::copy_device2device(const void* device_src, void* device_dest, const size_t num_bytes) {
    MUDOCK_CHECK(
        cudaMemcpyAsync(device_dest, device_src, num_bytes, cudaMemcpyDeviceToDevice, impl_->stream));
  };
  void queue_cuda::operator()() {
    // Setup default stream/device
  };

  void queue_cuda::synchronize() { MUDOCK_CHECK(cudaStreamSynchronize(impl_->stream)); };

} // namespace mudock
