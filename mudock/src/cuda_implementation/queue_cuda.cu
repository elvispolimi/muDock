#include <cstdlib>
#include <memory>
#include <mudock/compute/devices_memory.hpp>
#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/cuda_implementation/cuda_utils.cuh>
#include <mudock/cuda_implementation/queue_cuda.cuh>
#include <mudock/log.hpp>
#include <mudock/molecule.hpp>
#include <mudock/type_alias.hpp>
#include <mudock/utils.hpp>
#include <mutex>

namespace mudock {
#ifdef MUDOCK_KERNEL_LOCK
  namespace {
    constexpr int k_max_devices = 16;

    struct device_kernel_lock {
      std::mutex mutex;
      cudaEvent_t event{};
      bool event_created{false};
      bool has_previous_event{false};

      ~device_kernel_lock() noexcept(false) {
        if (event_created) {
          cudaEventDestroy(event);
        }
      }
    };

    device_memory_array<k_max_devices, device_kernel_lock> cuda_kernel_locks;

    device_kernel_lock* get_kernel_lock(const int dev) {
      cuda_kernel_locks.init(dev, std::function<std::unique_ptr<device_kernel_lock>()>([&]() {
                               MUDOCK_CHECK(cudaSetDevice(dev));
                               auto lock = std::make_unique<device_kernel_lock>();
                               MUDOCK_CHECK(cudaEventCreateWithFlags(&lock->event, cudaEventDisableTiming));
                               lock->event_created      = true;
                               lock->has_previous_event = false;
                               return lock;
                             }));
      return cuda_kernel_locks.v[dev].get_data();
    }
  } // namespace
#endif

  struct queue_cuda::impl {
    cudaStream_t stream;
    int device_id;
    impl(const int id) {
      device_id = id;
      MUDOCK_CHECK(cudaSetDevice(static_cast<int>(id)));
      MUDOCK_CHECK(cudaStreamCreate(&stream));
    };
    ~impl() noexcept(false) { MUDOCK_CHECK(cudaStreamDestroy(stream)); };
  };

  queue_cuda::queue_cuda(const int _id, const device_type _dev_type)
      : queue(_id, _dev_type), impl_(std::make_unique<impl>(_id)) {
    assert(dev_type == device_type::GPU && "CUDA supports only GPUs devices");
  };
  queue_cuda::~queue_cuda() = default; // unique_ptr will destroy Impl

  queue_cuda::queue_cuda(queue_cuda&&) noexcept = default;

  //TODO find a better way to specific grid and block dimensions
  //TODO what about shared mem?
  void queue_cuda::launch_kernel(void* f, void* args[], const index3D gridDim, const index3D blockDim) {
    assert(gridDim.size_x() > 0 && blockDim.size_x() > 0);
    assert(gridDim.size_y() > 0 && blockDim.size_y() > 0);
    assert(gridDim.size_z() > 0 && blockDim.size_z() > 0);

    const dim3 grid  = dim3{static_cast<unsigned int>(gridDim.size_x()),
                            static_cast<unsigned int>(gridDim.size_y()),
                            static_cast<unsigned int>(gridDim.size_z())};
    const dim3 block = dim3{static_cast<unsigned int>(blockDim.size_x()),
                            static_cast<unsigned int>(blockDim.size_y()),
                            static_cast<unsigned int>(blockDim.size_z())};

    // const std::size_t shared_mem =
    //     std::max(configuration.population_number, static_cast<std::size_t>(BLOCK_SIZE)) * sizeof(fp_type);
#ifdef MUDOCK_KERNEL_LOCK
    auto* lock = get_kernel_lock(impl_->device_id);
    std::unique_lock<std::mutex> guard(lock->mutex);
    if (lock->has_previous_event) {
      // Cross-stream dependency on the previous kernel recorded for this device.
      MUDOCK_CHECK(cudaStreamWaitEvent(impl_->stream, lock->event, 0));
    }
#endif

    MUDOCK_CHECK(cudaLaunchKernel(f, grid, block, args, 0, impl_->stream));
    MUDOCK_CHECK_KERNELCALL();

#ifdef MUDOCK_KERNEL_LOCK
    MUDOCK_CHECK(cudaEventRecord(lock->event, impl_->stream));
    lock->has_previous_event = true;
#endif
  };
  void queue_cuda::launch_kernel(void* f, void* args[], const int gridDim, const int blockDim) {
    assert(gridDim >= 0 && blockDim >= 0);
    launch_kernel(f, args, {gridDim, 1, 1}, {blockDim, 1, 1});
  };

  void queue_cuda::alloc(void** ptr, const size_t bytes) {
    // TOOD check if required
    queue_cuda::free(ptr);
    MUDOCK_CHECK(cudaMallocAsync(ptr, bytes, impl_->stream));
  };
  void queue_cuda::free(void** ptr) {
    if (*ptr != nullptr) {
      MUDOCK_CHECK(cudaFreeAsync(*ptr, impl_->stream));
      *ptr = nullptr;
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
