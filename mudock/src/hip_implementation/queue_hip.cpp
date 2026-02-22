#include <mudock/compute/devices_memory.hpp>
#include <mudock/compute/reorder_buffer.hpp>
#include <mudock/hip_implementation/hip_utils.hpp>
#include <mudock/hip_implementation/queue_hip.hpp>
#include <mudock/log.hpp>
#include <mudock/type_alias.hpp>
#include <mudock/utils.hpp>
#include <mutex>

namespace mudock {
#ifdef MUDOCK_KERNEL_LOCK
  namespace {
    constexpr int k_max_devices_kernel_lock = 16;

    struct device_kernel_lock {
      std::mutex mutex;
      hipEvent_t event{};
      bool event_created{false};
      bool has_event{false};

      ~device_kernel_lock() noexcept(false) {
        if (event_created) {
          (void) hipEventDestroy(event); // best-effort at teardown; runtime may already be shutdown
        }
      }
    };

    device_memory_array<k_max_devices_kernel_lock, device_kernel_lock> hip_kernel_locks;

    device_kernel_lock* get_kernel_lock(const int dev) {
      hip_kernel_locks.init(dev, std::function<std::unique_ptr<device_kernel_lock>()>([&]() {
                              MUDOCK_CHECK(hipSetDevice(dev));
                              auto lock = std::make_unique<device_kernel_lock>();
                              MUDOCK_CHECK(hipEventCreateWithFlags(&lock->event, hipEventDisableTiming));
                              lock->event_created = true;
                              lock->has_event     = false;
                              return lock;
                            }));
      return hip_kernel_locks.v[dev].get_data();
    }
  } // namespace
#endif

  struct queue_hip::impl {
    hipStream_t stream;
    int device_id;
    impl(const int id) {
      device_id = id;
      MUDOCK_CHECK(hipSetDevice(static_cast<int>(id)));
      MUDOCK_CHECK(hipStreamCreate(&stream));
    };
    ~impl() noexcept(false) { MUDOCK_CHECK(hipStreamDestroy(stream)); };
  };

  queue_hip::queue_hip(const int _id, const device_type dev_type)
      : queue(_id, dev_type), impl_(std::make_unique<impl>(_id)) {
    assert(dev_type == device_type::GPU && "HIP supports only GPUs devices");
  };
  queue_hip::~queue_hip() = default; // unique_ptr will destroy Impl

  queue_hip::queue_hip(queue_hip&&) noexcept = default;

  //TODO find a better way to specific grid and block dimensions
  //TODO what about shared mem?
  void queue_hip::launch_kernel(void* f, void* args[], const index3D gridDim, const index3D blockDim) {
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
    if (lock->has_event) {
      MUDOCK_CHECK(hipStreamWaitEvent(impl_->stream, lock->event, 0));
    }
#endif

    MUDOCK_CHECK(hipLaunchKernel(f, grid, block, args, 0, impl_->stream));
    MUDOCK_CHECK_KERNELCALL();

#ifdef MUDOCK_KERNEL_LOCK
    MUDOCK_CHECK(hipEventRecord(lock->event, impl_->stream));
    lock->has_event = true;
#endif
  };
  void queue_hip::launch_kernel(void* f, void* args[], const int gridDim, const int blockDim) {
    assert(gridDim >= 0 && blockDim >= 0);
    launch_kernel(f, args, {gridDim, 1, 1}, {blockDim, 1, 1});
  };

  void queue_hip::alloc(void** ptr, const size_t bytes) {
    // TOOD check if required
    queue_hip::free(ptr);
    MUDOCK_CHECK(hipMallocAsync(ptr, bytes, impl_->stream));
  };
  void queue_hip::free(void** ptr) {
    if (*ptr != nullptr) {
      MUDOCK_CHECK(hipFreeAsync(*ptr, impl_->stream));
      *ptr = nullptr;
    }
  };
  void queue_hip::set_to_value(void* ptr, const size_t num_bytes, const char value) {
    MUDOCK_CHECK(hipMemsetAsync(ptr, value, num_bytes, impl_->stream));
  };
  void queue_hip::copy_host2device(const void* host, void* device, const size_t num_bytes) {
    MUDOCK_CHECK(hipMemcpyAsync(device, host, num_bytes, hipMemcpyHostToDevice, impl_->stream));
  };
  void queue_hip::copy_device2host(const void* device, void* host, const size_t num_bytes) {
    MUDOCK_CHECK(hipMemcpyAsync(host, device, num_bytes, hipMemcpyDeviceToHost, impl_->stream));
  };
  void queue_hip::copy_device2device(const void* device_src, void* device_dest, const size_t num_bytes) {
    MUDOCK_CHECK(hipMemcpyAsync(device_dest, device_src, num_bytes, hipMemcpyDeviceToDevice, impl_->stream));
  };
  void queue_hip::operator()() {
    // Setup default stream/device
  };

  void queue_hip::synchronize() { MUDOCK_CHECK(hipStreamSynchronize(impl_->stream)); };

} // namespace mudock
