#include <mudock/log.hpp>
#include <omp.h>

namespace mudock {
  void set_device(const std::size_t device_id) {
    // TODO polygeist missing translation
    // TODO check linkinf with external dependencies using
    // https://github.com/ivanradanov/cpucuda_runtime
    omp_set_default_device(device_id);
    info("Worker CUDA on duty! Set affinity to device ", device_id);
  }
} // namespace mudock
