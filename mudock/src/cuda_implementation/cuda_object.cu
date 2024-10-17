#include <cuda_runtime.h>
#include <mudock/cpp_implementation/chromosome.hpp>
#include <mudock/cuda_implementation/cuda_check_error_macro.cuh>
#include <mudock/cuda_implementation/cuda_object.cuh>
#ifdef MUDOCK_ENABLE_POLY
  #include <polygeist/cuda_random.cuh>
#else
  #include <mudock/cuda_implementation/cuda_random.cuh>
#endif
#include <mudock/type_alias.hpp>

namespace mudock {
  template<class T>
  cuda_object<T>::cuda_object(cuda_object&& other): stream(other.stream) {
    dev_ptr       = other.dev_ptr;
    size          = other.size;
    other.dev_ptr = nullptr;
    other.size    = 0;
  }

// TODO destructor not supported by Polygeist
#ifndef MUDOCK_ENABLE_POLY
  template<class T>
  cuda_object<T>::~cuda_object() noexcept(false) {
    if (dev_ptr != nullptr)
      MUDOCK_CHECK(cudaFreeAsync(dev_ptr, stream));
  }
#endif

  template<class T>
  void cuda_object<T>::alloc(const size_t num_elements) {
    if (size < num_elements) {
      if (dev_ptr != nullptr)
        MUDOCK_CHECK(cudaFreeAsync(dev_ptr, stream));
      MUDOCK_CHECK(cudaMallocAsync(&dev_ptr, sizeof(T) * num_elements, stream));
    }
    size = num_elements;
  }
  template<class T>
  void cuda_object<T>::set_to_value(const int value) {
    MUDOCK_CHECK(cudaMemsetAsync(dev_ptr, value, sizeof(T) * size, stream));
  }

  template<class T>
  void cuda_object<T>::copy_host2device(const T* const host) {
    MUDOCK_CHECK(cudaMemcpyAsync(dev_ptr, host, sizeof(T) * size, cudaMemcpyHostToDevice, stream));
  }
  template<class T>
  void cuda_object<T>::copy_device2host(T* const host) const {
    MUDOCK_CHECK(cudaMemcpyAsync(host, dev_ptr, sizeof(T) * size, cudaMemcpyDeviceToHost, stream));
  }

  template<class T>
  [[nodiscard]] T* cuda_object<T>::dev_pointer() const {
    return dev_ptr;
  }
  template<class T>
  [[nodiscard]] std::size_t cuda_object<T>::num_elements() const {
    return size;
  }
  template<class T>
  [[nodiscard]] const cudaStream_t& cuda_object<T>::get_stream() const {
    return stream;
  }

  template class cuda_object<int>;
  template class cuda_object<fp_type>;
  template class cuda_object<fp_type*>;
#ifdef MUDOCK_ENABLE_POLY
  template class cuda_object<XORWOWState>;
#else
  template class cuda_object<curandState>;
  template class cuda_object<cudaTextureObject_t>;
#endif
  template class cuda_object<chromosome>;
} // namespace mudock
