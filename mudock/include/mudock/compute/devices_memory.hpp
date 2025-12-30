#pragma once

#include <functional>
#include <memory>
#include <mutex>

namespace mudock {
  template<typename T>
  struct device_memory {
    std::unique_ptr<T> data;
    std::once_flag once;

    inline void init(std::function<std::unique_ptr<T>()> f) {
      std::call_once(once, [&] { data = f(); });
    };
    template<typename... Args>
    void init(Args... args) {
      std::call_once(once, [&] { data = std::make_unique<T>(args...); });
    }
    T* get_data() { return data.get(); }

    T& operator()() { return *data; }
  };

  template<int max_num_devices, typename T>
  struct device_memory_array {
    std::array<device_memory<T>, max_num_devices> v;

    inline void init(const int id, std::function<std::unique_ptr<T>()> f) { v[id].init(f); };
    template<typename... Args>
    void init(const int id, Args... args) {
      v[id].init(id, args...);
    }
  };

} // namespace mudock
