#pragma once

#include <string>

namespace mudock {
  // Enumerates the kernels type that are support by muDock
  // NOTE:  NONE is only used to get the end fo the Enum
  //        it is actually a trick to enable iterating over the enum
  enum class kernel_type : int { CPP = 0, NONE };

  // From a string you get the kernel type, if available
  // Otherwise it throws a runtime exception 
  kernel_type parse_kernel(const std::string& kernel_kind);

} // namespace mudock
