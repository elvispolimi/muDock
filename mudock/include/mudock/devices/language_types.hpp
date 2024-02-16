#pragma once

#include <string>

namespace mudock {
  enum class language_types : int { CPP = 0, NONE };

  language_types parse_language(std::string);

} // namespace mudock
