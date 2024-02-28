#pragma once

#include <string>

namespace mudock {
  // Enumerates the languages that are support by muDock
  // NOTE:  NONE is only used to get the end fo the Enum
  //        it is actually a trick to enable iterating over the enum
  enum class language_types : int { CPP = 0, NONE };

  // From a string you get the language type, if available
  // Otherwise it throws a runtime exception 
  language_types parse_language(std::string);

} // namespace mudock
