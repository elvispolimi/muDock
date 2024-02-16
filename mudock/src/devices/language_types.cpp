#include <array>
#include <mudock/devices/language_types.hpp>
#include <stdexcept>

namespace mudock {

  static std::array<std::string_view, 1> language_kind_names = {{
      "CPP",
  }};

  language_types parse_language(std::string language_kind) {
    for (std::size_t i{0}; i < language_kind_names.size(); ++i) {
      if (language_kind == language_kind_names[i]) {
        return static_cast<language_types>(i);
      }
    }
    throw std::runtime_error{"Invalid language  " + language_kind + " !"};
  }
} // namespace mudock
