#include <algorithm>
#include <stdexcept>
#include <mudock/format/format_id.hpp>

namespace mudock {
  supported_format parse_supported_format(const std::string_view extension) {
    const auto element_it = std::find_if(std::begin(FORMAT_EXTENSIONS),
                                         std::end(FORMAT_EXTENSIONS),
                                         [&extension](const auto& e) { return e.extension == extension; });
    if (element_it != std::end(FORMAT_EXTENSIONS))
      return element_it->format;
    else
      throw std::runtime_error("Missing extension");
  }

  std::string_view parse_supported_format(const supported_format format) {
    const auto element_it = std::find_if(std::begin(FORMAT_EXTENSIONS),
                                         std::end(FORMAT_EXTENSIONS),
                                         [&format](const auto& e) { return e.format == format; });
    if (element_it != std::end(FORMAT_EXTENSIONS))
      return element_it->extension;
    else
      throw std::runtime_error("Missing format");
  }
} // namespace mudock
