#include <cctype>
#include <mudock/compute/parse_ids.hpp>
#include <stdexcept>

namespace mudock {
  std::vector<std::size_t> parse_ids(std::string_view description) {
    std::vector<std::size_t> ids;
    while (!description.empty()) {
      // get the group description
      const auto end_group_description_index = description.find(',');
      const auto group_description           = description.substr(0, end_group_description_index);
      if (group_description.empty()) [[unlikely]] {
        throw std::runtime_error(std::string{"ID description has an empty group ("} +
                                 std::string{description} + std::string{")"});
      }

      // check if we are talking about a range
      const auto dash_index = group_description.find('-');
      if (dash_index != std::string_view::npos) {
        const auto first_number = std::stoull(std::string(group_description.substr(0, dash_index)));
        const auto second_number =
            std::stoull(std::string(group_description.substr(dash_index + std::size_t{1})));
        for (std::size_t i = first_number; i <= second_number; ++i) { ids.emplace_back(i); }
      } else { // or if we have a plain number
        ids.emplace_back(std::stoull(std::string(group_description)));
      }

      // remove the current group from the description
      description = description.substr(end_group_description_index);
    }
    return ids;
  }
} // namespace mudock
