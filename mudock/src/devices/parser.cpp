#include <cstdint>
#include <exception>
#include <mudock/devices/parser.hpp>
#include <sstream>
#include <stdexcept>
#include <string>

namespace mudock {
  enum class conf_fsm { language, device, id, finish };

  device_conf parse_conf(const std::string conf) {
    device_types d_t{device_types::CPU};
    language_types l_t{language_types::CPP};
    std::vector<std::size_t> ids;

    std::istringstream line_tokens{conf};
    std::string token;
    conf_fsm state = conf_fsm::language;
    while (std::getline(line_tokens, token, ':') && state != conf_fsm::finish) {
      switch (state) {
        case conf_fsm::language: {
          l_t   = parse_language(token);
          state = conf_fsm::device;
          break;
        }
        case conf_fsm::device: {
          d_t   = parse_device(token);
          state = conf_fsm::id;
          break;
        }
        case conf_fsm::id: {
          std::istringstream id_tokens{token};
          std::string id_token;
          while (std::getline(id_tokens, id_token, ',')) {
            if (id_token.find('-') != std::string::npos) {
              std::istringstream range(token);
              int start, end;
              char dash;
              if (range >> start >> dash >> end) {
                for (int i = start; i <= end; ++i) { ids.push_back(i); }
              }
            } else {
              ids.push_back(std::stoi(id_token));
            }
          }
          state = conf_fsm::finish;
          break;
        }
        default: throw std::runtime_error{"Error configuration parse state!"};
      };
    }

    if (state != conf_fsm::finish)
      throw std::runtime_error{"Cannot parse input devices configuration:" + conf + "!"};

    return {d_t, l_t, ids};
  };
} // namespace mudock
