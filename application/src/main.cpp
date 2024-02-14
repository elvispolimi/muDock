#include "command_line_args.hpp"

#include <iostream>
#include <mudock/mudock.hpp>

int main(int argc, char* argv[]) {
  const auto args = parse_command_line_arguments(argc, argv);

  return EXIT_SUCCESS;
}
