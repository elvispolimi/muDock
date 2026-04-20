#include "command_line_args.hpp"
#include "tbb_stream_entry.hpp"

int main(int argc, char** argv) {
  const auto args = parse_command_line_arguments(argc, argv);
  return run_tbb_stream_entry(args);
}
