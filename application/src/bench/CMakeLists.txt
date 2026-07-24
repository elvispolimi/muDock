add_executable_mudock(adt_score_bench
                      "${CMAKE_CURRENT_SOURCE_DIR}/adt_score_bench.cpp"
                      "${CMAKE_CURRENT_SOURCE_DIR}/../command_line_args.cpp")
mudock_configure_application_target(adt_score_bench "${CMAKE_CURRENT_SOURCE_DIR}/..")
mudock_link_application_target(adt_score_bench)
