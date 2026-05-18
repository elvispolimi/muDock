add_executable_mudock(local_search
                      "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp"
                      "${CMAKE_CURRENT_SOURCE_DIR}/../command_line_args.cpp")
mudock_configure_application_target(local_search "${CMAKE_CURRENT_SOURCE_DIR}/..")
mudock_link_application_target(local_search)
