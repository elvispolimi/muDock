# headers
set(header_path "${CMAKE_CURRENT_SOURCE_DIR}/..")
list(APPEND mpi_header_files "${header_path}/command_line_args.hpp")

# sources
set(source_path "${CMAKE_CURRENT_SOURCE_DIR}")
list(APPEND mpi_source_files "${header_path}/command_line_args.cpp"
    "${header_path}/tbb/tbb_stream_entry.cpp"
    "${source_path}/main_mpi.cpp")

# ##############################################################################
# ### Declare the compile step
# ##############################################################################

# declare the MPI executable
add_executable_mudock("${CMAKE_PROJECT_NAME}_mpi" "${mpi_header_files}"
                    "${mpi_source_files}")
mudock_configure_application_target("${CMAKE_PROJECT_NAME}_mpi" "${header_path}")

# add the dependencies
mudock_link_application_target("${CMAKE_PROJECT_NAME}_mpi")
target_link_libraries("${CMAKE_PROJECT_NAME}_mpi" PRIVATE MPI::MPI_CXX)
