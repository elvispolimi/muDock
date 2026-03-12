# headers
set(header_path "${CMAKE_CURRENT_SOURCE_DIR}/..")
list(APPEND mpi_header_files "${header_path}/command_line_args.hpp")

# sources
set(source_path "${CMAKE_CURRENT_SOURCE_DIR}")
list(APPEND mpi_source_files "${header_path}/command_line_args.cpp"
    "${source_path}/main_mpi.cpp")

# ##############################################################################
# ### Declare the compile step
# ##############################################################################

# declare the MPI executable
set(files ${mpi_header_files} ${mpi_source_files})
add_executable_mudock("${CMAKE_PROJECT_NAME}_mpi" "${mpi_header_files}"
                    "${mpi_source_files}")
target_include_directories("${CMAKE_PROJECT_NAME}_mpi" PRIVATE "${header_path}")
set_target_properties(
"${CMAKE_PROJECT_NAME}_mpi"
PROPERTIES CXX_STANDARD 20
            CXX_STANDARD_REQUIRED ON
            CXX_EXTENSIONS OFF
            POSITION_INDEPENDENT_CODE FALSE)

# Same flags as main executable
target_compile_options("${CMAKE_PROJECT_NAME}_mpi" PRIVATE -fno-plt)
target_link_options("${CMAKE_PROJECT_NAME}_mpi" PRIVATE -Wl,-z,now)
if(NOT MUDOCK_ENABLE_OMP)
target_compile_options("${CMAKE_PROJECT_NAME}_mpi" PRIVATE -Wno-unknown-pragmas)
endif()

# add the dependencies
if(MUDOCK_ENABLE_SYCL)
string(FIND "${ACPP_TARGETS}" "omp" SUBSTRING_INDEX)
if(SUBSTRING_INDEX GREATER -1)
    target_link_libraries("${CMAKE_PROJECT_NAME}_mpi" omp)
endif()
target_link_libraries(
    "${CMAKE_PROJECT_NAME}_mpi"
    libmudock
    Boost::program_options
    Boost::context
    pthread
    Boost::fiber
    m
    stdc++
    MPI::MPI_CXX)
else()
target_link_libraries(
    "${CMAKE_PROJECT_NAME}_mpi"
    libmudock
    Boost::program_options
    Boost::context
    pthread
    m
    stdc++
    MPI::MPI_CXX)
endif()

if(MUDOCK_ENABLE_GH)
target_link_libraries("${CMAKE_PROJECT_NAME}_mpi" hwy::hwy)
endif()

if(MUDOCK_ENABLE_OMP OR MUDOCK_ENABLE_OMP_OFFLOAD)
target_link_libraries("${CMAKE_PROJECT_NAME}_mpi" OpenMP::OpenMP_CXX)
endif()

if(MUDOCK_ENABLE_OMP_OFFLOAD)
target_compile_options("${CMAKE_PROJECT_NAME}_mpi"
                        PRIVATE "-fopenmp-targets=nvptx64")
target_link_options("${CMAKE_PROJECT_NAME}_mpi"
                    PRIVATE "-fopenmp-targets=nvptx64")
endif()
