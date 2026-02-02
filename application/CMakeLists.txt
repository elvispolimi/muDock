# ##############################################################################
# ### Declare the application source files
# ##############################################################################

# headers
set(header_path "${CMAKE_CURRENT_SOURCE_DIR}/src")
list(APPEND header_files "${header_path}/command_line_args.hpp")

# sources
set(source_path "${CMAKE_CURRENT_SOURCE_DIR}/src")
list(APPEND source_files "${header_path}/command_line_args.cpp"
     "${source_path}/main.cpp")

# ##############################################################################
# ### Declare the compile step
# ##############################################################################

# declare the application executable
set(files ${header_files} ${source_files})
add_executable_mudock("${CMAKE_PROJECT_NAME}" "${header_files}"
                      "${source_files}")
target_include_directories("${CMAKE_PROJECT_NAME}" PRIVATE "${header_path}")
set_target_properties(
  "${CMAKE_PROJECT_NAME}"
  PROPERTIES CXX_STANDARD 20
             CXX_STANDARD_REQUIRED ON
             CXX_EXTENSIONS OFF
             POSITION_INDEPENDENT_CODE FALSE)
# Required by OpenBabel dynamic plugin/libraries loading
target_compile_options("${CMAKE_PROJECT_NAME}" PRIVATE -fno-plt)
target_link_options("${CMAKE_PROJECT_NAME}" PRIVATE -Wl,-z,now)
if(NOT MUDOCK_ENABLE_OMP)
  target_compile_options("${CMAKE_PROJECT_NAME}" PRIVATE -Wno-unknown-pragmas)
endif()

# add the dependencies
if(MUDOCK_ENABLE_SYCL)
  string(FIND "${ACPP_TARGETS}" "omp" SUBSTRING_INDEX)
  if(SUBSTRING_INDEX GREATER -1)
    target_link_libraries("${CMAKE_PROJECT_NAME}" omp)
  endif()
  target_link_libraries(
    ${CMAKE_PROJECT_NAME}
    libmudock
    Boost::program_options
    Boost::context
    pthread
    Boost::fiber
    m
    stdc++)
else()
  target_link_libraries(
    ${CMAKE_PROJECT_NAME}
    libmudock
    Boost::program_options
    Boost::context
    pthread
    m
    stdc++)
endif()

if(MUDOCK_ENABLE_GH)
  target_link_libraries("${CMAKE_PROJECT_NAME}" hwy::hwy)
endif()

if(MUDOCK_ENABLE_OMP OR MUDOCK_ENABLE_OMP_OFFLOAD)
  target_link_libraries("${CMAKE_PROJECT_NAME}" OpenMP::OpenMP_CXX)
endif()

if(MUDOCK_ENABLE_OMP_OFFLOAD)
  target_compile_options("${CMAKE_PROJECT_NAME}"
                         PRIVATE "-fopenmp-targets=nvptx64")
  # target_link_options("${CMAKE_PROJECT_NAME}" PRIVATE "-fopenmp")
  target_link_options("${CMAKE_PROJECT_NAME}" PRIVATE
                      "-fopenmp-targets=nvptx64")
endif()

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/src/bench" bench)
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/src/tests" tests)

# if(MUDOCK_ENABLE_LIKWID)
#   # target_include_directories("${CMAKE_PROJECT_NAME}" PRIVATE
#   # ${LIKWID_INCLUDE_DIR})
#   target_link_libraries("${CMAKE_PROJECT_NAME}" likwid::likwid)
#   target_link_libraries(test_score likwid::likwid)
#
#   add_definitions(-DLIKWID_PERFMON)
# endif()

# ##############################################################################
# ### Declare the install step
# ##############################################################################
add_executable_mudock(converter "${CMAKE_CURRENT_SOURCE_DIR}/src/converter.cpp")
set_target_properties(
  converter
  PROPERTIES CXX_STANDARD 20
             CXX_STANDARD_REQUIRED ON
             CXX_EXTENSIONS OFF)
target_link_libraries(converter libmudock Boost::program_options)
target_compile_options(converter PRIVATE -fno-plt)
target_link_options(converter PRIVATE -Wl,-z,now)

# ##############################################################################
# ### Declare the install step
# ##############################################################################

# define the install step
install(TARGETS "${CMAKE_PROJECT_NAME}"
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
