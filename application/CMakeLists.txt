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
add_executable("${CMAKE_PROJECT_NAME}" ${header_files} ${source_files})
target_include_directories("${CMAKE_PROJECT_NAME}" PRIVATE "${header_path}")
set_target_properties(
  "${CMAKE_PROJECT_NAME}"
  PROPERTIES CXX_STANDARD 20
             CXX_STANDARD_REQUIRED ON
             CXX_EXTENSIONS OFF
             POSITION_INDEPENDENT_CODE FALSE)
target_compile_options("${CMAKE_PROJECT_NAME}" PRIVATE -fno-plt)
target_link_options("${CMAKE_PROJECT_NAME}" PRIVATE -Wl,-z,now)

# add the dependencies
if(MUDOCK_ENABLE_SYCL)
  string(FIND "${ACPP_TARGETS}" "omp" SUBSTRING_INDEX)
  if(SUBSTRING_INDEX GREATER -1)
    target_link_libraries("${CMAKE_PROJECT_NAME}" omp)
  endif()
  target_link_libraries(
    "${CMAKE_PROJECT_NAME}"
    libmudock
    Boost::program_options
    Boost::context
    pthread
    Boost::fiber
    m
    stdc++)
else()
  target_link_libraries(
    "${CMAKE_PROJECT_NAME}"
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

if(MUDOCK_ENABLE_OMP)
  target_compile_options("${CMAKE_PROJECT_NAME}"
                         PRIVATE "-fopenmp-targets=nvptx64")
  target_link_options("${CMAKE_PROJECT_NAME}" PRIVATE "-fopenmp")
  target_link_options("${CMAKE_PROJECT_NAME}" PRIVATE
                      "-fopenmp-targets=nvptx64")
endif()

if(MUDOCK_ENABLE_LIKWID)
  target_include_directories("${CMAKE_PROJECT_NAME}"
                             PRIVATE ${LIKWID_INCLUDE_DIR})
  target_link_libraries("${CMAKE_PROJECT_NAME}" likwid::likwid)
  add_definitions(-DLIKWID_PERFMON)
endif()

if(MUDOCK_ENABLE_TEST)
  add_executable(test_grid
                 "${CMAKE_CURRENT_SOURCE_DIR}/src/tests/grid_test.cpp")
  target_include_directories(test_grid
                             PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src/")
  set_target_properties(
    test_grid
    PROPERTIES CXX_STANDARD 20
               CXX_STANDARD_REQUIRED ON
               CXX_EXTENSIONS OFF)
  target_link_libraries(test_grid PRIVATE libmudock Boost::program_options)

  add_executable(test_score
                 "${CMAKE_CURRENT_SOURCE_DIR}/src/tests/score_test.cpp")
  target_include_directories(test_score
                             PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src/")
  set_target_properties(
    test_score
    PROPERTIES CXX_STANDARD 20
               CXX_STANDARD_REQUIRED ON
               CXX_EXTENSIONS OFF)
  target_link_libraries(test_score PRIVATE libmudock Boost::program_options)
endif()

add_executable(converter "${CMAKE_CURRENT_SOURCE_DIR}/src/converter.cpp")
set_target_properties(
  converter
  PROPERTIES CXX_STANDARD 20
             CXX_STANDARD_REQUIRED ON
             CXX_EXTENSIONS OFF)
target_link_libraries(converter PRIVATE libmudock Boost::program_options)
target_compile_options(converter PRIVATE -fno-plt)
target_link_options(converter PRIVATE -Wl,-z,now)

# ##############################################################################
# ### Declare the install step
# ##############################################################################
add_executable(converter "${CMAKE_CURRENT_SOURCE_DIR}/src/converter.cpp")
set_target_properties(
  converter
  PROPERTIES CXX_STANDARD 20
             CXX_STANDARD_REQUIRED ON
             CXX_EXTENSIONS OFF)
target_link_libraries(converter PRIVATE libmudock Boost::program_options)
target_compile_options(converter PRIVATE -fno-plt)
target_link_options(converter PRIVATE -Wl,-z,now)

# ##############################################################################
# ### Declare the install step
# ##############################################################################

# define the install step
install(TARGETS "${CMAKE_PROJECT_NAME}"
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
