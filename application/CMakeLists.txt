# ##############################################################################
# ### Declare the application source files
# ##############################################################################

# headers
set(header_path "${CMAKE_CURRENT_SOURCE_DIR}/src")
list(APPEND header_files "${header_path}/command_line_args.hpp")

# sources
set(source_path "${CMAKE_CURRENT_SOURCE_DIR}/src")
list(APPEND source_files "${header_path}/command_line_args.cpp" "${source_path}/main.cpp")

# ##############################################################################
# ### Declare the compile step
# ##############################################################################

function(mudock_configure_application_target target include_dir)
  target_include_directories("${target}" PRIVATE "${include_dir}")
  set_target_properties(
    "${target}"
    PROPERTIES CXX_STANDARD 20
               CXX_STANDARD_REQUIRED ON
               CXX_EXTENSIONS OFF
               POSITION_INDEPENDENT_CODE FALSE)
  target_compile_options("${target}" PRIVATE -fno-plt)
  target_link_options("${target}" PRIVATE -Wl,-z,now)
  if(NOT MUDOCK_ENABLE_OMP)
    target_compile_options("${target}" PRIVATE -Wno-unknown-pragmas)
  endif()
  if(MUDOCK_ENABLE_GH)
    target_link_libraries("${target}" PRIVATE hwy::hwy)
  endif()
  if(MUDOCK_ENABLE_OMP OR MUDOCK_ENABLE_OMP_OFFLOAD)
    target_link_libraries("${target}" PRIVATE OpenMP::OpenMP_CXX)
  endif()
  if(MUDOCK_ENABLE_OMP_OFFLOAD)
    target_compile_options("${target}" PRIVATE "-fopenmp-targets=nvptx64")
    target_link_options("${target}" PRIVATE "-fopenmp-targets=nvptx64")
  endif()
endfunction()

function(mudock_link_application_target target)
  target_link_libraries("${target}" PRIVATE libmudock Boost::program_options)
endfunction()

# declare the application executable
add_executable_mudock("${CMAKE_PROJECT_NAME}" "${header_files}"
                      "${source_files}")
mudock_configure_application_target("${CMAKE_PROJECT_NAME}" "${header_path}")
mudock_link_application_target("${CMAKE_PROJECT_NAME}")
if(MUDOCK_ENABLE_MPI)
  target_link_libraries("${CMAKE_PROJECT_NAME}" PRIVATE MPI::MPI_CXX)
endif()

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/src/bench" bench)
if(MUDOCK_ENABLE_TEST)
  add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/src/tests" tests)
  add_executable_mudock(test_grid
                        "${CMAKE_CURRENT_SOURCE_DIR}/src/tests/grid_test.cpp")
  target_include_directories(test_grid
                             PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src/")
  set_target_properties(
    test_grid
    PROPERTIES CXX_STANDARD 20
               CXX_STANDARD_REQUIRED ON
               CXX_EXTENSIONS OFF)
  target_link_libraries(test_grid libmudock Boost::program_options)

  add_executable_mudock(
    precomputed_test_score_adt "${CMAKE_CURRENT_SOURCE_DIR}/src/tests/precomputed_adt_score_test.cpp")
  target_include_directories(precomputed_test_score_adt
                             PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src/")
  set_target_properties(
    precomputed_test_score_adt
    PROPERTIES CXX_STANDARD 20
               CXX_STANDARD_REQUIRED ON
               CXX_EXTENSIONS OFF)
  target_link_libraries(precomputed_test_score_adt libmudock Boost::program_options)


  add_executable_mudock(
    test_score_adt "${CMAKE_CURRENT_SOURCE_DIR}/src/tests/adt_score_test.cpp")
  target_include_directories(test_score_adt
                             PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src/")
  set_target_properties(
    test_score_adt
    PROPERTIES CXX_STANDARD 20
               CXX_STANDARD_REQUIRED ON
               CXX_EXTENSIONS OFF)
  target_link_libraries(test_score_adt libmudock Boost::program_options)

  add_executable_mudock(test_score_adt_impl
                        "${CMAKE_CURRENT_SOURCE_DIR}/src/tests/adt_impl_score_test.cpp")
  target_include_directories(test_score_adt_impl
                             PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src/")
  set_target_properties(
    test_score_adt_impl
    PROPERTIES CXX_STANDARD 20
               CXX_STANDARD_REQUIRED ON
               CXX_EXTENSIONS OFF)
  target_link_libraries(test_score_adt_impl libmudock Boost::program_options)
endif()

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
mudock_configure_application_target(converter "${header_path}")
mudock_link_application_target(converter)

# ##############################################################################
# ### Declare the install step
# ##############################################################################

# define the install step
install(TARGETS "${CMAKE_PROJECT_NAME}"
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
