add_executable_mudock(adt_score_bench
                      "${CMAKE_CURRENT_SOURCE_DIR}/adt_score_bench.cpp"
                      "${CMAKE_CURRENT_SOURCE_DIR}/../command_line_args.cpp")
target_include_directories(adt_score_bench
                           PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/..")
set_target_properties(
  adt_score_bench
  PROPERTIES CXX_STANDARD 20
             CXX_STANDARD_REQUIRED ON
             CXX_EXTENSIONS OFF)
if(MUDOCK_ENABLE_SYCL)
  string(FIND "${ACPP_TARGETS}" "omp" SUBSTRING_INDEX)
  if(SUBSTRING_INDEX GREATER -1)
    target_link_libraries(adt_score_bench omp)
  endif()
  target_link_libraries(
    adt_score_bench
    libmudock
    Boost::program_options
    Boost::context
    pthread
    Boost::fiber
    m
    stdc++)
else()
  target_link_libraries(
    adt_score_bench
    libmudock
    Boost::program_options
    Boost::context
    pthread
    m
    stdc++)
endif()
if(MUDOCK_ENABLE_GH)
  target_link_libraries(adt_score_bench hwy::hwy)
endif()
if(MUDOCK_ENABLE_OMP OR MUDOCK_ENABLE_OMP_OFFLOAD)
  target_link_libraries(adt_score_bench OpenMP::OpenMP_CXX)
endif()
if(MUDOCK_ENABLE_OMP_OFFLOAD)
  target_compile_options(adt_score_bench PRIVATE "-fopenmp-targets=nvptx64")
  target_link_options(adt_score_bench PRIVATE "-fopenmp-targets=nvptx64")
endif()
