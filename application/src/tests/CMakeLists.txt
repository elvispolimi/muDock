add_executable_mudock(test_grid
                      "${CMAKE_CURRENT_SOURCE_DIR}/grid_test.cpp")
target_include_directories(test_grid
                            PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/..")
target_compile_definitions(
  test_grid PRIVATE MUDOCK_CTEST_SKIP_RETURN_CODE=${MUDOCK_CTEST_SKIP_RETURN_CODE})
set_target_properties(
  test_grid
  PROPERTIES CXX_STANDARD 20
              CXX_STANDARD_REQUIRED ON
              CXX_EXTENSIONS OFF)
target_link_libraries(test_grid libmudock Boost::program_options)

  add_executable_mudock(
    test_score_adt "${CMAKE_CURRENT_SOURCE_DIR}/adt_score_test.cpp")
  target_include_directories(test_score_adt
                              PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/..")
  set_target_properties(
    test_score_adt
    PROPERTIES CXX_STANDARD 20
                CXX_STANDARD_REQUIRED ON
                CXX_EXTENSIONS OFF)
  target_link_libraries(test_score_adt libmudock Boost::program_options)

  add_executable_mudock(test_score_adt_impl
                        "${CMAKE_CURRENT_SOURCE_DIR}/adt_impl_score_test.cpp")
  target_include_directories(test_score_adt_impl
                              PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/..")
  set_target_properties(
    test_score_adt_impl
    PROPERTIES CXX_STANDARD 20
                CXX_STANDARD_REQUIRED ON
                CXX_EXTENSIONS OFF)
  target_link_libraries(test_score_adt_impl libmudock Boost::program_options)

  add_executable_mudock(test_score_x "${CMAKE_CURRENT_SOURCE_DIR}/x_score_test.cpp")
  target_include_directories(test_score_x
                              PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/..")
  target_compile_definitions(
    test_score_x PRIVATE MUDOCK_CTEST_SKIP_RETURN_CODE=${MUDOCK_CTEST_SKIP_RETURN_CODE})
  set_target_properties(
    test_score_x
    PROPERTIES CXX_STANDARD 20
                CXX_STANDARD_REQUIRED ON
                CXX_EXTENSIONS OFF)
  target_link_libraries(test_score_x libmudock Boost::program_options)