list(APPEND sycl_compilers "oneAPI" "AdapativeCPP")

set(SYCL_COMPILER
    "AdaptiveCpp"
    CACHE STRING "SYCL compiler")
set_property(CACHE SYCL_COMPILER PROPERTY STRINGS ${sycl_compilers})

# ##############################################################################
# AdaptiveCPP Helper Functions #
# ##############################################################################

function(add_sycl_files SYCL_SOURCES HEADER_PATH HEADER_FILES OUTPUT)
  set(GENERATED_CPP_SOURCES "")
  message(STATUS "${HEADER_PATH}")
  string(TOUPPER "${CMAKE_BUILD_TYPE}" BUILD_TYPE_UPPER)
  get_target_property(MUDOCK_DEFINES libmudock COMPILE_DEFINITIONS)
  string(REPLACE ";" ";-D" MUDOCK_DEFINES "${MUDOCK_DEFINES}")
  set(MUDOCK_DEFINES "-D${MUDOCK_DEFINES}")
  set(COMPILE_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${BUILD_TYPE_UPPER}}")
  string(REPLACE " " ";" CXX_FLAGS_LIST "${COMPILE_FLAGS}")
  set(CXX_FLAGS_LIST "${CXX_FLAGS_LIST};${MUDOCK_DEFINES}")

  foreach(SYCL_FILE ${SYCL_SOURCES})
    get_filename_component(BASENAME ${SYCL_FILE} NAME_WE)
    if(SYCL_COMPILER STREQUAL "AdaptiveCpp")
      set(GENERATED_FILE "${CMAKE_CURRENT_BINARY_DIR}/${BASENAME}.acpp.o")

      add_custom_command(
        OUTPUT ${GENERATED_FILE}
        COMMAND
          ${ACPP_COMPILER} --acpp-targets=${SYCL_TARGETS} ${CXX_FLAGS_LIST}
          --std=c++20 -o ${GENERATED_FILE} -c ${SYCL_FILE} -I${HEADER_PATH}
        DEPENDS "${SYCL_FILE}" "${HEADER_FILES}"
        COMMENT "Compiling SYCL source ${SYCL_FILE} with acpp")
    elseif(SYCL_COMPILER STREQUAL "oneAPI")
      set(GENERATED_FILE "${CMAKE_CURRENT_BINARY_DIR}/${BASENAME}.oneapi.o")

      add_custom_command(
        OUTPUT ${GENERATED_FILE}
        COMMAND
          ${LLVM_TOOLS_BINARY_DIR}/clang++ -fsycl -fsycl-targets=${SYCL_TARGETS}
          ${CXX_FLAGS_LIST} --std=c++20 -o ${GENERATED_FILE} -c ${SYCL_FILE}
          ${BOOST_INCLUDE_DIRS} -I${HEADER_PATH} -I${Boost_INCLUDE_DIRS}
          -I${LLVM_INCLUDE_DIRS}
        DEPENDS "${SYCL_FILE}" "${HEADER_FILES}"
        COMMENT "Compiling SYCL source ${SYCL_FILE} with dpcpp")
    endif()

    list(APPEND GENERATED_CPP_SOURCES ${GENERATED_FILE})
  endforeach()
  add_custom_target(sycl_targets ALL DEPENDS "${GENERATED_CPP_SOURCES}")
  set(${OUTPUT}
      "${GENERATED_CPP_SOURCES}"
      PARENT_SCOPE)
endfunction()
