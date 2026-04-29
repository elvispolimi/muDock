# ##############################################################################
# oneAPI SYCL Helper Functions #
# ##############################################################################

function(add_sycl_files SYCL_SOURCES HEADER_PATH HEADER_FILES OUTPUT)
  set(GENERATED_CPP_SOURCES "")
  set(SYCL_SYSTEM_INCLUDE_FLAGS "")
  string(TOUPPER "${CMAKE_BUILD_TYPE}" BUILD_TYPE_UPPER)
  get_target_property(MUDOCK_DEFINES libmudock COMPILE_DEFINITIONS)
  string(REPLACE ";" ";-D" MUDOCK_DEFINES "${MUDOCK_DEFINES}")
  set(MUDOCK_DEFINES "-D${MUDOCK_DEFINES}")
  set(COMPILE_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${BUILD_TYPE_UPPER}}")
  string(REPLACE " " ";" CXX_FLAGS_LIST "${COMPILE_FLAGS}")
  set(CXX_FLAGS_LIST "${CXX_FLAGS_LIST};${MUDOCK_DEFINES}")

  foreach(INCLUDE_DIR IN LISTS BOOST_INCLUDE_DIRS Boost_INCLUDE_DIRS
                                 LLVM_INCLUDE_DIRS)
    if(INCLUDE_DIR)
      string(REGEX REPLACE "^-I" "" INCLUDE_DIR "${INCLUDE_DIR}")
      list(APPEND SYCL_SYSTEM_INCLUDE_FLAGS "-isystem" "${INCLUDE_DIR}")
    endif()
  endforeach()

  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(SYCL_EXTRA_FLAGS "-O0" "-g")
  elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    set(SYCL_EXTRA_FLAGS "-O2" "-g")
  elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(SYCL_EXTRA_FLAGS "-O3" "-DNDEBUG")
  elseif(CMAKE_BUILD_TYPE STREQUAL "MinSizeRel")
    set(SYCL_EXTRA_FLAGS "-Os" "-DNDEBUG")
  endif()

  foreach(SYCL_FILE ${SYCL_SOURCES})
    get_filename_component(BASENAME ${SYCL_FILE} NAME_WE)
    set(GENERATED_FILE "${CMAKE_CURRENT_BINARY_DIR}/${BASENAME}.oneapi.o")
  
    #TOOD fix the -Wno-sign-conversion
    add_custom_command(
      OUTPUT ${GENERATED_FILE}
      COMMAND
        ${LLVM_TOOLS_BINARY_DIR}/clang++ -fsycl -fsycl-targets=${SYCL_TARGETS}
        ${SYCL_BACKEND_FLAGS_COMPILE} ${CXX_FLAGS_LIST} ${global_c_cxx_flags} -Wno-sign-conversion --std=c++20 -MMD -MF ${GENERATED_FILE}.d -o ${GENERATED_FILE} -c ${SYCL_FILE}
        -I${HEADER_PATH} ${SYCL_SYSTEM_INCLUDE_FLAGS}
      DEPENDS "${SYCL_FILE}"
      DEPFILE "${GENERATED_FILE}.d"
      COMMENT "Compiling SYCL source ${SYCL_FILE} with dpcpp")

    list(APPEND GENERATED_CPP_SOURCES ${GENERATED_FILE})
  endforeach()
  add_custom_target(sycl_targets ALL DEPENDS "${GENERATED_CPP_SOURCES}")
  set(${OUTPUT}
      "${GENERATED_CPP_SOURCES}"
      PARENT_SCOPE)
endfunction()
