# ##############################################################################
# HIP Helper Functions #
# ##############################################################################

function(hip_separable_compilation target)
  if(CMAKE_HIP_PLATFORM STREQUAL "nvidia")
    set(MUDOCK_LIB "$<TARGET_FILE:libmudock>")
    set(OBJECTS $<TARGET_OBJECTS:${target}>)
    set(DEVICE_LINK_OBJ ${CMAKE_CURRENT_BINARY_DIR}/${target}.o)

    separate_arguments(flags UNIX_COMMAND PROGRAM ${CMAKE_HIP_FLAGS})

    string(REPLACE ";" " " OBJECTS_CLEAN "${OBJECTS}")
    separate_arguments(OBJECTS_ARGS UNIX_COMMAND PROGRAM ${OBJECTS_CLEAN})
    # FIXME check if -D defintions are correctly set
    add_custom_command(
      OUTPUT ${DEVICE_LINK_OBJ}
      COMMAND ${CMAKE_HIP_COMPILER} ${flags} -shared -dlink -o
              ${DEVICE_LINK_OBJ} ${OBJECTS_ARGS} ${MUDOCK_LIB}
      DEPENDS ${OBJECTS_ARGS} libmudock
      VERBATIM)
    add_custom_target(${target}_device_link_obj ALL DEPENDS ${DEVICE_LINK_OBJ})
    add_dependencies(${target} ${target}_device_link_obj)
    target_sources("${target}" PRIVATE ${DEVICE_LINK_OBJ})
    # target_link_libraries("${target}" PRIVATE ${DEVICE_LINK_OBJ})
  endif()
endfunction()

function(add_executable_mudock target)
  add_executable(${target} ${ARGN})
  if(MUDOCK_ENABLE_HIP)
    hip_separable_compilation(${target})
  endif()
endfunction()
