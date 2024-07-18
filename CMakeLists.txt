cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
project(muDock VERSION 1.0)
enable_language(CXX)

#########################################################################
#### Define the compile options
#########################################################################

option(MUDOCK_ENABLE_CUDA "Enable the CUDA kernel implementation" OFF)
if (MUDOCK_ENABLE_CUDA)
    enable_language(CUDA)
endif()

#########################################################################
#### Tweak the default configuration
#########################################################################

# the idea is that we want to build the optimized version, unless the user says otherwise
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "The type of build" FORCE)
  message(STATUS "Setting build type to '${CMAKE_BUILD_TYPE}' as none was specified")
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo")
endif()

# export also compile commands to have a better integration with the ide
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

#########################################################################
#### Handle third party dependenies
#########################################################################

# Third-party libraries that must be manually installed
find_package(Boost CONFIG REQUIRED COMPONENTS program_options graph)
find_package(rdkit REQUIRED)

# Third-party libraries that we download and compile
include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/prerequisite.cmake")

#########################################################################
#### Enable global compiler flags
#########################################################################
# NOTE: we don't want to add those flags as target dependency, so we
#       declare them as global, which affect only this directory and
#       all the child ones, but not the parents

# add more warnings (these flags are compatible with gcc and clang)
list(APPEND global_c_cxx_flags -Wall -Wextra -Wpedantic -Wshadow -Wdouble-promotion)

# if we are using ninja, force to have colored output
if ("${CMAKE_GENERATOR}" STREQUAL "Ninja")
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        list(APPEND global_c_cxx_flags -fdiagnostics-color=always)
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        list(APPEND global_c_cxx_flags -fcolor-diagnostics)
    endif()
endif()

# add all the previously defined flags to each target language
foreach(flag IN LISTS global_c_cxx_flags)
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:${flag}>)
endforeach()

#########################################################################
#### Add each component that compose this project
#########################################################################

add_subdirectory(application)
add_subdirectory(mudock)

#########################################################################
#### Print a small recap of the project configuration
#########################################################################

message(STATUS "" )
message(STATUS "//===---------------------------===//" )
message(STATUS "//       CONFIGURATION RECAP       //" )
message(STATUS "//===---------------------------===//" )
message(STATUS "Build type..........................: ${CMAKE_BUILD_TYPE} (CMAKE_BUILD_TYPE)")
message(STATUS "Install prefix......................: ${CMAKE_INSTALL_PREFIX} (CMAKE_INSTALL_PREFIX)")
message(STATUS "Enable CUDA.........................: ${MUDOCK_ENABLE_CUDA} (MUDOCK_ENABLE_CUDA)")
if (MUDOCK_ENABLE_CUDA)
    message(STATUS "Target CUDA architecture............: ${CMAKE_CUDA_ARCHITECTURES} (CMAKE_CUDA_ARCHITECTURES)")
endif()
message(STATUS "//===---------------------------===//" )
message(STATUS "" )
