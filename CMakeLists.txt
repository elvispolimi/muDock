cmake_minimum_required(VERSION 3.25 FATAL_ERROR)
project(muDock VERSION 1.0)
enable_language(CXX)

# ########################################################################
# ### Define the compile options
# ########################################################################
set(MUDOCK_CPU_ARCHITECTURES "none" CACHE STRING "Optimize for specific target CPU architecture (generic to disable)")
set(MUDOCK_GPU_ARCHITECTURES "none" CACHE STRING "Optimize for specific target GPU architecture (none to disable, otherwise cuda:sm_$$ or hip:$$)")

option(MUDOCK_ENABLE_CUDA "Enable the CUDA kernel implementation" OFF)

if(MUDOCK_ENABLE_CUDA)
    enable_language(CUDA)
endif()

option(MUDOCK_ENABLE_SYCL "Enable the SYCL kernel implementation" OFF)

if(MUDOCK_ENABLE_SYCL)
    # handle AdaptiveCPP
    # check if we need to enable the openmp backend for SYCL
    if(NOT MUDOCK_CPU_ARCHITECTURES STREQUAL "none")
        if(NOT MUDOCK_CPU_ARCHITECTURES STREQUAL "generic")
            list(APPEND sycl_targets "omp")
        else()
            list(APPEND sycl_targets "${MUDOCK_CPU_ARCHITECTURES}")
        endif()
    endif()

    if(NOT MUDOCK_GPU_ARCHITECTURES STREQUAL "none")
        list(APPEND sycl_targets "${MUDOCK_GPU_ARCHITECTURES}")
    endif()

    # check if we need to enable the AMD or NVIDIA backend for SYCL

    # prepare the option for the AdaptiveCPP package finder
    list(JOIN sycl_targets ";" SYCL_TARGETS)
    set(ACPP_TARGETS "${SYCL_TARGETS}" CACHE STRING "AdaptiveCPP targets")
    find_package(AdaptiveCpp REQUIRED)

    # TODO OneAPI
endif()

# Check if specific target architectures should be used
# TODO check this
if(NOT MUDOCK_CPU_ARCHITECTURE STREQUAL "generic")
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        list(APPEND release_flags "-march=native")

        if(CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "ppc64le")
            list(APPEND release_flags "-mcpu=native")
        endif()
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
        list(APPEND release_flags "-march=native")
    else()
        message(WARNING "Unable to optimize for the native architecture: unknown compiler \"${CMAKE_CXX_COMPILER_ID}\"")
    endif()
endif()

# ########################################################################
# ### Tweak the default configuration
# ########################################################################

# the idea is that we want to build the optimized version, unless the user says otherwise
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "The type of build" FORCE)
    message(STATUS "Setting build type to '${CMAKE_BUILD_TYPE}' as none was specified")
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo")
endif()

# export also compile commands to have a better integration with the ide
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ########################################################################
# ### Handle third party dependenies
# ########################################################################

# Third-party libraries that must be manually installed
find_package(Boost CONFIG REQUIRED COMPONENTS program_options graph)
find_package(rdkit REQUIRED)

# Third-party libraries that we download and compile
include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/prerequisite.cmake")

# ########################################################################
# ### Enable global compiler flags
# ########################################################################
# NOTE: we don't want to add those flags as target dependency, so we
# declare them as global, which affect only this directory and
# all the child ones, but not the parents

# add more warnings (these flags are compatible with gcc and clang)
list(APPEND global_c_cxx_flags -Wall -Wextra -Wpedantic -Wshadow -Wdouble-promotion)

# if we are using ninja, force to have colored output
if("${CMAKE_GENERATOR}" STREQUAL "Ninja")
    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        list(APPEND global_c_cxx_flags -fdiagnostics-color=always)
    elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        list(APPEND global_c_cxx_flags -fcolor-diagnostics)
    endif()
endif()

# add all the previously defined flags to each target language
foreach(flag IN LISTS global_c_cxx_flags)
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:${flag}>)
endforeach()

# ########################################################################
# ### Add each component that compose this project
# ########################################################################
add_subdirectory(application)
add_subdirectory(mudock)

# ########################################################################
# ### Print a small recap of the project configuration
# ########################################################################
message(STATUS "")
message(STATUS "//===---------------------------===//")
message(STATUS "//       CONFIGURATION RECAP       //")
message(STATUS "//===---------------------------===//")
message(STATUS "Build type..........................: ${CMAKE_BUILD_TYPE} (CMAKE_BUILD_TYPE)")
message(STATUS "Install prefix......................: ${CMAKE_INSTALL_PREFIX} (CMAKE_INSTALL_PREFIX)")

message(STATUS "Enable CUDA.........................: ${MUDOCK_ENABLE_CUDA} (MUDOCK_ENABLE_CUDA)")

if(MUDOCK_ENABLE_CUDA)
    # TODO use instead the MUDOCK_GPU_ARCHITECTURES
    message(STATUS "Target CUDA architecture............: ${CMAKE_CUDA_ARCHITECTURES} (CMAKE_CUDA_ARCHITECTURES)")
endif()

message(STATUS "Enable SYCL.........................: ${MUDOCK_ENABLE_SYCL} (MUDOCK_ENABLE_SYCL)")

if(MUDOCK_ENABLE_SYCL)
    message(STATUS "Target SYCL architecture............: ${MUDOCK_GPU_ARCHITECTURES} (MUDOCK_GPU_ARCHITECTURES)")
endif()

message(STATUS "//===---------------------------===//")
message(STATUS "")
