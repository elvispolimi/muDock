cmake_minimum_required(VERSION 3.25 FATAL_ERROR)
cmake_policy(SET CMP0167 NEW)

# Utilities
# define the function to remove the version from the end of shared object
function(remove_version_from_filename input_filename version output_filename)
    # Extract the directory path and filename without extension
    get_filename_component(directory_path ${input_filename} DIRECTORY)
    get_filename_component(filename_without_extension ${input_filename} NAME_WE)

    # Reconstruct the modified filename with the directory path
    set(${output_filename} "${directory_path}/${filename_without_extension}.so" PARENT_SCOPE)
endfunction()

# ########################################################################
# ### Define the compile options
# ########################################################################
# Accelerator
set(MUDOCK_GPU_ARCHITECTURES "none" CACHE STRING "Optimize for specific target GPU architecture (none to disable, otherwise cuda:sm_$$ or hip:$$)")
set(MUDOCK_CPU_ARCHITECTURES "native" CACHE STRING "Optimize for specific target GPU architecture (default is native, none to disable optimizations, or specify the CPU architecture)")


if(NOT MUDOCK_GPU_ARCHITECTURES STREQUAL "none")
    # Find the position of the colon
    string(FIND ${MUDOCK_GPU_ARCHITECTURES} ":" COLON_POS)

    # Separate the two parts
    if(COLON_POS GREATER_EQUAL 0)
        # Get the part before the colon
        string(SUBSTRING ${MUDOCK_GPU_ARCHITECTURES} 0 ${COLON_POS} GPU_PLATFORM)

        # Get the part after the colon
        string(SUBSTRING ${MUDOCK_GPU_ARCHITECTURES} ${COLON_POS} 1000 GPU_TEMP)
        string(SUBSTRING ${GPU_TEMP} 1 1000 GPU_ARCHITECTURE) # Skip the colon
    else()
        message(FATAL_ERROR "${MUDOCK_GPU_ARCHITECTURES} cannot be parsed.")
    endif()
endif()

# HIP
option(MUDOCK_ENABLE_HIP "Enable the HIP kernel implementation" OFF)

if(MUDOCK_ENABLE_HIP)
    if(NOT MUDOCK_ENABLE_POLY AND NOT MUDOCK_ENABLE_SYCL AND NOT MUDOCK_ENABLE_CUDA AND NOT MUDOCK_ENABLE_OMP)
        # enable_language(HIP)
        find_program(HIP_COMPILER hipcc)
        set(CMAKE_CXX_COMPILER "${HIP_COMPILER}")

        if(GPU_PLATFORM STREQUAL "cuda")
            set(HIP_PLATFORM "nvidia")
            add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-arch=${GPU_ARCHITECTURE}>)
        elseif(GPU_PLATFORM STREQUAL "amd")
            set(HIP_PLATFORM "${GPU_PLATFORM}")

            # TODO target arch
        else()
            # TODO add other HIP devices support
            message(ERROR "Wrong target for HIP implementation")
        endif()

        set(HIP_ARCHITECTURES "${GPU_ARCHITECTURE}")
        find_package(HIP REQUIRED)
        find_package(hiprand REQUIRED)
    else()
        message(FATAL_ERROR "Enable only one parallel framework at a time.")
    endif()
endif()

# Polygeist
option(MUDOCK_ENABLE_POLY "Enable the CUDA kernel transpilation using Polygeist" OFF)
set(POLY_TARGET "cpu" CACHE STRING "Polygeist transpilation target")
set(MUDOCK_POLY_LIB_PATH "" CACHE STRING "muDock Polygeist library folder")

if(MUDOCK_ENABLE_POLY)
    if(MUDOCK_POLY_LIB_PATH STREQUAL "")
        message(FATAL_ERROR "muDock Polygeist folder lib not defined, define it with MUDOCK_POLY_LIB_PATH")
    endif()
    if(NOT MUDOCK_ENABLE_HIP AND NOT MUDOCK_ENABLE_SYCL AND NOT MUDOCK_ENABLE_CUDA)
        # TODO set the compiler to be the clang from polygeist
        if(NOT MUDOCK_GPU_ARCHITECTURES STREQUAL "none")
            # TODO try on AMD GPUs
            # TODO change POLY_TARGET
            message(FATAL_ERROR "muDock does not support Polygeist transpilation to GPUs")
        endif()
    else()
        message(FATAL_ERROR "Enable only one parallel framework at a time.")
    endif()
endif()

# CPU
# TODO check this
project(muDock VERSION 1.0)
enable_language(CXX)

# Check if specific target architectures should be used
# TODO check this
if(MUDOCK_CPU_ARCHITECTURES STREQUAL "native")
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
        if(MUDOCK_ENABLE_HIP)
            list(APPEND release_flags -Xcompiler -march=native)
        else()
            list(APPEND release_flags "-march=native")
        endif()
    else()
        message(WARNING "Unable to optimize for the native architecture: unknown compiler \"${CMAKE_CXX_COMPILER_ID}\"")
    endif()
endif()

# CUDA
option(MUDOCK_ENABLE_CUDA "Enable the CUDA kernel implementation" OFF)

if(MUDOCK_ENABLE_CUDA)
    if(NOT MUDOCK_ENABLE_HIP AND NOT MUDOCK_ENABLE_SYCL AND NOT MUDOCK_ENABLE_POLY AND NOT MUDOCK_ENABLE_OMP)
        enable_language(CUDA)

        # Find the position of the _
        string(FIND ${GPU_ARCHITECTURE} "_" UNDER_POS)

        # Separate the two parts
        if(UNDER_POS GREATER_EQUAL 0)
            # Get the part after the _
            string(SUBSTRING ${GPU_ARCHITECTURE} ${UNDER_POS} 1000 CMAKE_CUDA_TEMP)
            string(SUBSTRING ${CMAKE_CUDA_TEMP} 1 1000 CMAKE_CUDA_ARCHITECTURES) # Skip the colon
        else()
            message(FATAL_ERROR "${MUDOCK_GPU_ARCHITECTURES} cannot get CUDA GPUs.")
        endif()
    else()
        message(FATAL_ERROR "Enable only one parallel framework at a time.")
    endif()
endif()


# SYCL
option(MUDOCK_ENABLE_SYCL "Enable the SYCL kernel implementation" OFF)

if(MUDOCK_ENABLE_SYCL)
    if(NOT MUDOCK_ENABLE_HIP AND NOT MUDOCK_ENABLE_POLY AND NOT MUDOCK_ENABLE_CUDA AND NOT MUDOCK_ENABLE_OMP)
        # handle AdaptiveCPP
        # check if we need to enable the openmp backend for SYCL
        if(MUDOCK_CPU_ARCHITECTURES STREQUAL "native") # TODO what to do with native?
            if(MUDOCK_CPU_ARCHITECTURES STREQUAL "none")
                list(APPEND sycl_targets "omp")
            endif()
        else()
            list(APPEND sycl_targets "${MUDOCK_CPU_ARCHITECTURES}")
        endif()

        if(NOT MUDOCK_GPU_ARCHITECTURES STREQUAL "none")
            list(APPEND sycl_targets "${MUDOCK_GPU_ARCHITECTURES}")
        elseif()
            list(APPEND sycl_targets "omp")
        endif()

        # check if we need to enable the AMD or NVIDIA backend for SYCL
        
        # prepare the option for the AdaptiveCPP package finder
        set(SYCL_TARGETS ${sycl_targets})
        set(ACPP_TARGETS "${SYCL_TARGETS}" CACHE STRING "AdaptiveCPP targets")
        find_package(AdaptiveCpp CONFIG REQUIRED)

    # TODO OneAPI
    else()
        message(FATAL_ERROR "Enable only one parallel framework at a time.")
    endif()
endif()

# OpenMP
option(MUDOCK_ENABLE_OMP "Enable the OpenMP kernel implementation" OFF)

if(MUDOCK_ENABLE_OMP)
    if(NOT MUDOCK_ENABLE_HIP AND NOT MUDOCK_ENABLE_CUDA AND NOT MUDOCK_ENABLE_SYCL)
        find_package(OpenMP REQUIRED)
        if(GPU_PLATFORM STREQUAL "cuda" OR GPU_PLATFORM STREQUAL "amd")
            set(OMP_TARGET "${MUDOCK_GPU_ARCHITECTURES}")
        else()
            # Set to use the CPU
            set(OMP_TARGET "cpu")
        endif()

    else()
        message(FATAL_ERROR "Enable only one parallel framework at a time.")
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
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_definitions(DEBUG_MODE)
endif()

# export also compile commands to have a better integration with the ide
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ########################################################################
# ### Handle third party dependenies
# ########################################################################

# Third-party libraries that must be manually installed
find_package(Boost CONFIG REQUIRED COMPONENTS program_options graph fiber)
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
# TODO check here
if(NOT MUDOCK_ENABLE_HIP)
    list(APPEND global_c_cxx_flags -Wall -Wextra -Wpedantic -Wshadow -Wdouble-promotion)
endif()

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
message(STATUS "CXX compiler........................: ${CMAKE_CXX_COMPILER} (CMAKE_CXX_COMPILER)")
message(STATUS "Build type..........................: ${CMAKE_BUILD_TYPE} (CMAKE_BUILD_TYPE)")
message(STATUS "Install prefix......................: ${CMAKE_INSTALL_PREFIX} (CMAKE_INSTALL_PREFIX)")
message(STATUS "CPU Architecture....................: ${MUDOCK_CPU_ARCHITECTURES} (MUDOCK_CPU_ARCHITECTURES)")

message(STATUS "Enable CUDA.........................: ${MUDOCK_ENABLE_CUDA} (MUDOCK_ENABLE_CUDA)")

if(MUDOCK_ENABLE_CUDA)
    # TODO use instead the MUDOCK_GPU_ARCHITECTURES
    message(STATUS "Target CUDA architecture............: ${CMAKE_CUDA_ARCHITECTURES} (CMAKE_CUDA_ARCHITECTURES)")
endif()

message(STATUS "Enable CUDA transpilation...........: ${MUDOCK_ENABLE_POLY} (MUDOCK_ENABLE_POLY)")

if(MUDOCK_ENABLE_POLY)
    # TODO use instead the MUDOCK_GPU_ARCHITECTURES
    message(STATUS "Target Transpilation architecture...: ${POLY_TARGET} (POLY_TARGET)")
endif()

message(STATUS "Enable HIP..........................: ${MUDOCK_ENABLE_HIP} (MUDOCK_ENABLE_HIP)")

if(MUDOCK_ENABLE_HIP)
    message(STATUS "Target HIP platform.................: ${HIP_PLATFORM} (HIP_PLATFORM)")
    message(STATUS "Target HIP targets..................: ${HIP_ARCHITECTURES} (HIP_ARCHITECTURES)")
endif()

message(STATUS "Enable SYCL.........................: ${MUDOCK_ENABLE_SYCL} (MUDOCK_ENABLE_SYCL)")

if(MUDOCK_ENABLE_SYCL)
    message(STATUS "Target SYCL architecture............: ${ACPP_TARGETS} (ACPP_TARGETS)")
endif()

message(STATUS "Enable OpenMP.......................: ${MUDOCK_ENABLE_OMP} (MUDOCK_ENABLE_OMP)")

if(MUDOCK_ENABLE_OMP)
    message(STATUS "Target OpenMP architecture..........: ${OMP_TARGET} (OMP_TARGET)")
endif()

message(STATUS "//===---------------------------===//")
message(STATUS "")
