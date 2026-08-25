cmake_minimum_required(VERSION 3.25 FATAL_ERROR)
if(POLICY CMP0144)
  cmake_policy(SET CMP0144 NEW)
endif()
if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()

# Utilities define the function to remove the version from the end of shared
# object
function(remove_version_from_filename input_filename version output_filename)
  # Extract the directory path and filename without extension
  get_filename_component(directory_path ${input_filename} DIRECTORY)
  get_filename_component(filename_without_extension ${input_filename} NAME_WE)

  # Reconstruct the modified filename with the directory path
  set(${output_filename}
      "${directory_path}/${filename_without_extension}.so"
      PARENT_SCOPE)
endfunction()

# Inlude compialtios details
include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/mudock.cmake")

# ##############################################################################
# ### Define the compile options
# ##############################################################################
# add the flag -g for IDE debugging
add_compile_options(-g) #-O0 to be tested but takes a lot of time
# Accelerator
set(MUDOCK_GPU_ARCHITECTURES
    "none"
    CACHE
      STRING
      "Optimize for target GPU architecture (none, nvidia:sm_$$, amd:gfx$$$, intel:$$; cuda:sm_$$ kept as alias for nvidia)"
)
set(MUDOCK_CPU_TARGET
    "native"
    CACHE
      STRING
      "Optimize for specific target CPU (default is empty, no optimizations)")
set(MUDOCK_CPU_ARCHITECTURES
    "native"
    CACHE STRING
          "Optimize for specific target CPU (default is empty, no optimizations"
)
set(MUDOCK_CPU_TUNE
    "native"
    CACHE STRING
          "Tune for specific target CPU (default is empty, no optimizations")

include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/sycl.cmake")

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

    # Backward-compatibility alias: treat "cuda:*" as NVIDIA target naming.
    if(GPU_PLATFORM STREQUAL "cuda")
      set(GPU_PLATFORM "nvidia")
    endif()
  else()
    message(FATAL_ERROR "${MUDOCK_GPU_ARCHITECTURES} cannot be parsed.")
  endif()
endif()

# CPU
project(muDock VERSION 1.0)
enable_language(CXX)

function(mudock_normalize_stage_bucket_policy out_var input_policy)
  set(policy "${input_policy}")
  set(${out_var} "${policy}" PARENT_SCOPE)
endfunction()

set(MUDOCK_DEFAULT_ATOM_CLUSTER_LEVEL "OFF")
set(MUDOCK_DEFAULT_STAGE_BUCKET_POLICY "MAX_UTILIZATION")
set(MUDOCK_DEFAULT_GENETIC_BUCKET_COMBINE "MIN")
set(MUDOCK_DEFAULT_UNROLL_FACTOR "8")

set(MUDOCK_ATOM_CLUSTER_LEVEL "${MUDOCK_DEFAULT_ATOM_CLUSTER_LEVEL}" CACHE STRING
    "Ligand atom-cluster granularity: OFF|MEDIUM|LARGE|EXTREME")
set_property(CACHE MUDOCK_ATOM_CLUSTER_LEVEL PROPERTY STRINGS OFF MEDIUM LARGE EXTREME)
option(MUDOCK_ENABLE_IPO "Enable IPO" OFF)
if(MUDOCK_ENABLE_IPO)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT ipo_supported OUTPUT ipo_output)
  if(ipo_supported)
    message(STATUS "IPO is enabled.")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
  else()
    message(WARNING "IPO is not supported: ${ipo_output}")
  endif()
endif()

option(
  OVERRIDE_PRECOMPUTE_OPT
  "Compile selected precompute-heavy sources with -O2"
  OFF)

#MPI
option(MUDOCK_ENABLE_MPI "Enable MPI support" OFF)
if(MUDOCK_ENABLE_MPI)
  find_package(MPI REQUIRED COMPONENTS CXX)
endif()

# HIP
option(MUDOCK_ENABLE_HIP "Enable the HIP kernel implementation" OFF)

if(MUDOCK_ENABLE_HIP)
  if(NOT MUDOCK_ENABLE_SYCL
     AND NOT MUDOCK_ENABLE_CUDA
     AND NOT MUDOCK_ENABLE_OMP_OFFLOAD)
    find_program(HIP_COMPILER hipcc)

    if(GPU_PLATFORM STREQUAL "nvidia")
      set(HIP_PLATFORM "nvidia")
      find_package(CUDAToolkit REQUIRED COMPONENTS curand)
    elseif(GPU_PLATFORM STREQUAL "amd")
      set(HIP_PLATFORM "${GPU_PLATFORM}")
      find_package(rocrand REQUIRED)

      # TODO target arch
    else()
      # TODO add other HIP devices support
      message(FATAL_ERROR "Wrong target for HIP implementation")
    endif()

    set(HIP_ARCHITECTURE "${GPU_ARCHITECTURE}")
    if(HIP_PLATFORM STREQUAL "nvidia")
      string(REGEX REPLACE "^sm_" "" HIP_ARCHITECTURE "${HIP_ARCHITECTURE}")
    endif()
    set(CMAKE_HIP_ARCHITECTURES "${HIP_ARCHITECTURE}")
    find_package(hip REQUIRED)
    enable_language(HIP)
    find_package(hiprand REQUIRED)
  else()
    message(FATAL_ERROR "Enable only one parallel framework at a time.")
  endif()
endif()
# Check if specific target architectures should be used TODO check this
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU"
   OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang"
   OR CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
  if(MUDOCK_ENABLE_HIP)
    list(APPEND release_flags -march=native)
  else()
    list(APPEND release_flags "-march=${MUDOCK_CPU_ARCHITECTURES}")
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64")
      list(APPEND release_flags "-mcpu=${MUDOCK_CPU_TARGET}")
    endif()
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
      list(APPEND release_flags "-mtune=${MUDOCK_CPU_TUNE}")
    endif()
  endif()
else()
  message(
    WARNING
      "Unable to optimize for the given architecture: unknown compiler \"${CMAKE_CXX_COMPILER_ID}\""
  )
endif()

# Google Highway
option(MUDOCK_ENABLE_GH "Enable Google Highway for CPU vectorization" OFF)
if(MUDOCK_ENABLE_GH)
  find_package(HWY REQUIRED)
endif()

# XSimd
option(MUDOCK_ENABLE_XSIMD "Enable XSimd for CPU vectorization" OFF)
if(MUDOCK_ENABLE_XSIMD)
  find_package(xsimd REQUIRED CONFIG)
  set(XSIMD_INCLUDE_DIR "${xsimd_DIR}/../../../include")
endif()

# CUDA
option(MUDOCK_ENABLE_CUDA "Enable the CUDA kernel implementation" OFF)
option(MUDOCK_ENABLE_KERNEL_LOCK "Serialize kernel launches per device" OFF)

if(MUDOCK_ENABLE_CUDA)
  if(NOT MUDOCK_ENABLE_HIP
     AND NOT MUDOCK_ENABLE_SYCL
     AND NOT MUDOCK_ENABLE_OMP_OFFLOAD)
    enable_language(CUDA)
    # Compile commands for CUDA includes
    # https://discourse.cmake.org/t/cmake-target-include-directories-doesnt-export-to-compile-commands-json-when-using-cuda/10072/9
    set(CMAKE_CUDA_USE_RESPONSE_FILE_FOR_INCLUDES 0)
    set(CMAKE_CUDA_USE_RESPONSE_FILE_FOR_LIBRARIES 0)
    set(CMAKE_CUDA_USE_RESPONSE_FILE_FOR_OBJECTS 0)

    # Find the position of the _
    string(FIND ${GPU_ARCHITECTURE} "_" UNDER_POS)

    # Separate the two parts
    if(UNDER_POS GREATER_EQUAL 0)
      # Get the part after the _
      string(SUBSTRING ${GPU_ARCHITECTURE} ${UNDER_POS} 1000 CMAKE_CUDA_TEMP)
      string(SUBSTRING ${CMAKE_CUDA_TEMP} 1 1000 CUDA_ARCHITECTURES)
      set(CMAKE_CUDA_ARCHITECTURES
          "${CUDA_ARCHITECTURES}"
          CACHE STRING "CUDA architectures to target (e.g., 75 for Turing)"
                FORCE)
    else()
      message(FATAL_ERROR "${MUDOCK_GPU_ARCHITECTURES} cannot get CUDA GPUs.")
    endif()
  else()
    message(FATAL_ERROR "Enable only one parallel framework at a time.")
  endif()
endif()

# SYCL
option(MUDOCK_ENABLE_SYCL "Enable the SYCL kernel implementation" OFF)
set(MUDOCK_SYCL_EXTRA_FLAGS "" CACHE STRING "Extra SYCL flags appended to both oneAPI compile and link steps")
set(MUDOCK_SYCL_EXTRA_COMPILE_FLAGS "" CACHE STRING "Extra SYCL flags appended only to oneAPI compile steps")
set(MUDOCK_SYCL_EXTRA_LINK_FLAGS "" CACHE STRING "Extra SYCL link flags appended to the clang++ link step")
set(MUDOCK_SYC_EXTRA_LINK_FLAGS "" CACHE STRING "Deprecated alias for MUDOCK_SYCL_EXTRA_LINK_FLAGS")
set(MUDOCK_SYCL_WG_SIZE "" CACHE STRING "SYCL work-group size override for kernels (e.g. 32 or 64)")

set(MUDOCK_STAGE_BUCKET_OVERRIDE "" CACHE STRING
    "Override stage bucket size (optional integer)")
set(MUDOCK_STAGE_BUCKET_MULTIPLE_OVERRIDE "" CACHE STRING
    "Override stage bucket scaling factor (optional real > 0)")
set(MUDOCK_STAGE_BUCKET_POLICY "${MUDOCK_DEFAULT_STAGE_BUCKET_POLICY}" CACHE STRING
    "Stage bucket policy: DEVICE_ALIGNED|SM_ALIGNED|MAX_UTILIZATION")
set_property(CACHE MUDOCK_STAGE_BUCKET_POLICY PROPERTY STRINGS DEVICE_ALIGNED
                                                          SM_ALIGNED
                                                          MAX_UTILIZATION)
set(MUDOCK_GENETIC_BUCKET_COMBINE "${MUDOCK_DEFAULT_GENETIC_BUCKET_COMBINE}" CACHE STRING "GENETIC score/geom combine policy: MIN|SCORE_ONLY|GEOM_ONLY|LCM")
set_property(CACHE MUDOCK_GENETIC_BUCKET_COMBINE PROPERTY STRINGS MIN SCORE_ONLY GEOM_ONLY LCM)
option(MUDOCK_ENABLE_STAGE_BUCKET_TRACE
       "Enable detailed runtime logging for stage bucket selection"
       OFF)
option(MUDOCK_DISABLE_UNROLL "Disable kernel loop unrolling" OFF)
set(MUDOCK_UNROLL_FACTOR "${MUDOCK_DEFAULT_UNROLL_FACTOR}" CACHE STRING "Kernel loop unroll factor when unrolling is enabled")

if(MUDOCK_ENABLE_SYCL)
  if(NOT MUDOCK_ENABLE_HIP
     AND NOT MUDOCK_ENABLE_CUDA
     AND NOT MUDOCK_ENABLE_OMP_OFFLOAD)
    if(MUDOCK_SYCL_WG_SIZE STREQUAL "")
      if(MUDOCK_GPU_ARCHITECTURES MATCHES "amd")
        set(MUDOCK_SYCL_WG_SIZE 64)
      elseif(MUDOCK_GPU_ARCHITECTURES MATCHES "intel")
        set(MUDOCK_SYCL_WG_SIZE 32)
      else()
        # NVIDIA and others default to 32
        set(MUDOCK_SYCL_WG_SIZE 32)
      endif()
    endif()
    find_package(LLVM REQUIRED)
    set(CMAKE_LINKER
        clang++
        CACHE STRING "C++ linker" FORCE)
    # Optional user-provided compile flags for SYCL source compilation.
    # If passed as -DSYCL_BACKEND_FLAGS_COMPILE="...", append to auto backend flags.
    set(SYCL_BACKEND_FLAGS_COMPILE_USER "${SYCL_BACKEND_FLAGS_COMPILE}")
    set(MUDOCK_SYCL_EXTRA_FLAGS_USER "${MUDOCK_SYCL_EXTRA_FLAGS}")
    set(MUDOCK_SYCL_EXTRA_COMPILE_FLAGS_USER "${MUDOCK_SYCL_EXTRA_COMPILE_FLAGS}")
    set(SYCL_BACKEND_FLAGS_LINK "")
    set(SYCL_BACKEND_FLAGS_COMPILE "")
    set(MUDOCK_SYCL_EXTRA_LINK_FLAGS_EFFECTIVE "${MUDOCK_SYCL_EXTRA_LINK_FLAGS}")
    if(NOT MUDOCK_SYC_EXTRA_LINK_FLAGS STREQUAL "")
      if(NOT MUDOCK_SYCL_EXTRA_LINK_FLAGS_EFFECTIVE STREQUAL "")
        string(APPEND MUDOCK_SYCL_EXTRA_LINK_FLAGS_EFFECTIVE " ")
      endif()
      string(APPEND MUDOCK_SYCL_EXTRA_LINK_FLAGS_EFFECTIVE
                    "${MUDOCK_SYC_EXTRA_LINK_FLAGS}")
      message(
        WARNING
          "MUDOCK_SYC_EXTRA_LINK_FLAGS is deprecated. Use MUDOCK_SYCL_EXTRA_LINK_FLAGS instead."
      )
    endif()
    if(NOT MUDOCK_SYCL_EXTRA_FLAGS STREQUAL "")
      if(NOT MUDOCK_SYCL_EXTRA_LINK_FLAGS_EFFECTIVE STREQUAL "")
        string(APPEND MUDOCK_SYCL_EXTRA_LINK_FLAGS_EFFECTIVE " ")
      endif()
      string(APPEND MUDOCK_SYCL_EXTRA_LINK_FLAGS_EFFECTIVE
                    "${MUDOCK_SYCL_EXTRA_FLAGS}")
    endif()
    string(FIND "${MUDOCK_GPU_ARCHITECTURES}" ":" COLON_INDEX)

    # Check that a colon was found
    if(NOT COLON_INDEX EQUAL -1)
      math(EXPR VALUE_INDEX "${COLON_INDEX} + 1")
      string(SUBSTRING "${MUDOCK_GPU_ARCHITECTURES}" ${VALUE_INDEX} -1
                       VALUE_AFTER_COLON)
      if(MUDOCK_GPU_ARCHITECTURES MATCHES "nvidia")
        list(APPEND sycl_targets "nvptx64-nvidia-cuda")
        if(VALUE_AFTER_COLON MATCHES "^sm_")
          set(CUDA_ARCH "${VALUE_AFTER_COLON}")
        else()
          set(CUDA_ARCH "sm_${VALUE_AFTER_COLON}")
        endif()
        set(SYCL_BACKEND_FLAGS_LINK
            "-Xsycl-target-backend=nvptx64-nvidia-cuda;--cuda-gpu-arch=${CUDA_ARCH};-Xcuda-ptxas;-v")
        set(SYCL_BACKEND_FLAGS_COMPILE "${SYCL_BACKEND_FLAGS_LINK}")
      elseif(MUDOCK_GPU_ARCHITECTURES MATCHES "amd")
        list(APPEND sycl_targets "amdgcn-amd-amdhsa")
        set(SYCL_BACKEND_FLAGS_LINK
            "-Xsycl-target-backend=amdgcn-amd-amdhsa;--offload-arch=${VALUE_AFTER_COLON};-Xsycl-target-frontend=amdgcn-amd-amdhsa;-Rpass-analysis=kernel-resource-usage")
        set(SYCL_BACKEND_FLAGS_COMPILE "${SYCL_BACKEND_FLAGS_LINK}")
      elseif(MUDOCK_GPU_ARCHITECTURES MATCHES "intel")
        list(APPEND sycl_targets "intel_gpu_${VALUE_AFTER_COLON}")
        set(SYCL_BACKEND_FLAGS_COMPILE "${SYCL_BACKEND_FLAGS_LINK}")
      endif()
      if(NOT SYCL_BACKEND_FLAGS_COMPILE_USER STREQUAL "")
        separate_arguments(SYCL_BACKEND_FLAGS_COMPILE_USER_LIST UNIX_COMMAND
                           "${SYCL_BACKEND_FLAGS_COMPILE_USER}")
        list(APPEND SYCL_BACKEND_FLAGS_COMPILE ${SYCL_BACKEND_FLAGS_COMPILE_USER_LIST})
      endif()
      if(NOT MUDOCK_SYCL_EXTRA_FLAGS_USER STREQUAL "")
        separate_arguments(MUDOCK_SYCL_EXTRA_FLAGS_USER_LIST UNIX_COMMAND
                           "${MUDOCK_SYCL_EXTRA_FLAGS_USER}")
        list(APPEND SYCL_BACKEND_FLAGS_COMPILE
             ${MUDOCK_SYCL_EXTRA_FLAGS_USER_LIST})
      endif()
      if(NOT MUDOCK_SYCL_EXTRA_COMPILE_FLAGS_USER STREQUAL "")
        separate_arguments(MUDOCK_SYCL_EXTRA_COMPILE_FLAGS_USER_LIST UNIX_COMMAND
                           "${MUDOCK_SYCL_EXTRA_COMPILE_FLAGS_USER}")
        list(APPEND SYCL_BACKEND_FLAGS_COMPILE
             ${MUDOCK_SYCL_EXTRA_COMPILE_FLAGS_USER_LIST})
      endif()
      string(JOIN " " SYCL_BACKEND_FLAGS_LINK_STR ${SYCL_BACKEND_FLAGS_LINK})
      set(CMAKE_CXX_LINK_EXECUTABLE
          "${LLVM_TOOLS_BINARY_DIR}/clang++ -fsycl -fsycl-targets=${sycl_targets} <FLAGS> <CMAKE_CXX_LINK_FLAGS> <LINK_FLAGS> ${SYCL_BACKEND_FLAGS_LINK_STR} ${MUDOCK_SYCL_EXTRA_LINK_FLAGS_EFFECTIVE} <OBJECTS> -o <TARGET> <LINK_LIBRARIES>"
      )

    else()
      message(WARNING "No colon found in string.")
    endif()

    set(SYCL_TARGETS ${sycl_targets})
  else()
    message(FATAL_ERROR "Enable only one parallel framework at a time.")
  endif()
endif()

mudock_normalize_stage_bucket_policy(MUDOCK_STAGE_BUCKET_POLICY
                                     "${MUDOCK_STAGE_BUCKET_POLICY}")
if(NOT MUDOCK_STAGE_BUCKET_POLICY STREQUAL "MAX_UTILIZATION"
   AND NOT MUDOCK_STAGE_BUCKET_POLICY STREQUAL "SM_ALIGNED"
   AND NOT MUDOCK_STAGE_BUCKET_POLICY STREQUAL "DEVICE_ALIGNED")
  message(FATAL_ERROR "MUDOCK_STAGE_BUCKET_POLICY must be one of: DEVICE_ALIGNED, SM_ALIGNED, MAX_UTILIZATION")
endif()
if(NOT MUDOCK_ATOM_CLUSTER_LEVEL STREQUAL "OFF"
   AND NOT MUDOCK_ATOM_CLUSTER_LEVEL STREQUAL "MEDIUM"
   AND NOT MUDOCK_ATOM_CLUSTER_LEVEL STREQUAL "LARGE"
   AND NOT MUDOCK_ATOM_CLUSTER_LEVEL STREQUAL "EXTREME")
  message(FATAL_ERROR "MUDOCK_ATOM_CLUSTER_LEVEL must be one of: OFF, MEDIUM, LARGE, EXTREME")
endif()
if(NOT MUDOCK_GENETIC_BUCKET_COMBINE STREQUAL "MIN"
   AND NOT MUDOCK_GENETIC_BUCKET_COMBINE STREQUAL "SCORE_ONLY"
   AND NOT MUDOCK_GENETIC_BUCKET_COMBINE STREQUAL "GEOM_ONLY"
   AND NOT MUDOCK_GENETIC_BUCKET_COMBINE STREQUAL "LCM")
  message(FATAL_ERROR "MUDOCK_GENETIC_BUCKET_COMBINE must be one of: MIN, SCORE_ONLY, GEOM_ONLY, LCM")
endif()

# OpenMP
option(MUDOCK_ENABLE_OMP "Enable OpenMP" OFF)
option(MUDOCK_ENABLE_OMP_OFFLOAD "Enable OpenMP kernel offload" OFF)
option(
  MUDOCK_ENABLE_FAST
  "Enable the various compilers optimizimization, like OpenMP SIMD optimizations and fast math"
  OFF)

set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "The type of build" FORCE)
  message(STATUS "Setting build type to '${CMAKE_BUILD_TYPE}' as none was specified")
endif()

# If user hasn't set it AND build type is Release, set it to ON
if(NOT MUDOCK_ENABLE_FAST AND CMAKE_BUILD_TYPE STREQUAL "Release")
  set(MUDOCK_ENABLE_FAST ON)
endif()

if(MUDOCK_ENABLE_OMP OR MUDOCK_ENABLE_OMP_OFFLOAD)
  find_package(OpenMP REQUIRED)
endif()

if(MUDOCK_ENABLE_OMP_OFFLOAD)
  if(NOT MUDOCK_ENABLE_HIP
     AND NOT MUDOCK_ENABLE_CUDA
     AND NOT MUDOCK_ENABLE_SYCL)
    if(GPU_PLATFORM STREQUAL "nvidia" OR GPU_PLATFORM STREQUAL "amd")
      set(OMP_TARGET "${MUDOCK_GPU_ARCHITECTURES}")
    else()
      # Set to use the CPU
      set(OMP_TARGET "cpu")
    endif()

  else()
    message(FATAL_ERROR "Enable only one parallel framework at a time.")
  endif()
endif()

# ##############################################################################
# ### Profiling
# ##############################################################################
option(MUDOCK_ENABLE_LIKWID "Enable profiling using likwid" OFF)
option(MUDOCK_LIKWID_CPP "Enable profiling CPP using likwid" OFF)
option(MUDOCK_LIKWID_CUDA "Enable profiling CUDA using likwid" OFF)
set(_MUDOCK_LIKWID_BACKENDS 0)
if(MUDOCK_ENABLE_LIKWID)
  find_package(likwid REQUIRED)
  if (MUDOCK_LIKWID_CPP)
    math(EXPR _MUDOCK_LIKWID_BACKENDS "${_MUDOCK_LIKWID_BACKENDS}+1")
  endif()
  if (MUDOCK_LIKWID_CUDA)
    message(FATAL_ERROR "CUDA LIKWID profiler not supported yet")
    math(EXPR _MUDOCK_LIKWID_BACKENDS "${_MUDOCK_LIKWID_BACKENDS}+1")
  endif()
  if (_MUDOCK_LIKWID_BACKENDS EQUAL 0)
    message(STATUS
    "No LIKWID backend defined, switching to CPP default"
  )
    set(MUDOCK_LIKWID_CPP ON)
  endif()
  if (_MUDOCK_LIKWID_BACKENDS GREATER 1)
    message(FATAL_ERROR
    "Only one LIKWID backend may be enabled at a time.\n"
    "Enable exactly one of:\n"
    "  MUDOCK_LIKWID_CPU\n"
    "  MUDOCK_LIKWID_CUDA"
  )
  endif()
endif()

# ##############################################################################
# ### Testing
# ##############################################################################
option(MUDOCK_ENABLE_TEST "Enabel muDock tests" OFF)
# Exit code reported by test binaries when CTest should mark them as skipped.
set(MUDOCK_CTEST_SKIP_RETURN_CODE 77)
if(MUDOCK_ENABLE_TEST)
  if(CMAKE_BUILD_TYPE STREQUAL "Release")
    message(
      FATAL_ERROR
        "Cannot enable Release mode and tests due to fast-math and O3 approximations!"
    )
  endif()
  enable_testing()
  add_subdirectory(test)
endif()

# ##############################################################################
# ### Tweak the default configuration
# ##############################################################################

# the idea is that we want to build the optimized version, unless the user says
# otherwise FIX ME native is not supported by mcpu
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE STRING "The type of build" FORCE)
  message(
    STATUS "Setting build type to '${CMAKE_BUILD_TYPE}' as none was specified")
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
                                               "RelWithDebInfo")
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_compile_definitions(DEBUG_MODE)
endif()

# export also compile commands to have a better integration with the ide
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ##############################################################################
# ### Handle third party dependenies
# ##############################################################################

# Third-party libraries that must be manually installed
find_package(Boost QUIET CONFIG COMPONENTS program_options graph)
if(NOT Boost_FOUND)
  find_package(Boost REQUIRED COMPONENTS program_options graph)
endif()
find_package(Eigen3 NO_MODULE) # pre-req of OpenBabel
find_package(OpenBabel3 REQUIRED)
set(CMAKE_ENABLE_EXPORTS ON)

# ##############################################################################
# ### Enable global compiler flags
# ##############################################################################
# NOTE: we don't want to add those flags as target dependency, so we declare
# them as global, which affect only this directory and all the child ones, but
# not the parents

# add more warnings (these flags are compatible with gcc and clang) TODO check
# here
if(NOT MUDOCK_ENABLE_HIP AND (("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
                              OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang"))
  list(
    APPEND
    global_c_cxx_flags
    -Wall
    -Wextra
    -Wshadow
    -Wdouble-promotion
    -Wconversion
    )
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
add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-Wpedantic>)
foreach(flag IN LISTS global_c_cxx_flags)
  add_compile_options($<$<COMPILE_LANGUAGE:CXX>:${flag}>)
  if(MUDOCK_ENABLE_CUDA)
    add_compile_options($<$<COMPILE_LANGUAGE:CUDA>:${flag}>)
  endif()
  if(MUDOCK_ENABLE_HIP)
    add_compile_options($<$<COMPILE_LANGUAGE:HIP>:${flag}>)
  endif()
endforeach()

# ##############################################################################
# ### Add each component that compose this project
# ##############################################################################
add_subdirectory(application)
add_subdirectory(mudock)

# ##############################################################################
# ### Print a small recap of the project configuration
# ##############################################################################
function(mudock_print_config label value var_name)
  string(LENGTH "${label}" label_len)
  math(EXPR dots_len "38 - ${label_len}")
  if(dots_len LESS 1)
    set(dots_len 1)
  endif()
  string(REPEAT "." ${dots_len} dots)
  message(STATUS "${label}${dots}: ${value} (${var_name})")
endfunction()

message(STATUS "")
message(STATUS "//===---------------------------===//")
message(STATUS "//       CONFIGURATION RECAP       //")
message(STATUS "//===---------------------------===//")
mudock_print_config("CXX compiler" "${CMAKE_CXX_COMPILER}" "CMAKE_CXX_COMPILER")
mudock_print_config("Build type" "${CMAKE_BUILD_TYPE}" "CMAKE_BUILD_TYPE")
mudock_print_config("Install prefix" "${CMAKE_INSTALL_PREFIX}"
                    "CMAKE_INSTALL_PREFIX")
if(NOT MUDOCK_ATOM_CLUSTER_LEVEL STREQUAL "${MUDOCK_DEFAULT_ATOM_CLUSTER_LEVEL}")
  mudock_print_config("Atom cluster level" "${MUDOCK_ATOM_CLUSTER_LEVEL}"
                      "MUDOCK_ATOM_CLUSTER_LEVEL")
endif()
mudock_print_config("CPU Architecture" "${MUDOCK_CPU_ARCHITECTURES}"
                    "MUDOCK_CPU_ARCHITECTURES")
if(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64")
  mudock_print_config("CPU Target" "${MUDOCK_CPU_TARGET}" "MUDOCK_CPU_TARGET")
endif()
if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
  mudock_print_config("CPU Tune" "${MUDOCK_CPU_TUNE}" "MUDOCK_CPU_TUNE")
endif()

mudock_print_config("Fast mode" "${MUDOCK_ENABLE_FAST}" "MUDOCK_ENABLE_FAST")
if(OVERRIDE_PRECOMPUTE_OPT)
  mudock_print_config("Override precompute opt" "${OVERRIDE_PRECOMPUTE_OPT}"
                      "OVERRIDE_PRECOMPUTE_OPT")
endif()

if(MUDOCK_ENABLE_TEST)
  mudock_print_config("Enable Test" "${MUDOCK_ENABLE_TEST}"
                      "MUDOCK_ENABLE_TEST")
endif()

mudock_print_config("IPO Optimizations" "${MUDOCK_ENABLE_IPO}"
                    "MUDOCK_ENABLE_IPO")
mudock_print_config("Likwid Profiling" "${MUDOCK_ENABLE_LIKWID}"
                    "MUDOCK_ENABLE_LIKWID")
if(MUDOCK_DISABLE_UNROLL)
  mudock_print_config("Disable kernel unrolling" "${MUDOCK_DISABLE_UNROLL}"
                      "MUDOCK_DISABLE_UNROLL")
elseif(NOT MUDOCK_UNROLL_FACTOR STREQUAL "${MUDOCK_DEFAULT_UNROLL_FACTOR}")
  mudock_print_config("Kernel unroll factor" "${MUDOCK_UNROLL_FACTOR}"
                      "MUDOCK_UNROLL_FACTOR")
endif()
if(MUDOCK_ENABLE_STAGE_BUCKET_TRACE)
  mudock_print_config("Stage bucket trace" "${MUDOCK_ENABLE_STAGE_BUCKET_TRACE}"
                      "MUDOCK_ENABLE_STAGE_BUCKET_TRACE")
endif()
if(NOT MUDOCK_STAGE_BUCKET_OVERRIDE STREQUAL "")
  mudock_print_config("Stage bucket override" "${MUDOCK_STAGE_BUCKET_OVERRIDE}"
                      "MUDOCK_STAGE_BUCKET_OVERRIDE")
endif()
if(NOT MUDOCK_STAGE_BUCKET_MULTIPLE_OVERRIDE STREQUAL "")
  mudock_print_config("Stage bucket multiple override"
                      "${MUDOCK_STAGE_BUCKET_MULTIPLE_OVERRIDE}"
                      "MUDOCK_STAGE_BUCKET_MULTIPLE_OVERRIDE")
endif()
if(NOT MUDOCK_STAGE_BUCKET_POLICY STREQUAL "${MUDOCK_DEFAULT_STAGE_BUCKET_POLICY}")
  mudock_print_config("Stage bucket policy" "${MUDOCK_STAGE_BUCKET_POLICY}"
                      "MUDOCK_STAGE_BUCKET_POLICY")
endif()
if(NOT MUDOCK_GENETIC_BUCKET_COMBINE STREQUAL
       "${MUDOCK_DEFAULT_GENETIC_BUCKET_COMBINE}")
  mudock_print_config("Genetic stage bucket combine"
                      "${MUDOCK_GENETIC_BUCKET_COMBINE}"
                      "MUDOCK_GENETIC_BUCKET_COMBINE")
endif()
if(MUDOCK_ENABLE_LIKWID)
  mudock_print_config("Likwid CPP Profiling" "${MUDOCK_LIKWID_CPP}"
                      "MUDOCK_LIKWID_CPP")
  mudock_print_config("Likwid CUDA Profiling" "${MUDOCK_LIKWID_CUDA}"
                      "MUDOCK_LIKWID_CUDA")
endif()

mudock_print_config("Enable Google Highway" "${MUDOCK_ENABLE_GH}"
                    "MUDOCK_ENABLE_GH")
mudock_print_config("Enable XSimd" "${MUDOCK_ENABLE_XSIMD}"
                    "MUDOCK_ENABLE_XSIMD")
mudock_print_config("Enable CUDA" "${MUDOCK_ENABLE_CUDA}"
                    "MUDOCK_ENABLE_CUDA")
if(MUDOCK_ENABLE_KERNEL_LOCK)
  mudock_print_config("Kernel lock" "${MUDOCK_ENABLE_KERNEL_LOCK}"
                      "MUDOCK_ENABLE_KERNEL_LOCK")
endif()

if(MUDOCK_ENABLE_CUDA)
  mudock_print_config("Target CUDA architecture" "${CMAKE_CUDA_ARCHITECTURES}"
                      "CMAKE_CUDA_ARCHITECTURES")
endif()

mudock_print_config("Enable HIP" "${MUDOCK_ENABLE_HIP}" "MUDOCK_ENABLE_HIP")

if(MUDOCK_ENABLE_HIP)
  mudock_print_config("Target HIP platform" "${HIP_PLATFORM}" "HIP_PLATFORM")
  mudock_print_config("Target HIP targets" "${CMAKE_HIP_ARCHITECTURES}"
                      "CMAKE_HIP_ARCHITECTURES")
endif()

mudock_print_config("Enable SYCL" "${MUDOCK_ENABLE_SYCL}" "MUDOCK_ENABLE_SYCL")

if(MUDOCK_ENABLE_SYCL)
  mudock_print_config("Target SYCL architecture" "${SYCL_TARGETS}"
                      "SYCL_TARGETS")
  mudock_print_config("SYCL Work-group size" "${MUDOCK_SYCL_WG_SIZE}"
                      "MUDOCK_SYCL_WG_SIZE")
  mudock_print_config("oneAPI Compiler Path" "${LLVM_TOOLS_BINARY_DIR}"
                      "LLVM_TOOLS_BINARY_DIR")
  if(NOT SYCL_BACKEND_FLAGS_LINK STREQUAL "")
    string(JOIN " " SYCL_BACKEND_FLAGS_PRINT ${SYCL_BACKEND_FLAGS_LINK})
    mudock_print_config("oneAPI Backend Flags" "${SYCL_BACKEND_FLAGS_PRINT}"
                        "SYCL_BACKEND_FLAGS_LINK")
  endif()
  if(NOT SYCL_BACKEND_FLAGS_COMPILE_USER STREQUAL "")
    mudock_print_config("oneAPI User Compile Flags"
                        "${SYCL_BACKEND_FLAGS_COMPILE_USER}"
                        "SYCL_BACKEND_FLAGS_COMPILE")
  endif()
  if(NOT MUDOCK_SYCL_EXTRA_FLAGS STREQUAL "")
    mudock_print_config("oneAPI Extra SYCL Flags" "${MUDOCK_SYCL_EXTRA_FLAGS}"
                        "MUDOCK_SYCL_EXTRA_FLAGS")
  endif()
  if(NOT MUDOCK_SYCL_EXTRA_COMPILE_FLAGS STREQUAL "")
    mudock_print_config("oneAPI Extra Compile Flags"
                        "${MUDOCK_SYCL_EXTRA_COMPILE_FLAGS}"
                        "MUDOCK_SYCL_EXTRA_COMPILE_FLAGS")
  endif()
  if(NOT MUDOCK_SYCL_EXTRA_LINK_FLAGS STREQUAL "")
    mudock_print_config("oneAPI Extra Link Flags"
                        "${MUDOCK_SYCL_EXTRA_LINK_FLAGS}"
                        "MUDOCK_SYCL_EXTRA_LINK_FLAGS")
  endif()
  if(NOT MUDOCK_SYCL_EXTRA_LINK_FLAGS_EFFECTIVE STREQUAL "")
    mudock_print_config("oneAPI Effective Link Flags"
                        "${MUDOCK_SYCL_EXTRA_LINK_FLAGS_EFFECTIVE}"
                        "MUDOCK_SYCL_EXTRA_LINK_FLAGS_EFFECTIVE")
  endif()
  if(NOT MUDOCK_SYC_EXTRA_LINK_FLAGS STREQUAL "")
    mudock_print_config("oneAPI Deprecated Link Flags"
                        "${MUDOCK_SYC_EXTRA_LINK_FLAGS}"
                        "MUDOCK_SYC_EXTRA_LINK_FLAGS")
  endif()
endif()

mudock_print_config("Enable OpenMP" "${MUDOCK_ENABLE_OMP}" "MUDOCK_ENABLE_OMP")
mudock_print_config("Enable OpenMP Offload" "${MUDOCK_ENABLE_OMP_OFFLOAD}"
                    "MUDOCK_ENABLE_OMP_OFFLOAD")

if(MUDOCK_ENABLE_OMP_OFFLOAD)
  mudock_print_config("Target OpenMP architecture" "${OMP_TARGET}" "OMP_TARGET")
endif()

message(
  STATUS
    "Enable MPI............................: ${MUDOCK_ENABLE_MPI} (MUDOCK_ENABLE_MPI)"
)

message(STATUS "//===---------------------------===//")
message(STATUS "")
