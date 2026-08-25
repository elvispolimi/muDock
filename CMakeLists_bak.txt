cmake_minimum_required(VERSION 3.25 FATAL_ERROR)
cmake_policy(SET CMP0167 NEW)

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

# Inlude compialtios details and Third-party libraries that we download and
# compile
include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/prerequisite.cmake")
include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/mudock.cmake")

# ##############################################################################
# ### Define the compile options
# ##############################################################################
# Accelerator
set(MUDOCK_GPU_ARCHITECTURES
    "none"
    CACHE
      STRING
      "Optimize for specific target GPU architecture (none to disable, otherwise cuda:sm_$$ or hip:$$)"
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
  else()
    message(FATAL_ERROR "${MUDOCK_GPU_ARCHITECTURES} cannot be parsed.")
  endif()
endif()

# CPU
project(muDock VERSION 1.0)
enable_language(CXX)
set(MUDOCK_BUCKET_LEVEL "OFF" CACHE STRING "Bucket level: OFF|MEDIUM|LARGE|EXTREME")
set_property(CACHE MUDOCK_BUCKET_LEVEL PROPERTY STRINGS OFF MEDIUM LARGE EXTREME)
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

if(MUDOCK_ENABLE_KERNEL_LOCK)
  add_compile_definitions(MUDOCK_KERNEL_LOCK)
endif()

# SYCL
option(MUDOCK_ENABLE_SYCL "Enable the SYCL kernel implementation" OFF)
set(MUDOCK_SYCL_EXTRA_LINK_FLAGS "" CACHE STRING "Extra SYCL link flags appended to the clang++ link step")
set(MUDOCK_SYCL_WG_SIZE "" CACHE STRING "SYCL work-group size override for kernels (e.g. 32 or 64)")
set(MUDOCK_ADT_BUCKET_OVERRIDE "" CACHE STRING "Override ADT score bucket size (optional integer)")
set(MUDOCK_ADT_BUCKET_MULTIPLE_OVERRIDE "" CACHE STRING "Override ADT score bucket multiple (optional integer)")
option(MUDOCK_DISABLE_UNROLL "Disable kernel loop unrolling" OFF)
set(MUDOCK_UNROLL_FACTOR 8 CACHE STRING "Kernel loop unroll factor when unrolling is enabled")

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
    # Only GPU targets supported for now prepare the option for the AdaptiveCPP
    # package finder
    if(${SYCL_COMPILER} STREQUAL "AdaptiveCpp")
      # TODO add ACPP generic pass support and check performance
      list(APPEND sycl_targets "${MUDOCK_GPU_ARCHITECTURES}")
      find_package(AdaptiveCpp CONFIG REQUIRED)
    elseif(${SYCL_COMPILER} STREQUAL "oneAPI")
      find_package(LLVM REQUIRED)
      set(CMAKE_LINKER
          clang++
          CACHE STRING "C++ linker" FORCE)
      set(SYCL_BACKEND_FLAGS_LINK "")
      set(SYCL_BACKEND_FLAGS_COMPILE "")
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
        string(JOIN " " SYCL_BACKEND_FLAGS_LINK_STR ${SYCL_BACKEND_FLAGS_LINK})
        set(CMAKE_CXX_LINK_EXECUTABLE
            "${LLVM_TOOLS_BINARY_DIR}/clang++ -fsycl -fsycl-targets=${sycl_targets} <FLAGS> <CMAKE_CXX_LINK_FLAGS> <LINK_FLAGS> ${SYCL_BACKEND_FLAGS_LINK_STR} ${MUDOCK_SYCL_EXTRA_LINK_FLAGS} <OBJECTS> -o <TARGET> <LINK_LIBRARIES>"
        )

      else()
        message(WARNING "No colon found in string.")
      endif()
    endif()

    set(SYCL_TARGETS ${sycl_targets})
  else()
    message(FATAL_ERROR "Enable only one parallel framework at a time.")
  endif()
endif()

if(NOT MUDOCK_ADT_BUCKET_OVERRIDE STREQUAL "")
  add_compile_definitions(MUDOCK_ADT_BUCKET_OVERRIDE=${MUDOCK_ADT_BUCKET_OVERRIDE})
endif()
if(NOT MUDOCK_ADT_BUCKET_MULTIPLE_OVERRIDE STREQUAL "")
  add_compile_definitions(MUDOCK_ADT_BUCKET_MULTIPLE_OVERRIDE=${MUDOCK_ADT_BUCKET_MULTIPLE_OVERRIDE})
endif()
if(NOT MUDOCK_UNROLL_FACTOR STREQUAL "")
  add_compile_definitions(MUDOCK_UNROLL_FACTOR=${MUDOCK_UNROLL_FACTOR})
endif()
if(MUDOCK_DISABLE_UNROLL)
  add_compile_definitions(MUDOCK_DISABLE_UNROLL)
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
find_package(Boost CONFIG REQUIRED COMPONENTS program_options graph fiber)
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
    -Wpedantic
    -Wshadow
    -Wdouble-promotion)
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

# ##############################################################################
# ### Add each component that compose this project
# ##############################################################################
add_subdirectory(application)
add_subdirectory(mudock)

# ##############################################################################
# ### Print a small recap of the project configuration
# ##############################################################################
message(STATUS "")
message(STATUS "//===---------------------------===//")
message(STATUS "//       CONFIGURATION RECAP       //")
message(STATUS "//===---------------------------===//")
message(
  STATUS
    "CXX compiler........................: ${CMAKE_CXX_COMPILER} (CMAKE_CXX_COMPILER)"
)
message(
  STATUS
    "Build type..........................: ${CMAKE_BUILD_TYPE} (CMAKE_BUILD_TYPE)"
)
message(
  STATUS
    "Install prefix......................: ${CMAKE_INSTALL_PREFIX} (CMAKE_INSTALL_PREFIX)"
)
message(
  STATUS
    "Bucket level........................: ${MUDOCK_BUCKET_LEVEL} (MUDOCK_BUCKET_LEVEL)"
)
message(
  STATUS
    "CPU Architecture....................: ${MUDOCK_CPU_ARCHITECTURES} (MUDOCK_CPU_ARCHITECTURES)"
)
if(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64")
  message(
    STATUS
      "CPU Target..........................: ${MUDOCK_CPU_TARGET} (MUDOCK_CPU_TARGET)"
  )
endif()
if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
  message(
    STATUS
      "CPU Tune............................: ${MUDOCK_CPU_TUNE} (MUDOCK_CPU_TUNE)"
  )
endif()

message(
  STATUS
    "Fast mode ..........................: ${MUDOCK_ENABLE_FAST} (MUDOCK_ENABLE_FAST)"
)

if(MUDOCK_ENABLE_TEST)
  message(
    STATUS
      "Enable Test.........................: ${MUDOCK_ENABLE_TEST} (MUDOCK_ENABLE_TEST)"
  )
endif()

message(
  STATUS
    "IPO Optimizations...................: ${MUDOCK_ENABLE_IPO} (MUDOCK_ENABLE_IPO)"
)

message(
  STATUS
    "Likwid Profiling....................: ${MUDOCK_ENABLE_LIKWID} (MUDOCK_ENABLE_LIKWID)"
)
message(
  STATUS
    "Disable kernel unrolling............: ${MUDOCK_DISABLE_UNROLL} (MUDOCK_DISABLE_UNROLL)"
)
if(NOT MUDOCK_DISABLE_UNROLL)
  message(
    STATUS
      "Kernel unroll factor.................: ${MUDOCK_UNROLL_FACTOR} (MUDOCK_UNROLL_FACTOR)"
  )
endif()
if(NOT MUDOCK_ADT_BUCKET_OVERRIDE STREQUAL "")
  message(
    STATUS
      "ADT bucket override.................: ${MUDOCK_ADT_BUCKET_OVERRIDE} (MUDOCK_ADT_BUCKET_OVERRIDE)"
  )
endif()
if(NOT MUDOCK_ADT_BUCKET_MULTIPLE_OVERRIDE STREQUAL "")
  message(
    STATUS
      "ADT bucket multiple override........: ${MUDOCK_ADT_BUCKET_MULTIPLE_OVERRIDE} (MUDOCK_ADT_BUCKET_MULTIPLE_OVERRIDE)"
  )
endif()
if(MUDOCK_ENABLE_LIKWID)
  message(
  STATUS
    "Likwid CPP Profiling................: ${MUDOCK_LIKWID_CPP} (MUDOCK_LIKWID_CPP)"
)
  message(
  STATUS
    "Likwid CUDA Profiling...............: ${MUDOCK_LIKWID_CUDA} (MUDOCK_LIKWID_CUDA)"
)
endif()

message(
  STATUS
    "Enable Google Highway...............: ${MUDOCK_ENABLE_GH} (MUDOCK_ENABLE_GH)"
)

message(
  STATUS
    "Enable XSimd........................: ${MUDOCK_ENABLE_XSIMD} (MUDOCK_ENABLE_XSIMD)"
)

message(
  STATUS
    "Enable CUDA.........................: ${MUDOCK_ENABLE_CUDA} (MUDOCK_ENABLE_CUDA)"
)
message(
  STATUS
    "Kernel lock..........................: ${MUDOCK_ENABLE_KERNEL_LOCK} (MUDOCK_ENABLE_KERNEL_LOCK)"
)

if(MUDOCK_ENABLE_CUDA)
  # TODO use instead the MUDOCK_GPU_ARCHITECTURES
  message(
    STATUS
      "Target CUDA architecture............: ${CMAKE_CUDA_ARCHITECTURES} (CMAKE_CUDA_ARCHITECTURES)"
  )
endif()

message(
  STATUS
    "Enable HIP..........................: ${MUDOCK_ENABLE_HIP} (MUDOCK_ENABLE_HIP)"
)

if(MUDOCK_ENABLE_HIP)
  message(
    STATUS
      "Target HIP platform.................: ${HIP_PLATFORM} (HIP_PLATFORM)")
  message(
    STATUS
      "Target HIP targets..................: ${CMAKE_HIP_ARCHITECTURES} (CMAKE_HIP_ARCHITECTURES)"
  )
endif()

message(
  STATUS
    "Enable SYCL.........................: ${MUDOCK_ENABLE_SYCL} (MUDOCK_ENABLE_SYCL)"
)

if(MUDOCK_ENABLE_SYCL)
  message(
    STATUS
      "Target SYCL architecture............: ${SYCL_TARGETS} (SYCL_TARGETS)")
  message(
    STATUS
      "SYCL Compiler.......................: ${SYCL_COMPILER} (SYCL_COMPILER)")
  message(
    STATUS
      "SYCL Work-group size................: ${MUDOCK_SYCL_WG_SIZE} (MUDOCK_SYCL_WG_SIZE)"
  )
  if(${SYCL_COMPILER} STREQUAL "oneAPI")
    message(
      STATUS
        "oneAPI Compiler Path................: ${LLVM_TOOLS_BINARY_DIR} (LLVM_TOOLS_BINARY_DIR)"
    )
    if(NOT SYCL_BACKEND_FLAGS_LINK STREQUAL "")
      string(JOIN " " SYCL_BACKEND_FLAGS_PRINT ${SYCL_BACKEND_FLAGS_LINK})
      message(
        STATUS
          "oneAPI Backend Flags................: ${SYCL_BACKEND_FLAGS_PRINT} (SYCL_BACKEND_FLAGS_LINK)"
      )
    endif()
    if(NOT MUDOCK_SYCL_EXTRA_LINK_FLAGS STREQUAL "")
      message(
        STATUS
          "oneAPI Extra Link Flags.............: ${MUDOCK_SYCL_EXTRA_LINK_FLAGS} (MUDOCK_SYCL_EXTRA_LINK_FLAGS)"
      )
    endif()
  endif()
endif()

message(
  STATUS
    "Enable OpenMP.......................: ${MUDOCK_ENABLE_OMP} (MUDOCK_ENABLE_OMP)"
)
message(
  STATUS
    "Enable OpenMP Offload...............: ${MUDOCK_ENABLE_OMP_OFFLOAD} (MUDOCK_ENABLE_OMP_OFFLOAD)"
)

if(MUDOCK_ENABLE_OMP_OFFLOAD)
  message(
    STATUS "Target OpenMP architecture..........: ${OMP_TARGET} (OMP_TARGET)")
endif()

message(
  STATUS
    "Enable MPI..........................: ${MUDOCK_ENABLE_MPI} (MUDOCK_ENABLE_MPI)"
)

message(STATUS "//===---------------------------===//")
message(STATUS "")
