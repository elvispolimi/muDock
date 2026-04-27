# ##############################################################################
# ### Declare the application source files
# ##############################################################################

# headers
set(header_path "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(header_prefix "mudock")
list(
  APPEND
  header_files
  "${header_path}/${header_prefix}/chem/assign_autodock_babel_types.hpp"
  "${header_path}/${header_prefix}/chem/autodock_babel_types.hpp"
  "${header_path}/${header_prefix}/chem/autodock_grid_types.hpp"
  "${header_path}/${header_prefix}/chem/autodock_layer.hpp"
  "${header_path}/${header_prefix}/chem/autodock_ligand.hpp"
  "${header_path}/${header_prefix}/chem/autodock_parameters.hpp"
  "${header_path}/${header_prefix}/chem/autodock_protein.hpp"
  "${header_path}/${header_prefix}/chem/autodock_types.hpp"
  "${header_path}/${header_prefix}/chem/bond_types.hpp"
  "${header_path}/${header_prefix}/chem/cleanup.hpp"
  "${header_path}/${header_prefix}/chem/elements.hpp"
  "${header_path}/${header_prefix}/chem/geom_ligand.hpp"
  "${header_path}/${header_prefix}/chem/grid_const.hpp"
  "${header_path}/${header_prefix}/chem/mehler_solmajer.hpp"
  "${header_path}/${header_prefix}/chem/molecule_layer.hpp"
  "${header_path}/${header_prefix}/compute/safe_stack.hpp"
  "${header_path}/${header_prefix}/compute/parse_ids.hpp"
  "${header_path}/${header_prefix}/compute/reorder_buffer.hpp"
  "${header_path}/${header_prefix}/compute/threadpool.hpp"
  "${header_path}/${header_prefix}/compute/worker.hpp"
  "${header_path}/${header_prefix}/compute/queue.hpp"
  "${header_path}/${header_prefix}/compute/object.hpp"
  "${header_path}/${header_prefix}/compute/scratchpad.hpp"
  "${header_path}/${header_prefix}/compute/manager.hpp"
  "${header_path}/${header_prefix}/compute/buffer.hpp"
  "${header_path}/${header_prefix}/compute/stage.hpp"
  "${header_path}/${header_prefix}/compute/scoring.hpp"
  "${header_path}/${header_prefix}/compute/adt_score.hpp"
  "${header_path}/${header_prefix}/compute/pipeline.hpp"
  "${header_path}/${header_prefix}/compute/docking.hpp"
  "${header_path}/${header_prefix}/compute/genetic.hpp"
  "${header_path}/${header_prefix}/compute/transform.hpp"
  "${header_path}/${header_prefix}/compute/geometric_transform.hpp"
  "${header_path}/${header_prefix}/cpp_implementation/chromosome.hpp"
  "${header_path}/${header_prefix}/cpp_implementation/vectorization.hpp"
  "${header_path}/${header_prefix}/cpp_implementation/center_of_mass.hpp"
  "${header_path}/${header_prefix}/cpp_implementation/non_bonds.hpp"
  "${header_path}/${header_prefix}/cpp_implementation/buffer_cpp.hpp"
  "${header_path}/${header_prefix}/cpp_implementation/queue_cpp.hpp"
  "${header_path}/${header_prefix}/cpp_implementation/adt_score_cpp.hpp"
  "${header_path}/${header_prefix}/cpp_implementation/geom_transform_cpp.hpp"
  "${header_path}/${header_prefix}/cpp_implementation/genetic_cpp.hpp"
  "${header_path}/${header_prefix}/cpp_implementation/cpp_implementation.hpp"
  "${header_path}/${header_prefix}/cpp_implementation/mutate_cpp.hpp"
  "${header_path}/${header_prefix}/cpp_implementation/mutate.hpp"
  "${header_path}/${header_prefix}/format/concept.hpp"
  "${header_path}/${header_prefix}/format/dot.hpp"
  "${header_path}/${header_prefix}/format/mol2.hpp"
  "${header_path}/${header_prefix}/format/pdb.hpp"
  "${header_path}/${header_prefix}/format/supported_format.hpp"
  "${header_path}/${header_prefix}/format/adt_mol2.hpp"
  "${header_path}/${header_prefix}/format/writer.hpp"
  "${header_path}/${header_prefix}/format/pdbqt.hpp"
  "${header_path}/${header_prefix}/format/reader.hpp"
  "${header_path}/${header_prefix}/format/ob_wrapper.hpp"
  "${header_path}/${header_prefix}/grid/pi.hpp"
  "${header_path}/${header_prefix}/grid/mdindex.hpp"
  "${header_path}/${header_prefix}/grid/mdspan.hpp"
  "${header_path}/${header_prefix}/grid/space_grid.hpp"
  "${header_path}/${header_prefix}/grid/point3D.hpp"
  "${header_path}/${header_prefix}/molecule/bond.hpp"
  "${header_path}/${header_prefix}/molecule/constraints.hpp"
  "${header_path}/${header_prefix}/molecule/containers.hpp"
  "${header_path}/${header_prefix}/molecule/fragments.hpp"
  "${header_path}/${header_prefix}/molecule/graph.hpp"
  "${header_path}/${header_prefix}/molecule/properties.hpp"
  "${header_path}/${header_prefix}/molecule/property_table.hpp"
  "${header_path}/${header_prefix}/tbb_implementation/stream_filter.hpp"
  "${header_path}/${header_prefix}/tbb_implementation/parser_filter.hpp"
  "${header_path}/${header_prefix}/tbb_implementation/runtime_services.hpp"
  "${header_path}/${header_prefix}/tbb_implementation/tbb_pipeline.hpp"
  "${header_path}/${header_prefix}/batch.hpp"
  "${header_path}/${header_prefix}/chem.hpp"
  "${header_path}/${header_prefix}/compute.hpp"
  "${header_path}/${header_prefix}/format.hpp"
  "${header_path}/${header_prefix}/grid.hpp"
  "${header_path}/${header_prefix}/knobs.hpp"
  "${header_path}/${header_prefix}/log.hpp"
  "${header_path}/${header_prefix}/molecule.hpp"
  "${header_path}/${header_prefix}/mudock.hpp"
  "${header_path}/${header_prefix}/splitter.hpp"
  "${header_path}/${header_prefix}/type_alias.hpp"
  "${header_path}/${header_prefix}/utils.hpp")

# sources
set(source_path "${CMAKE_CURRENT_SOURCE_DIR}/src")
list(
  APPEND
  source_files
  "${source_path}/chem/assign_autodock_babel_types.cpp"
  "${source_path}/chem/assign_autodock_types.cpp"
  "${source_path}/chem/autodock_babel_types.cpp"
  "${source_path}/chem/autodock_protein.cpp"
  "${source_path}/chem/autodock_ligand.cpp"
  "${source_path}/chem/autodock_types.cpp"
  "${source_path}/chem/bond_types.cpp"
  "${source_path}/chem/autodock_grid_types.cpp"
  "${source_path}/chem/elements.cpp"
  "${source_path}/compute/parse_ids.cpp"
  "${source_path}/format/dot.cpp"
  "${source_path}/format/mol2.cpp"
  "${source_path}/format/adt_mol2.cpp"
  "${source_path}/format/ob_wrapper.cpp"
  "${source_path}/format/pdb.cpp"
  "${source_path}/format/pdbqt.cpp"
  "${source_path}/format/supported_format.cpp"
  "${source_path}/format/reader.cpp"
  "${source_path}/format/writer.cpp"
  "${source_path}/molecule/fragments.cpp"
  "${source_path}/molecule/graph.cpp"
  "${source_path}/molecule/property_table.cpp"
  "${source_path}/log.cpp"
  "${source_path}/cpp_implementation/mutate_cpp.cpp"
  "${source_path}/cpp_implementation/geom_transform_cpp.cpp"
  "${source_path}/cpp_implementation/adt_score_cpp.cpp"
  "${source_path}/cpp_implementation/genetic_cpp.cpp"
  "${source_path}/tbb_implementation/stream_filter.cpp"
  "${source_path}/tbb_implementation/parser_filter.cpp"
  "${source_path}/tbb_implementation/runtime_services.cpp"
  "${source_path}/tbb_implementation/tbb_pipeline.cpp"
)

if(MUDOCK_ENABLE_MPI)
  list(
    APPEND
    header_files
    "${header_path}/${header_prefix}/mpi_implementation/byte_range.hpp"
    "${header_path}/${header_prefix}/mpi_implementation/stream_partition.hpp"
    "${header_path}/${header_prefix}/mpi_implementation/distributed_ranges.hpp"
    "${header_path}/${header_prefix}/mpi_implementation/utilities.hpp"
    )
endif()    

if(MUDOCK_ENABLE_GH)
  list(APPEND source_files
       "${source_path}/gh_implementation/geom_transform_gh.cpp"
       "${source_path}/gh_implementation/mutate_gh.cpp"
       "${source_path}/gh_implementation/adt_score_gh.cpp")
endif()

if(MUDOCK_ENABLE_XSIMD)
  list(APPEND source_files
       "${source_path}/xsimd_implementation/geom_transform_xsimd.cpp"
       "${source_path}/xsimd_implementation/mutate_xsimd.cpp"
       "${source_path}/xsimd_implementation/adt_score_xsimd.cpp")
endif()

# declare the CUDA sources
if(MUDOCK_ENABLE_CUDA)
  list(
    APPEND
    header_files
    "${header_path}/${header_prefix}/cuda_implementation/queue_cuda.cuh"
    "${header_path}/${header_prefix}/cuda_implementation/adt_score_cuda.cuh"
    "${header_path}/${header_prefix}/cuda_implementation/genetic_cuda.cuh"
    "${header_path}/${header_prefix}/cuda_implementation/cuda_random.cuh"
    "${header_path}/${header_prefix}/cuda_implementation/cuda_utils.cuh"
    "${header_path}/${header_prefix}/cuda_implementation/mutate.cuh"
    "${header_path}/${header_prefix}/cuda_implementation/cuda_texture.cuh"
    "${header_path}/${header_prefix}/cuda_implementation/cuda_implementation.hpp"
    )
  list(
    APPEND
    source_files
    "${source_path}/cuda_implementation/queue_cuda.cu"
    "${source_path}/cuda_implementation/adt_score_cuda.cu"
    "${source_path}/cuda_implementation/genetic_cuda.cu"
    "${source_path}/cuda_implementation/cuda_random.cu"
    "${source_path}/cuda_implementation/geom_transform_cuda.cu"
  )
endif()

# declare the HIP sources
if(MUDOCK_ENABLE_HIP)
  list(
    APPEND
    header_files
    "${header_path}/${header_prefix}/hip_implementation/queue_hip.hpp"
    "${header_path}/${header_prefix}/hip_implementation/adt_score_hip.hpp"
    "${header_path}/${header_prefix}/hip_implementation/genetic_hip.hpp"
    "${header_path}/${header_prefix}/hip_implementation/hip_random.hpp"
    "${header_path}/${header_prefix}/hip_implementation/hip_texture.hpp"
    "${header_path}/${header_prefix}/hip_implementation/hip_utils.hpp"
    "${header_path}/${header_prefix}/hip_implementation/mutate.hpp"
    "${header_path}/${header_prefix}/hip_implementation/hip_implementation.hpp"
    )
  set(
    hip_source_files
    "${source_path}/hip_implementation/queue_hip.cpp"
    "${source_path}/hip_implementation/adt_score_hip.cpp"
    "${source_path}/hip_implementation/genetic_hip.cpp"
    "${source_path}/hip_implementation/hip_random.cpp"
    "${source_path}/hip_implementation/geom_transform_hip.cpp"
  )
  list(APPEND source_files ${hip_source_files})
  set_source_files_properties(${hip_source_files} PROPERTIES LANGUAGE HIP)
endif()

# declare the SYCL sources
if(MUDOCK_ENABLE_SYCL)
  list(
    APPEND
    header_files_sycl
    "${header_path}/${header_prefix}/sycl_implementation/queue_sycl.hpp"
    "${header_path}/${header_prefix}/sycl_implementation/invoke_kernel_sycl.hpp"
    "${header_path}/${header_prefix}/sycl_implementation/sycl_texture.hpp"
    "${header_path}/${header_prefix}/sycl_implementation/adt_score_sycl.hpp"
    "${header_path}/${header_prefix}/sycl_implementation/queue_sycl_impl.hpp"
    "${header_path}/${header_prefix}/sycl_implementation/sycl_kernel_lock.hpp"
    "${header_path}/${header_prefix}/sycl_implementation/invoke_kernel_sycl.hpp"
    "${header_path}/${header_prefix}/sycl_implementation/mutate.hpp"
    "${header_path}/${header_prefix}/sycl_implementation/sycl_random.hpp"
    "${header_path}/${header_prefix}/sycl_implementation/geom_transform_sycl.hpp"
    "${header_path}/${header_prefix}/sycl_implementation/genetic_sycl.hpp"
   )
  list(APPEND header_files "${header_files_sycl}")
  set(sycl_source_files
      "${source_path}/sycl_implementation/queue_sycl.cpp"
      "${source_path}/sycl_implementation/adt_score_sycl.cpp"
      "${source_path}/sycl_implementation/geom_transform_sycl.cpp"
      "${source_path}/sycl_implementation/sycl_random.cpp"
      "${source_path}/sycl_implementation/genetic_sycl.cpp"
      "${source_path}/sycl_implementation/sycl_kernel_lock.cpp"
  )
endif()

# declare the OpenMP sources
if(MUDOCK_ENABLE_OMP_OFFLOAD)
  list(
    APPEND
    header_files
    "${header_path}/${header_prefix}/omp_implementation/omp_manager.hpp"
    "${header_path}/${header_prefix}/omp_implementation/calc_energy.hpp"
    "${header_path}/${header_prefix}/omp_implementation/mutate.hpp"
    "${header_path}/${header_prefix}/omp_implementation/geometric_transformations.hpp"
    "${header_path}/${header_prefix}/omp_implementation/omp_random.hpp"
    "${header_path}/${header_prefix}/omp_implementation/omp_object.hpp"
    "${header_path}/${header_prefix}/omp_implementation/omp_wrapper.hpp"
    "${header_path}/${header_prefix}/omp_implementation/omp_worker.hpp"
    "${header_path}/${header_prefix}/omp_implementation/omp_batch_sizer.hpp"
    "${header_path}/${header_prefix}/omp_implementation/virtual_screen.hpp"
    "${header_path}/${header_prefix}/omp_implementation/evaluate_fitness.hpp")
  list(
    APPEND
    source_files
    "${source_path}/omp_implementation/omp_manager.cpp"
    "${source_path}/omp_implementation/omp_random.cpp"
    "${source_path}/omp_implementation/omp_worker.cpp"
    "${source_path}/omp_implementation/omp_batch_sizer.cpp"
    "${source_path}/omp_implementation/virtual_screen.cpp"
    "${source_path}/omp_implementation/evaluate_fitness.cpp")
endif()

# ##############################################################################
# ### Declare the compile step
# ##############################################################################

# declare the mudock library
find_package(TBB REQUIRED)
add_library(libmudock STATIC ${source_files} ${header_files})
add_library(mudock::libmudock ALIAS libmudock)
set_target_properties(libmudock PROPERTIES PUBLIC_HEADER "${headers}")
target_include_directories(libmudock PUBLIC "${header_path}")
target_include_directories(libmudock SYSTEM PUBLIC ${OpenBabel3_INCLUDE_DIRS})
set_target_properties(libmudock PROPERTIES PREFIX "" OUTPUT_NAME "libmudock")
set_target_properties(
  libmudock
  PROPERTIES CXX_STANDARD 20
             CXX_STANDARD_REQUIRED ON
             CXX_EXTENSIONS OFF
             POSITION_INDEPENDENT_CODE ON)
target_link_libraries(libmudock PUBLIC Boost::graph ${OpenBabel3_LIBRARIES} TBB::tbb)
if(NOT MUDOCK_ENABLE_OMP)
  target_compile_options(libmudock PRIVATE -Wno-unknown-pragmas)
endif()

target_compile_options(
  libmudock
  PUBLIC $<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANGUAGE:CXX>>:${release_flags}>)

# Ensure cache-driven compile definitions are applied to the target so custom
# SYCL compile commands can see them via COMPILE_DEFINITIONS.
if(MUDOCK_ENABLE_MPI)
  target_compile_definitions(libmudock PUBLIC MUDOCK_USE_MPI)
  target_link_libraries(libmudock PRIVATE MPI::MPI_CXX)
endif()
if(MUDOCK_ENABLE_KERNEL_LOCK)
  target_compile_definitions(libmudock PUBLIC MUDOCK_KERNEL_LOCK)
endif()
if(NOT MUDOCK_STAGE_BUCKET_OVERRIDE STREQUAL "")
  target_compile_definitions(libmudock PUBLIC MUDOCK_STAGE_BUCKET_OVERRIDE=${MUDOCK_STAGE_BUCKET_OVERRIDE})
endif()
if(NOT MUDOCK_STAGE_BUCKET_MULTIPLE_OVERRIDE STREQUAL "")
  target_compile_definitions(libmudock PUBLIC MUDOCK_STAGE_BUCKET_MULTIPLE_OVERRIDE=${MUDOCK_STAGE_BUCKET_MULTIPLE_OVERRIDE})
endif()
if(MUDOCK_STAGE_BUCKET_POLICY STREQUAL "MAX_UTILIZATION")
  target_compile_definitions(libmudock PUBLIC MUDOCK_STAGE_BUCKET_POLICY_MAX_UTILIZATION)
elseif(MUDOCK_STAGE_BUCKET_POLICY STREQUAL "SM_ALIGNED")
  target_compile_definitions(libmudock PUBLIC MUDOCK_STAGE_BUCKET_POLICY_SM_ALIGNED)
endif()
if(MUDOCK_GENETIC_BUCKET_COMBINE STREQUAL "SCORE_ONLY")
  target_compile_definitions(libmudock PUBLIC MUDOCK_GENETIC_BUCKET_COMBINE_SCORE_ONLY)
elseif(MUDOCK_GENETIC_BUCKET_COMBINE STREQUAL "GEOM_ONLY")
  target_compile_definitions(libmudock PUBLIC MUDOCK_GENETIC_BUCKET_COMBINE_GEOM_ONLY)
elseif(MUDOCK_GENETIC_BUCKET_COMBINE STREQUAL "LCM")
  target_compile_definitions(libmudock PUBLIC MUDOCK_GENETIC_BUCKET_COMBINE_LCM)
elseif(NOT MUDOCK_GENETIC_BUCKET_COMBINE STREQUAL "MIN")
  message(FATAL_ERROR "MUDOCK_GENETIC_BUCKET_COMBINE must be one of: MIN, SCORE_ONLY, GEOM_ONLY, LCM")
endif()
if(NOT MUDOCK_UNROLL_FACTOR STREQUAL "")
  target_compile_definitions(libmudock PUBLIC MUDOCK_UNROLL_FACTOR=${MUDOCK_UNROLL_FACTOR})
endif()
if(MUDOCK_DISABLE_UNROLL)
  target_compile_definitions(libmudock PUBLIC MUDOCK_DISABLE_UNROLL)
endif()
if(MUDOCK_ENABLE_STAGE_BUCKET_TRACE)
  target_compile_definitions(libmudock PUBLIC MUDOCK_ENABLE_STAGE_BUCKET_TRACE)
endif()

if(MUDOCK_ATOM_CLUSTER_LEVEL STREQUAL "EXTREME")
  target_compile_definitions(libmudock PUBLIC MUDOCK_ATOM_CLUSTER_LEVEL_EXTREME)
elseif(MUDOCK_ATOM_CLUSTER_LEVEL STREQUAL "LARGE")
  target_compile_definitions(libmudock PUBLIC MUDOCK_ATOM_CLUSTER_LEVEL_LARGE)
elseif(MUDOCK_ATOM_CLUSTER_LEVEL STREQUAL "MEDIUM")
  target_compile_definitions(libmudock PUBLIC MUDOCK_ATOM_CLUSTER_LEVEL_MEDIUM)
endif()

# Profiling
if(MUDOCK_ENABLE_LIKWID)
  target_link_libraries(libmudock PUBLIC likwid::likwid)
  target_include_directories(libmudock PUBLIC ${LIKWID_INCLUDE_DIR})
  if(MUDOCK_ENABLE_CUDA)
    target_compile_definitions(libmudock PUBLIC LIKWID_NVMON)
    if(MUDOCK_LIKWID_CUDA)
      target_compile_definitions(libmudock PUBLIC MUDOCK_LIKWID_CUDA)
    endif()
  endif()
  if(MUDOCK_LIKWID_CPP)
    target_compile_definitions(libmudock PUBLIC MUDOCK_LIKWID_CPP)
    target_compile_definitions(libmudock PUBLIC LIKWID_PERFMON)
  endif()
endif()

# Test
if(MUDOCK_ENABLE_TEST)
  target_compile_definitions(libmudock PUBLIC MUDOCK_TEST)
endif()

if(MUDOCK_ENABLE_GH)
  target_compile_definitions(libmudock PUBLIC MUDOCK_USE_GH)
  target_link_libraries(libmudock PRIVATE hwy::hwy)
endif()

if(MUDOCK_ENABLE_XSIMD)
  target_compile_definitions(libmudock PUBLIC MUDOCK_USE_XSIMD)
  target_include_directories(libmudock PUBLIC "${XSIMD_INCLUDE_DIR}")
endif()

if(MUDOCK_ENABLE_FAST)
  if(CMAKE_C_COMPILER_ID STREQUAL "GNU" OR CMAKE_C_COMPILER_ID STREQUAL "Clang")
    target_compile_options(
      libmudock PUBLIC $<$<COMPILE_LANGUAGE:CXX>:-fopenmp-simd -ffast-math>)
  else()
    message(
      STATUS
        "MuDock does not optimize SIMD auto-vectorization with ${CMAKE_C_COMPILER_ID}"
    )
  endif()

endif()

if(MUDOCK_ENABLE_HIP)
  target_compile_definitions(libmudock PUBLIC MUDOCK_USE_HIP)
  target_include_directories(
    libmudock SYSTEM PRIVATE $<$<COMPILE_LANGUAGE:HIP>:${hip_INCLUDE_DIRS}>)
  target_link_libraries(libmudock PUBLIC hip::hiprand hip::host)
  target_compile_options(
    libmudock
    PUBLIC
      $<$<AND:$<COMPILE_LANGUAGE:HIP>,$<HIP_COMPILER_ID:NVIDIA>>:--expt-relaxed-constexpr
      -rdc=true -Xptxas -v --resource-usage>)
  target_compile_options(
    libmudock
    PUBLIC $<$<AND:$<COMPILE_LANGUAGE:HIP>,$<HIP_COMPILER_ID:Clang>>:-fgpu-rdc -Rpass-analysis=kernel-resource-usage>)
  target_link_options(
    libmudock PUBLIC
    $<$<AND:$<COMPILE_LANGUAGE:HIP>,$<HIP_COMPILER_ID:Clang>>:-fgpu-rdc
    --hip-link -Rpass-analysis=kernel-resource-usage>)
  target_compile_options(
    libmudock
    PUBLIC
      $<$<AND:$<COMPILE_LANGUAGE:HIP>,$<CONFIG:RELEASE>,$<HIP_COMPILER_ID:NVIDIA>>:
      -use_fast_math>)
  # FIXME
  if(HIP_PLATFORM STREQUAL "nvidia")
    # FIX ME some platform requires -Wl,-rpath to link the correct CUDA rand lib
    target_link_libraries(libmudock PUBLIC CUDA::curand)
  elseif(HIP_PLATFORM STREQUAL "amd")
    if(DEFINED ROCRAND_INCLUDE_DIRS)
      target_include_directories(libmudock SYSTEM PRIVATE ${ROCRAND_INCLUDE_DIRS})
    endif()
    target_link_libraries(libmudock PUBLIC rocrand)
  endif()
endif()

# check if we need to activate the compilation flag that signals CUDA
if(MUDOCK_ENABLE_CUDA)
  target_compile_definitions(libmudock PUBLIC MUDOCK_USE_CUDA)
  if(CMAKE_BUILD_TYPE STREQUAL "Release")
    # CUDA_SEPARABLE_COMPILATION ON required for constant memory to be shared
    # across translation unit
    set_target_properties(
      libmudock
      PROPERTIES CUDA_STANDARD 20
                 CUDA_STANDARD_REQUIRED ON
                 CUDA_EXTENSIONS OFF
                 CUDA_SEPARABLE_COMPILATION ON)
  else()
    set_target_properties(
      libmudock
      PROPERTIES CUDA_STANDARD 20
                 CUDA_STANDARD_REQUIRED ON
                 CUDA_EXTENSIONS OFF
                 CUDA_SEPARABLE_COMPILATION ON)
  endif()
  # set(CMAKE_CUDA_ARCHITECTURES OFF)
  set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr")
  set(CMAKE_CUDA_FLAGS_RELEASE
      "${CMAKE_CUDA_FLAGS_RELEASE} -use_fast_math -Xptxas -v")
  set(CMAKE_CUDA_FLAGS_RELWITHDEBINFO
      "${CMAKE_CUDA_FLAGS_RELWITHDEBINFO} -G -Xptxas -v")
  set(CMAKE_CUDA_FLAGS_DEBUG "${CMAKE_CUDA_FLAGS_DEBUG} -G -Xptxas -v")
endif()

# check if we need to activate the compilation flag that signals OpenMP
if(MUDOCK_ENABLE_OMP_OFFLOAD)
  target_compile_definitions(libmudock PUBLIC MUDOCK_USE_OMP)
  target_link_libraries(libmudock PUBLIC OpenMP::OpenMP_CXX)

  # Enable offloading to AMD GPUs using Clang
  if(NOT OMP_TARGET EQUAL "cpu")
    if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
      if(GPU_PLATFORM STREQUAL "nvidia")
        target_compile_options(libmudock PRIVATE "-fopenmp")
        target_compile_options(libmudock PRIVATE "-fopenmp-targets=nvptx64")
        # target_compile_options(libmudock PRIVATE "-Xopenmp-target=nvptx64")
      elseif(GPU_PLATFORM STREQUAL "amd")
        target_compile_options(libmudock PRIVATE "-fopenmp")
        target_compile_options(libmudock
                               PRIVATE "-fopenmp-targets=amdgcn-amd-amdhsa")
        # target_compile_options(libmudock PRIVATE
        # "-Xopenmp-target=amdgcn-amd-amdhsa")
      endif()
      target_compile_options(libmudock
                             PRIVATE "--offload-arch=${GPU_ARCHITECTURE}")
    elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
      # TODO test with GCC, check also if GCC supports this
      message(FATAL_ERROR "muDock does not support GCC and OpenMP offloading")
      # if(GPU_PLATFORM STREQUAL "cuda") target_compile_options(libmudock
      # PRIVATE "-fopenmp") target_compile_options(libmudock PRIVATE
      # "-foffload=nvptx-none") target_compile_options(libmudock PRIVATE
      # "-foffload-options=nvptx-none=-march=${GPU_ARCHITECTURE}")
      # elseif(GPU_PLATFORM STREQUAL "amd") target_compile_options(libmudock
      # PRIVATE "-fopenmp") target_compile_options(libmudock PRIVATE
      # "-foffload=amdgcn-amdhsa") target_compile_options(libmudock PRIVATE
      # "-foffload-options=amdgcn-amdhsa=-march=${GPU_ARCHITECTURE}") endif()
    else()
      message(FATAL_ERROR "Not support OpenMP compiler")
    endif()
  endif()
endif()

# check if we need to activate the compilation flag that signals SYCL
if(MUDOCK_ENABLE_SYCL)
  target_compile_definitions(libmudock PUBLIC MUDOCK_USE_SYCL)
  if(NOT MUDOCK_SYCL_WG_SIZE STREQUAL "")
    target_compile_definitions(libmudock PUBLIC MUDOCK_SYCL_WG_SIZE=${MUDOCK_SYCL_WG_SIZE})
  endif()
  set(sycl_generated_source)
  set(sycl_objects"")
  add_sycl_files("${sycl_source_files}" ${header_path} "${header_files}"
                 sycl_objects)
  target_sources(libmudock PRIVATE "${sycl_objects}")
  add_dependencies(libmudock sycl_targets)
  target_link_directories(libmudock PUBLIC "${LLVM_LIBRARY_DIRS}")
  target_link_libraries(libmudock PRIVATE sycl)
endif()

if(MUDOCK_ENABLE_OMP)
  target_link_libraries(libmudock PRIVATE OpenMP::OpenMP_CXX)
  target_compile_options(libmudock PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-fopenmp>)
endif()

# ##############################################################################
# ### Declare the install step
# ##############################################################################

# at the moment wedon't deploy the library
