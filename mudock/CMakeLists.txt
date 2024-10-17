# ########################################################################
# ### Declare the application source files
# ########################################################################

# headers
set(header_path "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(header_prefix "mudock")
list(APPEND header_files
  "${header_path}/${header_prefix}/chem/apply_autodock_force_field.hpp"
  "${header_path}/${header_prefix}/chem/assign_autodock_babel_types.hpp"
  "${header_path}/${header_prefix}/chem/assign_autodock_types.hpp"
  "${header_path}/${header_prefix}/chem/autodock_babel_types.hpp"
  "${header_path}/${header_prefix}/chem/autodock_types.hpp"
  "${header_path}/${header_prefix}/chem/autodock_parameters.hpp"
  "${header_path}/${header_prefix}/chem/bond_types.hpp"
  "${header_path}/${header_prefix}/chem/elements.hpp"
  "${header_path}/${header_prefix}/compute/device_names.hpp"
  "${header_path}/${header_prefix}/compute/parse_ids.hpp"
  "${header_path}/${header_prefix}/compute/reorder_buffer.hpp"
  "${header_path}/${header_prefix}/compute/safe_stack.hpp"
  "${header_path}/${header_prefix}/compute/threadpool.hpp"
  "${header_path}/${header_prefix}/cpp_implementation/cpp_manager.hpp"
  "${header_path}/${header_prefix}/cpp_implementation/cpp_worker.hpp"
  "${header_path}/${header_prefix}/cpp_implementation/geometric_transformations.hpp"
  "${header_path}/${header_prefix}/cpp_implementation/virtual_screen.hpp"
  "${header_path}/${header_prefix}/cpp_implementation/trilinear_interpolation.hpp"
  "${header_path}/${header_prefix}/cpp_implementation/chromosome.hpp"
  "${header_path}/${header_prefix}/cpp_implementation/mutate.hpp"
  "${header_path}/${header_prefix}/cpp_implementation/center_of_mass.hpp"
  "${header_path}/${header_prefix}/cpp_implementation/weed_bonds.hpp"
  "${header_path}/${header_prefix}/cpp_implementation/calc_energy_cpp.hpp"
  "${header_path}/${header_prefix}/format/concept.hpp"
  "${header_path}/${header_prefix}/format/dot.hpp"
  "${header_path}/${header_prefix}/format/mol2.hpp"
  "${header_path}/${header_prefix}/format/pdb.hpp"
  "${header_path}/${header_prefix}/format/rdkit_wrapper.hpp"
  "${header_path}/${header_prefix}/grid/mdindex.hpp"
  "${header_path}/${header_prefix}/grid/grid_map.hpp"
  "${header_path}/${header_prefix}/grid/pi.hpp"
  "${header_path}/${header_prefix}/grid/point3D.hpp"
  "${header_path}/${header_prefix}/molecule/bond.hpp"
  "${header_path}/${header_prefix}/molecule/constraints.hpp"
  "${header_path}/${header_prefix}/molecule/containers.hpp"
  "${header_path}/${header_prefix}/molecule/fragments.hpp"
  "${header_path}/${header_prefix}/molecule/graph.hpp"
  "${header_path}/${header_prefix}/molecule/properties.hpp"
  "${header_path}/${header_prefix}/molecule/property_table.hpp"
  "${header_path}/${header_prefix}/batch.hpp"
  "${header_path}/${header_prefix}/chem.hpp"
  "${header_path}/${header_prefix}/compute.hpp"
  "${header_path}/${header_prefix}/cpp_implementation.hpp"
  "${header_path}/${header_prefix}/cuda_implementation.hpp"
  "${header_path}/${header_prefix}/sycl_implementation.hpp"
  "${header_path}/${header_prefix}/format.hpp"
  "${header_path}/${header_prefix}/grid.hpp"
  "${header_path}/${header_prefix}/knobs.hpp"
  "${header_path}/${header_prefix}/log.hpp"
  "${header_path}/${header_prefix}/molecule.hpp"
  "${header_path}/${header_prefix}/mudock.hpp"
  "${header_path}/${header_prefix}/splitter.hpp"
  "${header_path}/${header_prefix}/type_alias.hpp"
  "${header_path}/${header_prefix}/utils.hpp"
)

# sources
set(source_path "${CMAKE_CURRENT_SOURCE_DIR}/src")
list(APPEND source_files
  "${source_path}/chem/assign_autodock_babel_types.cpp"
  "${source_path}/chem/assign_autodock_types.cpp"
  "${source_path}/chem/autodock_babel_types.cpp"
  "${source_path}/chem/autodock_types.cpp"
  "${source_path}/chem/bond_types.cpp"
  "${source_path}/chem/elements.cpp"
  "${source_path}/compute/parse_ids.cpp"
  "${source_path}/compute/reorder_buffer.cpp"
  "${source_path}/cpp_implementation/cpp_manager.cpp"
  "${source_path}/cpp_implementation/cpp_worker.cpp"
  "${source_path}/cpp_implementation/mutate.cpp"
  "${source_path}/cpp_implementation/geometric_transformations.cpp"
  "${source_path}/cpp_implementation/virtual_screen.cpp"
  "${source_path}/cpp_implementation/weed_bonds.cpp"
  "${source_path}/cpp_implementation/calc_energy_cpp.cpp"
  "${source_path}/format/dot.cpp"
  "${source_path}/format/mol2.cpp"
  "${source_path}/format/rdkit_wrapper.cpp"
  "${source_path}/grid/grid_atom_map.cpp"
  "${source_path}/grid/grid_electrostatic_map.cpp"
  "${source_path}/grid/grid_desolvation_map.cpp"
  "${source_path}/molecule/fragments.cpp"
  "${source_path}/molecule/graph.cpp"
  "${source_path}/molecule/property_table.cpp"
  "${source_path}/log.cpp"
)

# declare the CUDA sources
if(MUDOCK_ENABLE_CUDA)
  list(APPEND header_files
    "${header_path}/${header_prefix}/cuda_implementation/calc_energy.cuh"
    "${header_path}/${header_prefix}/cuda_implementation/cuda_batch_sizer.cuh"
    "${header_path}/${header_prefix}/cuda_implementation/cuda_check_error_macro.cuh"
    "${header_path}/${header_prefix}/cuda_implementation/cuda_manager.hpp"
    "${header_path}/${header_prefix}/cuda_implementation/cuda_random.cuh"
    "${header_path}/${header_prefix}/cuda_implementation/cuda_object.cuh"
    "${header_path}/${header_prefix}/cuda_implementation/cuda_worker.hpp"
    "${header_path}/${header_prefix}/cuda_implementation/cuda_wrapper.cuh"
    "${header_path}/${header_prefix}/cuda_implementation/evaluate_fitness.cuh"
    "${header_path}/${header_prefix}/cuda_implementation/geometric_transformations.cuh"
    "${header_path}/${header_prefix}/cuda_implementation/map_textures.cuh"
    "${header_path}/${header_prefix}/cuda_implementation/mutate.cuh"
    "${header_path}/${header_prefix}/cuda_implementation/virtual_screen.cuh"
  )
  list(APPEND source_files
    "${source_path}/cuda_implementation/calc_energy.cu"
    "${source_path}/cuda_implementation/cuda_object.cu"
    "${source_path}/cuda_implementation/cuda_batch_sizer.cu"
    "${source_path}/cuda_implementation/cuda_manager.cpp"
    "${source_path}/cuda_implementation/cuda_worker.cpp"
    "${source_path}/cuda_implementation/cuda_random.cu"
    "${source_path}/cuda_implementation/evaluate_fitness.cu"
    "${source_path}/cuda_implementation/virtual_screen.cu"
  )
endif()

if(MUDOCK_ENABLE_POLY)
  set(poly_header_path "${CMAKE_CURRENT_SOURCE_DIR}/../polygeist/include")
  set(poly_source_path "${CMAKE_CURRENT_SOURCE_DIR}/../polygeist/src")
  list(APPEND header_files
    "${header_path}/${header_prefix}/cuda_implementation/cuda_batch_sizer.cuh"
    "${header_path}/${header_prefix}/cuda_implementation/cuda_check_error_macro.cuh"
    "${header_path}/${header_prefix}/cuda_implementation/cuda_manager.hpp"
    "${header_path}/${header_prefix}/cuda_implementation/cuda_object.cuh"
    "${header_path}/${header_prefix}/cuda_implementation/cuda_worker.hpp"
    "${header_path}/${header_prefix}/cuda_implementation/cuda_wrapper.cuh"
    "${header_path}/${header_prefix}/cuda_implementation/map_textures.cuh"
    "${poly_header_path}/polygeist/cuda_random.cuh"
    "${poly_header_path}/polygeist/virtual_screen.cuh"
  )
  list(APPEND source_files
    "${source_path}/cuda_implementation/cuda_manager.cpp"
    "${source_path}/cuda_implementation/cuda_worker.cpp"
  )
endif()

# declare the HIP sources
if(MUDOCK_ENABLE_HIP)
  list(APPEND header_files
    "${header_path}/${header_prefix}/hip_implementation/hip_batch_sizer.hpp"
    "${header_path}/${header_prefix}/hip_implementation/hip_check_error_macro.hpp"
    "${header_path}/${header_prefix}/hip_implementation/hip_manager.hpp"
    "${header_path}/${header_prefix}/hip_implementation/hip_random.hpp"
    "${header_path}/${header_prefix}/hip_implementation/hip_object.hpp"
    "${header_path}/${header_prefix}/hip_implementation/hip_worker.hpp"
    "${header_path}/${header_prefix}/hip_implementation/hip_wrapper.hpp"
    "${header_path}/${header_prefix}/hip_implementation/calc_energy.hpp"
    "${header_path}/${header_prefix}/hip_implementation/evaluate_fitness.hpp"
    "${header_path}/${header_prefix}/hip_implementation/geometric_transformations.hpp"
    "${header_path}/${header_prefix}/hip_implementation/mutate.hpp"
    "${header_path}/${header_prefix}/hip_implementation/virtual_screen.hpp"
  )
  list(APPEND source_files
    "${source_path}/hip_implementation/hip_batch_sizer.cpp"
    "${source_path}/hip_implementation/hip_manager.cpp"
    "${source_path}/hip_implementation/hip_worker.cpp"
    "${source_path}/hip_implementation/hip_random.cpp"
    "${source_path}/hip_implementation/evaluate_fitness.cpp"
    "${source_path}/hip_implementation/virtual_screen.cpp"
  )
endif()

# declare the SYCL sources
if(MUDOCK_ENABLE_SYCL)
  list(APPEND header_files
    "${header_path}/${header_prefix}/sycl_implementation/sycl_manager.hpp"
    "${header_path}/${header_prefix}/sycl_implementation/calc_energy.hpp"
    "${header_path}/${header_prefix}/sycl_implementation/mutate.hpp"
    "${header_path}/${header_prefix}/sycl_implementation/geometric_transformations.hpp"
    "${header_path}/${header_prefix}/sycl_implementation/sycl_random.hpp"
    "${header_path}/${header_prefix}/sycl_implementation/sycl_object.hpp"
    "${header_path}/${header_prefix}/sycl_implementation/sycl_wrapper.hpp"
    "${header_path}/${header_prefix}/sycl_implementation/sycl_worker.hpp"
    "${header_path}/${header_prefix}/sycl_implementation/sycl_batch_sizer.hpp"
    "${header_path}/${header_prefix}/sycl_implementation/virtual_screen.hpp"
    "${header_path}/${header_prefix}/sycl_implementation/evaluate_fitness.hpp"
  )
  list(APPEND source_files
    "${source_path}/sycl_implementation/sycl_manager.cpp"
    "${source_path}/sycl_implementation/sycl_random.cpp"
    "${source_path}/sycl_implementation/sycl_worker.cpp"
    "${source_path}/sycl_implementation/sycl_batch_sizer.cpp"
    "${source_path}/sycl_implementation/virtual_screen.cpp"
  )
endif()

# declare the OpenMP sources
if(MUDOCK_ENABLE_OMP)
  list(APPEND header_files
    "${header_path}/${header_prefix}/omp_implementation/omp_manager.hpp"
    "${header_path}/${header_prefix}/omp_implementation/calc_energy.hpp"
    "${header_path}/${header_prefix}/omp_implementation/mutate.hpp"
    "${header_path}/${header_prefix}/omp_implementation/geometric_transformations.hpp"
    # "${header_path}/${header_prefix}/omp_implementation/omp_random.hpp"
    "${header_path}/${header_prefix}/omp_implementation/omp_object.hpp"
    "${header_path}/${header_prefix}/omp_implementation/omp_wrapper.hpp"
    "${header_path}/${header_prefix}/omp_implementation/omp_worker.hpp"
    "${header_path}/${header_prefix}/omp_implementation/omp_batch_sizer.hpp"
    "${header_path}/${header_prefix}/omp_implementation/virtual_screen.hpp"
    "${header_path}/${header_prefix}/omp_implementation/evaluate_fitness.hpp"
  )
  list(APPEND source_files
    "${source_path}/omp_implementation/omp_manager.cpp"
    # "${source_path}/sycl_implementation/omp_random.cpp"
    "${source_path}/omp_implementation/omp_worker.cpp"
    "${source_path}/omp_implementation/omp_batch_sizer.cpp"
    "${source_path}/omp_implementation/virtual_screen.cpp"
  )
endif()

# ########################################################################
# ### Declare the compile step
# ########################################################################

# declare the mudock library
add_library(libmudock STATIC ${source_files} ${header_files})
add_library(mudock::libmudock ALIAS libmudock)
set_target_properties(libmudock PROPERTIES PUBLIC_HEADER "${headers}")
target_include_directories(libmudock PUBLIC "${header_path}")
if(MUDOCK_ENABLE_POLY)
  target_include_directories(libmudock PUBLIC "${poly_header_path}")
endif()
# Include Boost and RDKit directories for the target
target_include_directories(libmudock PRIVATE ${Boost_INCLUDE_DIRS})
target_include_directories(libmudock PRIVATE ${RDKit_INCLUDE_DIRS})
set_target_properties(libmudock PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(libmudock PROPERTIES PREFIX "" OUTPUT_NAME "libmudock")
set_target_properties(libmudock
  PROPERTIES
  CXX_STANDARD 20
  CXX_STANDARD_REQUIRED ON
  CXX_EXTENSIONS OFF
)
target_link_libraries(libmudock PUBLIC Microsoft.GSL::GSL)
target_compile_options(libmudock PUBLIC $<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANGUAGE:CXX>>:${release_flags}>)
add_rdkit_dep(libmudock)

# check if we need to activate the compilation flag that signals Polygeist transpilation
if(MUDOCK_ENABLE_POLY)
  target_compile_definitions(libmudock PUBLIC MUDOCK_USE_POLY)
  target_link_directories(libmudock PUBLIC ${MUDOCK_POLY_LIB_PATH})
  target_link_libraries(libmudock PRIVATE mudock_poly omp)
endif()

# check if we need to activate the compilation flag that signals HIP
if(MUDOCK_ENABLE_HIP)
  target_compile_definitions(libmudock PUBLIC MUDOCK_USE_HIP)
  target_link_libraries(libmudock PUBLIC hiprand)
  target_compile_options(libmudock PUBLIC $<$<AND:$<CONFIG:DEBUG>,$<COMPILE_LANGUAGE:CXX>>:-G -g -Xptxas -v --expt-relaxed-constexpr>)
  target_compile_options(libmudock PUBLIC $<$<AND:$<CONFIG:RELEASE>,$<COMPILE_LANGUAGE:CXX>>:-O3 -use_fast_math --expt-relaxed-constexpr>)
endif()

# check if we need to activate the compilation flag that signals CUDA
if(MUDOCK_ENABLE_CUDA)
  target_compile_definitions(libmudock PUBLIC MUDOCK_USE_CUDA)
  set_target_properties(libmudock
    PROPERTIES
    CUDA_STANDARD 20
    CUDA_STANDARD_REQUIRED ON
    CUDA_EXTENSIONS OFF
    CUDA_SEPARABLE_COMPILATION ON)
  set(CMAKE_CUDA_FLAGS "--expt-relaxed-constexpr")
  set(CMAKE_CUDA_FLAGS_RELEASE "-O3 -use_fast_math --expt-relaxed-constexpr")
  set(CMAKE_CUDA_FLAGS_DEBUG "-G -g -Xptxas -v --expt-relaxed-constexpr")
endif()

# check if we need to activate the compilation flag that signals OpenMP
if(MUDOCK_ENABLE_OMP)
  target_compile_definitions(libmudock PUBLIC MUDOCK_USE_OMP)
  target_link_libraries(libmudock PUBLIC OpenMP::OpenMP_CXX)

  # Enable offloading to AMD GPUs using Clang
  if(NOT OMP_TARGET EQUAL "cpu")
    if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
      if(GPU_PLATFORM STREQUAL "cuda")
        target_compile_options(libmudock PRIVATE "-fopenmp")
        target_compile_options(libmudock PRIVATE "-fopenmp-targets=nvptx64")
        # target_compile_options(libmudock PRIVATE "-Xopenmp-target=nvptx64")
      elseif(GPU_PLATFORM STREQUAL "amd")
        target_compile_options(libmudock PRIVATE "-fopenmp")
        target_compile_options(libmudock PRIVATE "-fopenmp-targets=amdgcn-amd-amdhsa")
        # target_compile_options(libmudock PRIVATE "-Xopenmp-target=amdgcn-amd-amdhsa")
      endif()
      target_compile_options(libmudock PRIVATE "--offload-arch=${GPU_ARCHITECTURE}")
    elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
      # TODO test with GCC, check also if GCC supports this
      message(ERROR "muDock does not support GCC and OpenMP offloading")
      # if(GPU_PLATFORM STREQUAL "cuda")
      #   target_compile_options(libmudock PRIVATE "-fopenmp")
      #   target_compile_options(libmudock PRIVATE "-foffload=nvptx-none")
      #   target_compile_options(libmudock PRIVATE "-foffload-options=nvptx-none=-march=${GPU_ARCHITECTURE}")
      # elseif(GPU_PLATFORM STREQUAL "amd")
      #   target_compile_options(libmudock PRIVATE "-fopenmp")
      #   target_compile_options(libmudock PRIVATE "-foffload=amdgcn-amdhsa")
      #   target_compile_options(libmudock PRIVATE "-foffload-options=amdgcn-amdhsa=-march=${GPU_ARCHITECTURE}")
      # endif()
    else()
      message(ERROR "Not support OpenMP compiler")
    endif()
  endif()
endif()

# check if we need to activate the compilation flag that signals SYCL
if(MUDOCK_ENABLE_SYCL)
  target_compile_definitions(libmudock PUBLIC MUDOCK_USE_SYCL)

  # Link runtime device libraries (if needed)
  string(FIND "${ACPP_TARGETS}" "cuda" position)

  if(${position} GREATER -1)
    target_link_libraries(libmudock PUBLIC cudart)
  endif()

  # TODO only with AdaptiveCpp
  # Manually adding acpp library
  # NOTE: the compiler should be a clang version with sycl
  # target_link_libraries(libmudock PUBLIC acpp-rt)
  # list(APPEND SYCL_FLAGS "-fsycl" "-fsycl-targets=${ACPP_TARGETS}")
  # target_compile_options(libmudock PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${SYCL_FLAGS}>)
  # target_link_options(libmudock PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${SYCL_FLAGS}>)

  # Automatic configuration
  # NOTE: the compiler will be set by default to the AdaptiveCpp one
  # NOTE: AdaptiveCpp Single pass compiler optimizations are not required
  add_sycl_to_target(TARGET libmudock SOURCES)
endif()

# ########################################################################
# ### Declare the install step
# ########################################################################

# at the moment we don't deploy the library
