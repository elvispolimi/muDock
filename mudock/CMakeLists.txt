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
  "${header_path}/${header_prefix}/chem/bond_types.hpp"
  "${header_path}/${header_prefix}/chem/elements.hpp"
  "${header_path}/${header_prefix}/format/concept.hpp"
  "${header_path}/${header_prefix}/format/dot.hpp"
  "${header_path}/${header_prefix}/format/mol2.hpp"
  "${header_path}/${header_prefix}/format/pdb.hpp"
  "${header_path}/${header_prefix}/format/rdkit_wrapper.hpp"
  "${header_path}/${header_prefix}/grid/mdindex.hpp"
  "${header_path}/${header_prefix}/grid/pi.hpp"
  "${header_path}/${header_prefix}/grid/point3D.hpp"
  "${header_path}/${header_prefix}/molecule/bond.hpp"
  "${header_path}/${header_prefix}/molecule/constraints.hpp"
  "${header_path}/${header_prefix}/molecule/containers.hpp"
  "${header_path}/${header_prefix}/molecule/fragments.hpp"
  "${header_path}/${header_prefix}/molecule/graph.hpp"
  "${header_path}/${header_prefix}/molecule/properties.hpp"
  "${header_path}/${header_prefix}/molecule/property_table.hpp"
  "${header_path}/${header_prefix}/chem.hpp"
  "${header_path}/${header_prefix}/format.hpp"
  "${header_path}/${header_prefix}/grid.hpp"
  "${header_path}/${header_prefix}/molecule.hpp"
  "${header_path}/${header_prefix}/mudock.hpp"
  "${header_path}/${header_prefix}/splitter.hpp"
  "${header_path}/${header_prefix}/type_alias.hpp"

  "${header_path}/${header_prefix}/controllers.hpp"
  "${header_path}/${header_prefix}/devices.hpp"
  "${header_path}/${header_prefix}/utils.hpp"

  "${header_path}/${header_prefix}/controllers/cpp_thread_pool.hpp"
  "${header_path}/${header_prefix}/controllers/squeue.hpp"
  "${header_path}/${header_prefix}/controllers/thread_pool.hpp"
  "${header_path}/${header_prefix}/controllers/worker.hpp"

  "${header_path}/${header_prefix}/devices/cpp/context.hpp"
  "${header_path}/${header_prefix}/devices/cpp/kernel.hpp"
  "${header_path}/${header_prefix}/devices/cpp/scratch.hpp"

  "${header_path}/${header_prefix}/devices/buffer.hpp"
  "${header_path}/${header_prefix}/devices/context.hpp"
  "${header_path}/${header_prefix}/devices/kernel.hpp"
  "${header_path}/${header_prefix}/devices/parser.hpp"
  "${header_path}/${header_prefix}/devices/device_types.hpp"
  "${header_path}/${header_prefix}/devices/kernel_types.hpp"
  "${header_path}/${header_prefix}/devices/scratch.hpp"
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
  "${source_path}/format/dot.cpp"
  "${source_path}/format/mol2.cpp"
  "${source_path}/format/rdkit_wrapper.cpp"
  "${source_path}/molecule/fragments.cpp"
  "${source_path}/molecule/graph.cpp"
  "${source_path}/molecule/property_table.cpp"

  "${source_path}/controllers/cpp_thread_pool.cpp"

  "${source_path}/devices/cpp/kernel.cpp"
  "${source_path}/devices/cpp/context.cpp"
  "${source_path}/devices/parser.cpp"
  "${source_path}/devices/device_types.cpp"
  "${source_path}/devices/kernel_types.cpp"
)

# ########################################################################
# ### Declare the compile step
# ########################################################################

# declare the mudock library
add_library(libmudock STATIC ${source_files} ${header_files})
add_library(mudock::libmudock ALIAS libmudock)
set_target_properties(libmudock PROPERTIES PUBLIC_HEADER "${headers}")
target_include_directories(libmudock PUBLIC "${header_path}")
set_target_properties(libmudock PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(libmudock PROPERTIES PREFIX "" OUTPUT_NAME "libmudock")
set_target_properties(libmudock
  PROPERTIES
  CXX_STANDARD 20
  CXX_STANDARD_REQUIRED ON
  CXX_EXTENSIONS OFF
)
target_link_libraries(libmudock PUBLIC Microsoft.GSL::GSL)
add_rdkit_dep(libmudock)

# ########################################################################
# ### Declare the install step
# ########################################################################

# at the moment we don't deploy the library
