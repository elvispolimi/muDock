#########################################################################
#### Declare the application source files
#########################################################################

# headers
set(header_path "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(header_prefix "mudock")
list(APPEND header_files
  "${header_path}/${header_prefix}/chem/bond_types.hpp"
  "${header_path}/${header_prefix}/chem/elements.hpp"
  "${header_path}/${header_prefix}/format/concept.hpp"
  "${header_path}/${header_prefix}/format/mol2.hpp"
  "${header_path}/${header_prefix}/format/pdb.hpp"
  "${header_path}/${header_prefix}/format/rdkit_wrapper.hpp"
  "${header_path}/${header_prefix}/grid/mdindex.hpp"
  "${header_path}/${header_prefix}/molecule/atom_coordinates.hpp"
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
)

# sources
set(source_path "${CMAKE_CURRENT_SOURCE_DIR}/src")
list(APPEND source_files
  "${source_path}/chem/bond_types.cpp"
  "${source_path}/chem/elements.cpp"
  "${source_path}/format/mol2.cpp"
  "${source_path}/format/rdkit_wrapper.cpp"
  "${source_path}/molecule/atom_coordinates.cpp"
  "${source_path}/molecule/fragments.cpp"
  "${source_path}/molecule/property_table.cpp"
  "${source_path}/molecule.cpp"
)

#########################################################################
#### Declare the compile step
#########################################################################

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

#########################################################################
#### Declare the install step
#########################################################################

# at the moment we don't deploy the library
