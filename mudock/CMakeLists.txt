#########################################################################
#### Declare the application source files
#########################################################################

# headers
set(header_path "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(header_prefix "mudock")
list(APPEND header_files
  "${header_path}/${header_prefix}/molecule.hpp"
  "${header_path}/${header_prefix}/type_alias.hpp"
)

# sources
set(source_path "${CMAKE_CURRENT_SOURCE_DIR}/src")
list(APPEND source_files
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

#########################################################################
#### Declare the install step
#########################################################################

# at the moment we don't deploy the library
