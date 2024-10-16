# ########################################################################
# ### Declare the application source files
# ########################################################################

# headers
set(header_path "${CMAKE_CURRENT_SOURCE_DIR}/src")
list(APPEND header_files
  "${header_path}/command_line_args.hpp"
)

# sources
set(source_path "${CMAKE_CURRENT_SOURCE_DIR}/src")
list(APPEND source_files
  "${header_path}/command_line_args.cpp"
  "${source_path}/main.cpp"
)

# ########################################################################
# ### Declare the compile step
# ########################################################################

# declare the application executable
add_executable("${CMAKE_PROJECT_NAME}" ${header_files} ${source_files})
target_include_directories("${CMAKE_PROJECT_NAME}" PRIVATE "${header_path}")
set_target_properties("${CMAKE_PROJECT_NAME}"
  PROPERTIES
  CXX_STANDARD 20
  CXX_STANDARD_REQUIRED ON
  CXX_EXTENSIONS OFF
)

# add the dependencies
if(MUDOCK_ENABLE_SYCL)
  string(FIND "${ACPP_TARGETS}" "omp" SUBSTRING_INDEX)
  if (SUBSTRING_INDEX GREATER -1)
    target_link_libraries("${CMAKE_PROJECT_NAME}" omp)
  endif()
  target_link_libraries("${CMAKE_PROJECT_NAME}" libmudock boost_program_options boost_context pthread boost_fiber m stdc++)
else()
  target_link_libraries("${CMAKE_PROJECT_NAME}" libmudock boost_program_options boost_context pthread m stdc++)
endif()

if(MUDOCK_ENABLE_OMP)
  target_compile_options("${CMAKE_PROJECT_NAME}" PRIVATE "-fopenmp-targets=nvptx64")
  target_link_options("${CMAKE_PROJECT_NAME}" PRIVATE "-fopenmp")
  target_link_options("${CMAKE_PROJECT_NAME}" PRIVATE "-fopenmp-targets=nvptx64")
endif()

# ########################################################################
# ### Declare the install step
# ########################################################################

# define the install step
install(TARGETS "${CMAKE_PROJECT_NAME}" RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
