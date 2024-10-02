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
  target_link_libraries("${CMAKE_PROJECT_NAME}" libmudock boost_program_options pthread boost_fiber m stdc++)
else()
  # target_link_libraries("${CMAKE_PROJECT_NAME}" libmudock Boost::program_options)
  target_link_libraries("${CMAKE_PROJECT_NAME}" libmudock boost_program_options pthread m stdc++)
endif()

# ########################################################################
# ### Declare the install step
# ########################################################################

# define the install step
install(TARGETS "${CMAKE_PROJECT_NAME}" RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
