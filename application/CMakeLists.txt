#########################################################################
#### Declare the application source files
#########################################################################

# headers
set(header_path "${CMAKE_CURRENT_SOURCE_DIR}/src")
list(APPEND header_files
)

# sources
set(source_path "${CMAKE_CURRENT_SOURCE_DIR}/src")
list(APPEND source_files
  "${source_path}/main.cpp"
)

#########################################################################
#### Declare the compile step
#########################################################################

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
target_link_libraries("${CMAKE_PROJECT_NAME}" libmudock)

#########################################################################
#### Declare the install step
#########################################################################

# define the install step
install(TARGETS "${CMAKE_PROJECT_NAME}" RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
