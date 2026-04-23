# Define a function to add tests for each subfolder
function(add_tests dir)
  # List all subdirectories in the given directory
  file(
    GLOB directories
    LIST_DIRECTORIES true
    ${dir}/*)

  # Iterate through the subdirectories
  foreach(subdir ${directories})
    if(IS_DIRECTORY ${subdir}) # Make sure it's a directory
      get_filename_component(subdir_name ${subdir} NAME) # Get the directory
      # Define the test executable for the subdirectory (e.g., test_subfolder)
      add_test(NAME test_grid_${subdir_name}
               COMMAND test_grid --pdbqt ${subdir}/${subdir_name}_pocket.pdbqt
                       --autogrid ${subdir}/${subdir_name}_pocket.maps.fld
      )# You can customize this
      set_tests_properties(test_grid_${subdir_name}
                           PROPERTIES SKIP_RETURN_CODE ${MUDOCK_CTEST_SKIP_RETURN_CODE})
    endif()
  endforeach()
endfunction()

add_tests("${CMAKE_SOURCE_DIR}/data")
