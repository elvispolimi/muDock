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
      # name
      message(STATUS "${subdir}")
      # Define the test executable for the subdirectory (e.g., test_subfolder)
      add_test(NAME test_score_${subdir_name}
               COMMAND test_score --dpf ${subdir}/${subdir_name}_epdb.dpf
      )# You can customize this
    endif()
  endforeach()
endfunction()

add_tests("${CMAKE_SOURCE_DIR}/data")

if(MUDOCK_ENABLE_GH)

endif()
