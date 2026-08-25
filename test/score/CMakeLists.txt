# Define a function to add tests for each subfolder
function(add_adt_tests dir)
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
      add_test(NAME test_score_adt_${subdir_name}
               COMMAND test_score_adt --dpf ${subdir}/${subdir_name}_epdb.dpf
      ) # You can customize this
    endif()
  endforeach()
endfunction()
add_adt_tests("${CMAKE_SOURCE_DIR}/data")

# X-Score tests: compare muDock's four ported terms and the affinity they predict against
# values from the original XScore v1.3 (see script/gen_xscore_reference.sh).
#
# Called on two directories.
# `data` (AutoDock tests, which contain extra files).
# `data_xscore` added purely for X-Score coverage (no need for extra test files, would fail if run on AutoDock's tests).
# A complex without a <name>_xscore.ref is reported as skipped, not failed.
function(add_x_score_tests dir)
  file(
    GLOB directories
    LIST_DIRECTORIES true
    ${dir}/*)

  foreach(subdir ${directories})
    if(IS_DIRECTORY ${subdir})
      get_filename_component(subdir_name ${subdir} NAME)
      add_test(
        NAME test_score_x_${subdir_name}
        COMMAND
          test_score_x -p ${subdir}/${subdir_name}_protein.pdb -l
          ${subdir}/${subdir_name}_ligand.mol2 -r
          ${subdir}/${subdir_name}_xscore.ref)
      set_tests_properties(test_score_x_${subdir_name}
                           PROPERTIES SKIP_RETURN_CODE ${MUDOCK_CTEST_SKIP_RETURN_CODE})
    endif()
  endforeach()
endfunction()

add_x_score_tests("${CMAKE_SOURCE_DIR}/data")
add_x_score_tests("${CMAKE_SOURCE_DIR}/data_xscore")

function(add_tests dir backend device id)
  # List all subdirectories in the given directory
  file(
    GLOB directories
    LIST_DIRECTORIES true
    ${dir}/*)

  # Iterate through the subdirectories
  foreach(subdir ${directories})
    if(IS_DIRECTORY ${subdir}) # Make sure it's a directory
      get_filename_component(subdir_name ${subdir} NAME) # Get the directory
      set(ligand_path "${subdir}/${subdir_name}_ligand.mol2")
      if(EXISTS "${subdir}/${subdir_name}_ligand.adtmol2")
        set(ligand_path "${subdir}/${subdir_name}_ligand.adtmol2")
      endif()
      # Define the test executable for the subdirectory (e.g., test_subfolder)
      add_test(
        NAME test_score_adt_impl_${backend}_${subdir_name}
        COMMAND
          test_score_adt_impl -p ${subdir}/${subdir_name}_protein.pdb -l
          ${ligand_path} --use ${backend}:${device}:${id})
    endif()
  endforeach()
endfunction()

if(MUDOCK_ENABLE_GH)
  add_tests("${CMAKE_SOURCE_DIR}/data" "GH" "CPU" "0")
endif()
if(MUDOCK_ENABLE_XSIMD)
  add_tests("${CMAKE_SOURCE_DIR}/data" "XSIMD" "CPU" "0")
endif()
if(MUDOCK_ENABLE_CUDA)
  add_tests("${CMAKE_SOURCE_DIR}/data" "CUDA" "GPU" "0")
endif()
if(MUDOCK_ENABLE_HIP)
  add_tests("${CMAKE_SOURCE_DIR}/data" "HIP" "GPU" "0")
endif()
if(MUDOCK_ENABLE_SYCL)
  add_tests("${CMAKE_SOURCE_DIR}/data" "SYCL" "GPU" "0")
endif()