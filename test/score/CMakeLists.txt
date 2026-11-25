# Define a function to add tests for each subfolder
function(add_score_tests dir)
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

      # --- ADT ----------------------------------------------------------------
      add_test(NAME test_score_adt_${subdir_name}
               COMMAND test_score_adt --dpf ${subdir}/${subdir_name}_epdb.dpf
      ) # You can customize this

      # -- ADT ----------------------------------------------------------------
      add_test(NAME test_score_ds_${subdir_name}
               COMMAND test_score_ds -p ${subdir}/${subdir_name}_protein.pdb
               -l ${subdir}/${subdir_name}_ligand.mol2
               -s ${subdir}/${subdir_name}_ligand_score.ds
      )
    endif()
  endforeach()
endfunction()
add_score_tests("${CMAKE_SOURCE_DIR}/data")

function(add_impl_tests dir backend device id)
  # List all subdirectories in the given directory
  file(
    GLOB directories
    LIST_DIRECTORIES true
    ${dir}/*)

  # Iterate through the subdirectories
  foreach(subdir ${directories})
    if(IS_DIRECTORY ${subdir}) # Make sure it's a directory
      get_filename_component(subdir_name ${subdir} NAME) # Get the directory
      
      # --- ADT ----------------------------------------------------------------
      set(adt_ligand_path "${subdir}/${subdir_name}_ligand.mol2")
      if(EXISTS "${subdir}/${subdir_name}_ligand.adtmol2")
        set(adt_ligand_path "${subdir}/${subdir_name}_ligand.adtmol2")
      endif()
      # Define the test executable for the subdirectory (e.g., test_subfolder)
      add_test(
        NAME test_score_adt_impl_${backend}_${subdir_name}
        COMMAND
          test_score_adt_impl -p ${subdir}/${subdir_name}_protein.pdb -l
          ${adt_ligand_path} --use ${backend}:${device}:${id}
      )

      # --- DS -----------------------------------------------------------------
      set(ds_ligand_path "${subdir}/${subdir_name}_ligand.mol2")
      add_test(
        NAME test_score_ds_impl_${backend}_${subdir_name}
        COMMAND
          test_score_ds_impl -p ${subdir}/${subdir_name}_protein.pdb -l
          ${ds_ligand_path} --use ${backend}:${device}:${id}
      )
    endif()
  endforeach()
endfunction()

if(MUDOCK_ENABLE_GH)
  add_impl_tests("${CMAKE_SOURCE_DIR}/data" "GH" "CPU" "0")
endif()
if(MUDOCK_ENABLE_XSIMD)
  add_impl_tests("${CMAKE_SOURCE_DIR}/data" "XSIMD" "CPU" "0")
endif()
if(MUDOCK_ENABLE_CUDA)
  add_impl_tests("${CMAKE_SOURCE_DIR}/data" "CUDA" "GPU" "0")
endif()
if(MUDOCK_ENABLE_HIP)
  add_impl_tests("${CMAKE_SOURCE_DIR}/data" "HIP" "GPU" "0")
endif()
if(MUDOCK_ENABLE_SYCL)
  add_impl_tests("${CMAKE_SOURCE_DIR}/data" "SYCL" "GPU" "0")
endif()
