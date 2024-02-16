#########################################################################
#### Grab the GSL library
#########################################################################

include(FetchContent)
FetchContent_Declare(GSL
    GIT_REPOSITORY "https://github.com/microsoft/GSL"
    GIT_TAG "v4.0.0"
    GIT_SHALLOW ON
)
FetchContent_MakeAvailable(GSL)


#########################################################################
#### Handle RDKit (it doesn't play nice with boost, we need to improvise)
#########################################################################

function(add_rdkit_dep target)
    get_target_property(RDKIT_LIB_FULLPATH RDKit::RDGeneral LOCATION)
    get_filename_component(RDKIT_LIB_DIRPATH "${RDKIT_LIB_FULLPATH}" DIRECTORY)
    cmake_path(GET RDKIT_LIB_DIRPATH PARENT_PATH RDKIT_INSTALL_PREFIX)
    target_include_directories(${target} SYSTEM PUBLIC "${RDKIT_INSTALL_PREFIX}/include/rdkit")
    target_link_directories(${target} PUBLIC "${RDKIT_INSTALL_PREFIX}/lib")
    target_link_directories(${target} PUBLIC "${RDKIT_INSTALL_PREFIX}/lib64")
    target_link_libraries(${target} PUBLIC
        RDKitFileParsers
        RDKitGraphMol
        RDKitRDGeneral
        RDKitSmilesParse
        RDKitSubstructMatch
    )
    message(STATUS "Found RDKit prefix: ${RDKIT_INSTALL_PREFIX}")
endfunction(add_rdkit_dep)
