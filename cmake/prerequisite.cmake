# ##############################################################################
# Grab the GSL library
# ##############################################################################

# TODO remove GSL depndencies
include(FetchContent)
FetchContent_Declare(
  GSL
  GIT_REPOSITORY "https://github.com/microsoft/GSL"
  GIT_TAG "v4.0.0"
  GIT_SHALLOW ON)
FetchContent_MakeAvailable(GSL)
