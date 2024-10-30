#pragma once

#ifdef MUDOCK_USE_SCOREP
  #include <scorep/SCOREP_User.h>
  #define SCOREP_DEFINE_MARKER(regionTag) SCOREP_USER_REGION_DEFINE(regionTag)
  #define SCOREP_MARKER_START(regionTag, name) \
    SCOREP_USER_REGION_BEGIN(regionTag, name, SCOREP_USER_REGION_TYPE_COMMON)
  #define SCOREP_MARKER_STOP(regionTag) SCOREP_USER_REGION_END(regionTag)
#else
  #define SCOREP_DEFINE_MARKER(regionTag)
  #define SCOREP_MARKER_START(regionTag, name)
  #define SCOREP_MARKER_STOP(regionTag)
#endif

SCOREP_DEFINE_MARKER(ga)
