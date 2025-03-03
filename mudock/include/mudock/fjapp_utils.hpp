#pragma once

#ifdef MUDOCK_USE_FJAPP
  #include <fj_tool/fapp.h>
  #define FJAPP_MARKER_START(regionTag) fapp_start(regionTag, 1, 0);
  #define FJAPP_MARKER_STOP(regionTag)  fapp_stop(regionTag, 1, 0);
#else
  #define FJAPP_MARKER_START(regionTag)
  #define FJAPP_MARKER_STOP(regionTag)
#endif
