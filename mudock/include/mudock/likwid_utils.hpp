#pragma once

#ifdef LIKWID_PERFMON
  #include <likwid-marker.h>
#else
  #define LIKWID_MARKER_INIT
  #define LIKWID_MARKER_THREADINIT
  #define LIKWID_MARKER_SWITCH
  #define LIKWID_MARKER_REGISTER(regionTag)
  #define LIKWID_MARKER_START(regionTag)
  #define LIKWID_MARKER_STOP(regionTag)
  #define LIKWID_MARKER_CLOSE
  #define LIKWID_MARKER_GET(regionTag, nevents, events, time, count)
#endif
