#pragma once

#if defined(LIKWID_PERFMON) || defined(LIKWID_NVMON)
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

  #define NVMON_MARKER_INIT
  #define NVMON_MARKER_THREADINIT
  #define NVMON_MARKER_SWITCH
  #define NVMON_MARKER_REGISTER(regionTag)
  #define NVMON_MARKER_START(regionTag)
  #define NVMON_MARKER_STOP(regionTag)
  #define NVMON_MARKER_CLOSE
  #define NVMON_MARKER_GET(regionTag, nevents, events, time, count)
#endif

#define MUDOCK_MARKER_INIT
#define MUDOCK_MARKER_CLOSE
#define MUDOCK_MARKER_THREADINIT
#ifdef MUDOCK_LIKWID_CPP
  #if defined(MUDOCK_LIKWID_CPP) && !defined(LIKWID_PERFMON)
    #error "Requested LIKWID CPP backend, but LIKWID_PERFMON is not set"
  #endif
  #undef MUDOCK_MARKER_INIT
  #undef MUDOCK_MARKER_CLOSE
  #undef MUDOCK_MARKER_THREADINIT

  #define MUDOCK_MARKER_INIT                    LIKWID_MARKER_INIT
  #define MUDOCK_MARKER_CLOSE                   LIKWID_MARKER_CLOSE
  #define MUDOCK_MARKER_THREADINIT              LIKWID_MARKER_THREADINIT
  #define MUDOCK_CPP_MARKER_SWITCH              LIKWID_MARKER_SWITCH
  #define MUDOCK_CPP_MARKER_REGISTER(regionTag) LIKWID_MARKER_REGISTER(regionTag)
  #define MUDOCK_CPP_MARKER_START(regionTag)    LIKWID_MARKER_START(regionTag)
  #define MUDOCK_CPP_MARKER_STOP(regionTag)     LIKWID_MARKER_STOP(regionTag)
  #define MUDOCK_CPP_MARKER_GET(regionTag, nevents, events, time, count) \
    LIKWID_MARKER_GET(regionTag, nevents, events, time, count)
#else
  #define MUDOCK_CPP_MARKER_SWITCH
  #define MUDOCK_CPP_MARKER_REGISTER(regionTag)
  #define MUDOCK_CPP_MARKER_START(regionTag)
  #define MUDOCK_CPP_MARKER_STOP(regionTag)
  #define MUDOCK_CPP_MARKER_GET(regionTag, nevents, events, time, count)
#endif
#ifdef MUDOCK_LIKWID_CUDA
  #if defined(MUDOCK_LIKWID_CUDA) && !defined(LIKWID_NVMON)
    #error "Requested LIKWID CUDA backend, but LIKWID_NVMON is not set"
  #endif
  #undef MUDOCK_MARKER_INIT
  #undef MUDOCK_MARKER_CLOSE
  #undef MUDOCK_MARKER_THREADINIT

  #define MUDOCK_MARKER_INIT                     NVMON_MARKER_INIT
  #define MUDOCK_MARKER_CLOSE                    NVMON_MARKER_CLOSE
  #define MUDOCK_MARKER_THREADINIT               NVMON_MARKER_THREADINIT
  #define MUDOCK_CUDA_MARKER_SWITCH              NVMON_MARKER_SWITCH
  #define MUDOCK_CUDA_MARKER_REGISTER(regionTag) NVMON_MARKER_REGISTER(regionTag)
  #define MUDOCK_CUDA_MARKER_START(regionTag)    NVMON_MARKER_START(regionTag)
  #define MUDOCK_CUDA_MARKER_STOP(regionTag)     NVMON_MARKER_STOP(regionTag)
  #define MUDOCK_CUDA_MARKER_GET(regionTag, nevents, events, time, count) \
    NVMON_MARKER_GET(regionTag, nevents, events, time, count)
#else
  #define MUDOCK_CUDA_MARKER_REGISTER(regionTag)
  #define MUDOCK_CUDA_MARKER_START(regionTag)
  #define MUDOCK_CUDA_MARKER_STOP(regionTag)
  #define MUDOCK_CUDA_MARKER_CLOSE
  #define MUDOCK_CUDA_MARKER_GET(regionTag, nevents, events, time, count)
#endif
