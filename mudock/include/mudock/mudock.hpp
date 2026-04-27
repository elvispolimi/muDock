#pragma once

#if defined(MUDOCK_USE_ALPAKA) &&                                                       \
    (defined(MUDOCK_ALPAKA_BACKEND_SERIAL) || defined(MUDOCK_ALPAKA_BACKEND_THREADS) || \
     defined(MUDOCK_ALPAKA_BACKEND_TBB) || defined(MUDOCK_ALPAKA_BACKEND_OMP2) ||       \
     defined(__CUDACC__) || defined(__HIPCC__))
  #include <mudock/alpaka_implementation.hpp>
#endif
#include <mudock/batch.hpp>
#include <mudock/chem.hpp>
#include <mudock/compute.hpp>
#include <mudock/format.hpp>
#include <mudock/grid.hpp>
#include <mudock/implementations.hpp>
#include <mudock/log.hpp>
#include <mudock/molecule.hpp>
#include <mudock/splitter.hpp>
#include <mudock/type_alias.hpp>
#include <mudock/utils.hpp>
