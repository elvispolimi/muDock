#pragma once

#include <mudock/compute/adt_score_args.hpp>
#include <mudock/compute/algorithm.hpp>
#include <mudock/cpp_implementation/adt_score_cpp.hpp>
#include <mudock/cpp_implementation/queue_cpp.hpp>

#ifdef MUDOCK_USE_GH
  #include <mudock/gh_implementation/adt_score_gh.hpp>
  #include <mudock/gh_implementation/queue_gh.hpp>
#endif

#ifdef MUDOCK_USE_XSIMD
  #include <mudock/xsimd_implementation/adt_score_xsimd.hpp>
  #include <mudock/xsimd_implementation/queue_xsimd.hpp>
#endif

namespace mudock {
  template<>
  struct genetic_scoring<queue_cpp, scoring_function::ADT> {
    static void score(const adt_score_args& args) { score_batch_cpp(args); }
  };

#ifdef MUDOCK_USE_GH
  template<>
  struct genetic_scoring<queue_gh, scoring_function::ADT> {
    static void score(const adt_score_args& args) { score_batch_gh(args); }
  };
#endif

#ifdef MUDOCK_USE_XSIMD
  template<>
  struct genetic_scoring<queue_xsimd, scoring_function::ADT> {
    static void score(const adt_score_args& args) { score_batch_xsimd(args); }
  };
#endif
} // namespace mudock
