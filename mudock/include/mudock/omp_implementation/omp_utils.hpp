#include <math.h>
#include <mudock/type_alias.hpp>

namespace mudock {
  typedef struct {
    mudock::fp_type value;
    int index;
  } min_index_pair;

#pragma omp declare reduction(min_index:min_index_pair : omp_out =                   \
                                  (omp_in.value < omp_out.value) ? omp_in : omp_out) \
    initializer(omp_priv = {INFINITY, -1})

} // namespace mudock
