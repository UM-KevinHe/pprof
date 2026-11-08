#ifndef PPROFR_OMP_H
#define PPROFR_OMP_H

// OpenMP is optional. Compilers without it (e.g. Apple clang) get no-op
// stand-ins and the parallel loops simply run on one thread.
#ifdef _OPENMP
#include <omp.h>
#else
inline void omp_set_num_threads(int) {}
inline int omp_get_max_threads() { return 1; }
#endif

#endif
