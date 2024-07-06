#include <stdlib.h>
#include <unistd.h>
#include <climits>
#include <cstdint>

bool enable_vect_dp_chaining = false;
uint64_t avg;
uint64_t minimizer_lookup_time, alignment_time, dp_time, rmq_time, rmq_t1, rmq_t2, rmq_t3, rmq_t4;