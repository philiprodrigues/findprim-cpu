#ifndef PROCESS_NAIVE_H
#define PROCESS_NAIVE_H

#include "constants.hh"
#include "ProcessingInfo.h"

void frugal_accum_update(int16_t& m, const int16_t s, int16_t& acc, const int16_t acclimit);

void
process_window_naive(ProcessingInfo& info);

#endif

/* Local Variables:  */
/* mode: c++         */
/* c-basic-offset: 4 */
/* End:              */
