#define SAMPLE_TYPE_SHORT

#ifdef SAMPLE_TYPE_SHORT
#define SAMPLE_TYPE short
#else
#define SAMPLE_TYPE float
#endif

#include "process_samples.h"
#include "run_algs_avx2.h"

#include <cstdio>

int main()
{
    // How many times to repeat the data in memory so as to make sure
    // it's bigger than the cache
    const int nrepeat=160;

    const bool transpose=false;
    const char* filename="waveforms-01397060-sigonly-collection";
    printf("Loading...\n");
    TPCData data(filename, nrepeat, transpose, false);
    data.setTaps(8, 100);
    printf("Running algorithm...\n");
    for(int i=0; i<10; ++i){
        run_algs_avx2<false, false, false, 8>(&data, 0, data.nchannels, 0);
    }
}
