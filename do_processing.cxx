#define SAMPLE_TYPE_SHORT

#ifdef SAMPLE_TYPE_SHORT
#define SAMPLE_TYPE short
#else
#define SAMPLE_TYPE float
#endif

#include "process_samples.h"
#include "benchmarking.h"
#include "run_algs_naive.h"
#include "run_algs_avx2.h"
#include "constants.h"

#include <cstdio>
#include <cstring> // For memset
#include <limits>
#include <iomanip>
#include <vector>

#include <stdint.h>
#include "../waveform-tools/read_samples.h"
#include "design_fir.h"


#include "boost/program_options.hpp"

namespace po = boost::program_options;

const std::vector<unsigned int> nthreads{1, 2, 4, 8, 16, 32, 64};

// Just print the first 8 values
void print256(__m256i var)
{
    int16_t *val = (int16_t*) &var;
    printf("% 5i % 5i % 5i % 5i % 5i % 5i % 5i % 5i ",
           val[0], val[1], val[2], val[3], val[4], val[5],
           val[6], val[7]);
    val+=8;
    printf("% 5i % 5i % 5i % 5i % 5i % 5i % 5i % 5i",
           val[0], val[1], val[2], val[3], val[4], val[5],
           val[6], val[7]);
}


    

//======================================================================
void saveNaiveHitsToFile(unsigned short* hits, const char* filename)
{
    std::ofstream fout(filename);
    unsigned short chan, hit_start, hit_charge, hit_tover;
    unsigned short* input_loc=hits;
    int nhit=0;
    while(true){
        chan        = *input_loc++;
        hit_start   = *input_loc++;
        hit_charge  = *input_loc++;
        hit_tover   = *input_loc++;
        // Four magic values indicates the end of hits
        if(chan==MAGIC && hit_start==MAGIC &&
           hit_charge==MAGIC && hit_tover==MAGIC) break;
        fout << std::setw(6) << chan << " "
             << std::setw(6) << hit_start << " "
             << std::setw(6) << hit_charge << " "
             << std::setw(6) << hit_tover << std::endl;
        ++nhit;
    }
    printf("naive: Saved %d hits\n", nhit);
}

//======================================================================
void saveIntrinHitsToFile(unsigned short* hits, const char* filename)
{
    std::ofstream fout(filename);
    unsigned short chan[16], hit_start[16], hit_charge[16], hit_tover[16];
    unsigned short* input_loc=hits;
    int nhit=0;

    while(true){
        for(int i=0; i<16; ++i) chan[i]       = *input_loc++;
        for(int i=0; i<16; ++i) hit_start[i]  = *input_loc++;
        for(int i=0; i<16; ++i) hit_charge[i] = *input_loc++;
        for(int i=0; i<16; ++i) hit_tover[i]  = *input_loc++;

        for(int i=0; i<16; ++i){
            if(hit_charge[i]){
                fout << std::setw(6) << chan[i] << " "
                     << std::setw(6) << hit_start[i] << " "
                     << std::setw(6) << hit_charge[i] << " "
                     << std::setw(6) << hit_tover[i] << std::endl;
                ++nhit;
            }
        }
        if(*input_loc==MAGIC) break;
    }
    printf("intrinsic: Saved %d hits\n", nhit);
}

//======================================================================
void saveOutputToFile(const SAMPLE_TYPE* samples, int nchannels, int nsamples,
                      const char* filename, bool transpose)
{
    std::ofstream fout(filename);
    if(transpose){
        for(int j=0; j<nchannels; ++j){
            for(int i=0; i<nsamples; ++i){
                fout << samples[nchannels*i+j] << " ";
            }
            fout << std::endl;
        }
    }
    else{
        for(int i=0; i<nchannels; ++i){
            for(int j=0; j<nsamples; ++j){
                fout << samples[nsamples*i+j] << " ";
            }
            fout << std::endl;
        }
    }
}

//======================================================================
int main(int argc, char** argv)
{
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("input,i", po::value<std::string>(), "input file name")
        ("tag,t",   po::value<std::string>(), "output tag")
        ("ntaps,n", po::value<int>()->default_value(8), "number of filter taps (8 or 32)")
        ("store,s", "store intermediate values")
        ("nrepeat,r", po::value<int>()->default_value(16), "number of repeats of the data")
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    

    if(vm.count("help") || vm.empty()) {
        std::cout << desc << "\n";
        return 1;
    }

    if(!vm.count("input")){
        std::cout << "No input file specified" << std::endl;
        std::cout << desc << std::endl;
        return 1;
    }

    bool store_intermediate=vm.count("store");

    std::string inputfile=vm["input"].as<std::string>();

    // How many times to repeat the data in memory so as to make sure
    // it's bigger than the cache
    const int nrepeat=16;

    const bool transpose=false;
    TPCData data(inputfile.c_str(), nrepeat, transpose, store_intermediate);
    data.setTaps(vm["ntaps"].as<int>(), 100);
    //----------------------------------------------------------------
    printf("Using %d samples, %d collection channels (%d unique) = %.1f APA*ms with type size %ld bytes. Total size %.1f MB\n",
           data.nsamples, data.nchannels, data.nchannels_uniq, data.APAmsData(), sizeof(SAMPLE_TYPE),
           data.dataSizeMB());

    // ---------------------------------------------------------
    // Run for benchmarking
    // ---------------------------------------------------------

    printf("% 20s     Threads\n", "");
    printf("% 20s ", "");
    for(size_t i=0; i<nthreads.size(); ++i){
        printf("% 8d ", nthreads[i]);
    }
    printf("\n");

    timefnNThreads("avx2", run_algs_avx2<false, false, false, 8>,
                   &data, nthreads);
    timefnNThreads("naive", run_algs_naive<false, false, false, 8>,
                   &data, nthreads);

    // ---------------------------------------------------------
    // Run again to get the output
    // ---------------------------------------------------------

#if NTAPS==32
#define TAG "doprocessing-32tap"
#else
#ifdef COHERENT_PEDSUB
#define TAG "doprocessing-cohpedsub"
#else
#define TAG "doprocessing"
#endif
#endif

    run_algs_naive<false, false, false, 8>(&data, 0, data.nchannels, 0);
    saveNaiveHitsToFile(data.hits, "hits-" TAG "-naive");

    if(store_intermediate){
        saveOutputToFile(data.pedsub, data.nchannels_uniq, data.nsamples, "pedsub-" TAG "-naive", false);
        saveOutputToFile(data.filtered, data.nchannels_uniq, data.nsamples, "filtered-" TAG "-naive", false);
    }

    run_algs_avx2<false, false, false, 8>(&data, 0, data.nchannels, 0);
    saveIntrinHitsToFile(data.hits, "hits-" TAG "-intrin");

    if(store_intermediate){
        saveOutputToFile(data.pedsub, data.nchannels_uniq, data.nsamples, "pedsub-" TAG "-intrin", true);
        saveOutputToFile(data.filtered, data.nchannels_uniq, data.nsamples, "filtered-" TAG "-intrin", true);
    }

}
