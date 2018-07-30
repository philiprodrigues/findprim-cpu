#ifndef BENCHMARKING_H
#define BENCHMARKING_H

#include "Timer.h"

#include <thread>

typedef void (*processing_fn)(TPCData*, int, int, int);

// ======================================================================
std::vector<std::pair<float, float>> timefnThreads(processing_fn fn,
                                                   TPCData* data,
                                                   int NTH)

{
    std::vector<std::pair<float, float>> ret;

    int nsize=data->nchannels * data->nsamples;

    Timer t;
    const float dataSizeGB=float(nsize)*sizeof(SAMPLE_TYPE)/(1024*1024*1024);
    const int nrun=10;
    long long totTime=0;

    // TODO: How to deal with this nicely?

    // Each thread has its own place to write output hits
    // unsigned short *hits[NTH];
    // for(int i=0; i<NTH; ++i) hits[i]=new unsigned short[nsize/10];

    // We process channels in "blocks" of 16 (because that's the
    // number of 16-bit shorts that fit in a 256-bit register)
    const int nblock=data->nchannels/16;
    const int blocks_per_thread=nblock/NTH;

    for(int j=0; j<nrun; ++j){
        // 2018-07-30: There should be no overwriting of the inputs any more

        // The processing overwrites the input, so give each run its own copy
        // SAMPLE_TYPE* my_src      = (SAMPLE_TYPE*)malloc(nsize*sizeof(SAMPLE_TYPE));
        // memcpy(my_src, src, nsize*sizeof(SAMPLE_TYPE));

        t.reset();
        std::vector<std::thread> threads;
        for(int ith=0; ith<NTH; ++ith){
            int first_chan=ith*blocks_per_thread*16;
            int last_chan=std::min(data->nchannels+1, (ith+1)*blocks_per_thread*16);
            int hit_offset=(nsize/(NTH+1))*ith;

            threads.emplace_back(std::thread{fn,
                        data, first_chan, last_chan, hit_offset});
        }
        for(int ith=0; ith<NTH; ++ith){
            threads[ith].join();
        }
        long long runTime=t.getTime();
        totTime+=runTime;
        const float speed=dataSizeGB/(1e-6*runTime);
        // printf("% 6.1fms, % 9.0f GB/s ", 1e-3*totTime/nrun, speed);
        // printf("\n");
        ret.push_back(std::make_pair(runTime, speed));
        // free(my_src);
        data->zeroOutput();
    }
    // for(int i=0; i<NTH; ++i) delete[] hits[i];

    return ret;
}

//======================================================================
void timefnNThreads(const char* name,
                    processing_fn fn,
                    TPCData* data,
                    const std::vector<unsigned int>& nthreads)
{
    // int nsize=data->nchannels*data->nsamples;
    // const float msData=data->msData();
    const float APAmsData=data->APAmsData();
    const float dataSizeGB=data->dataSizeGB();

    std::vector<std::pair<float, float> > times;
    std::vector<float> minTimes, maxTimes;
    printf("% 20s ", name);

    std::ofstream fout(std::string(name)+"-times");
    for(size_t i=0; i<nthreads.size(); ++i){
        fout << nthreads[i] << " ";
        std::vector<std::pair<float, float> > nthread_times=timefnThreads(fn, data, nthreads[i]);
        float totTime=0;
        float minTime=std::numeric_limits<float>::max();
        float maxTime=0;
        for(auto const& p: nthread_times){
            totTime+=p.first;
            fout << p.first << " ";
            minTime=std::min(minTime, p.first);
            maxTime=std::max(maxTime, p.first);
        }
        fout << std::endl;
        const unsigned int nrun=nthread_times.size();
        times.push_back(std::make_pair(1e-3*totTime/nrun, dataSizeGB/(1e-6*totTime/nrun)));
        minTimes.push_back(1e-3*minTime);
        maxTimes.push_back(1e-3*maxTime);
    }
    // ---------------------------------------------
    // Print min times
    for(size_t i=0; i<nthreads.size(); ++i){
        printf("% 8.1f ", minTimes[i]);
    }
    printf(" min ms\n");

    // ---------------------------------------------
    // Print average times
    printf("% 20s ", "");
    for(size_t i=0; i<nthreads.size(); ++i){
        printf("% 8.1f ", times[i].first);
    }
    printf(" avg ms\n");

    // ---------------------------------------------
    // Print max times
    printf("% 20s ", "");
    for(size_t i=0; i<nthreads.size(); ++i){
        printf("% 8.1f ", maxTimes[i]);
    }
    printf(" max ms\n");

    // ---------------------------------------------
    // Print "APA/server"
    printf("% 20s ", "");
    for(size_t i=0; i<nthreads.size(); ++i){
        printf("% 8.1f ", APAmsData/times[i].first);
    }
    printf(" APA/server\n");

    // ---------------------------------------------
    // Print average speeds
    printf("% 20s ", "");
    for(size_t i=0; i<nthreads.size(); ++i){
        printf("% 8.1f ", times[i].second);
    }
    printf(" GB/s\n");

    printf("\n");
}

#endif
