#define SAMPLE_TYPE_SHORT

#ifdef SAMPLE_TYPE_SHORT
#define SAMPLE_TYPE short
#else
#define SAMPLE_TYPE float
#endif

#ifndef NTAPS
#define NTAPS 8 // Keep this a power of two
#endif

#define "process_samples.h"

#include <immintrin.h>
#include <thread>
#include <cstdio>
#include <cstring> // For memset
#include <limits>
#include <iomanip>
#include <vector>

#include <stdint.h>
#include "../waveform-tools/read_samples.h"
#include "design_fir.h"
#include "Timer.h"

#include "boost/program_options.hpp"

namespace po = boost::program_options;

const unsigned short MAGIC = std::numeric_limits<unsigned short>::max();
const int threshold=600;

const int nnthreads=7;
const int nthreads[nnthreads]={1, 2, 4, 8, 16, 32, 64};

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
// std::vector<std::pair<float, float>> timefnThreads(processing_fn fn,
//                                                    TPCData* data,
//                                                    int NTH)

// {
//     std::vector<std::pair<float, float>> ret;

//     int nsize=data->nchannels * data->nsamples;

//     Timer t;
//     const float dataSizeGB=float(nsize)*sizeof(SAMPLE_TYPE)/(1024*1024*1024);
//     const int nrun=10;
//     long long totTime=0;
//     // Each thread has its own place to write output hits
//     unsigned short *hits[NTH];
//     for(int i=0; i<NTH; ++i) hits[i]=new unsigned short[nsize/10];

//     // We process channels in "blocks" of 16 (because that's the
//     // number of 16-bit shorts that fit in a 256-bit register)
//     const int nblock=nchannels/16;
//     const int blocks_per_thread=nblock/NTH;

//     for(int j=0; j<nrun; ++j){
//         // The processing overwrites the input, so give each run its own copy
//         SAMPLE_TYPE* my_src      = (SAMPLE_TYPE*)malloc(nsize*sizeof(SAMPLE_TYPE));
//         memcpy(my_src, src, nsize*sizeof(SAMPLE_TYPE));
//         t.reset();
//         std::vector<std::thread> threads;
//         for(int ith=0; ith<NTH; ++ith){
//             int first_chan=ith*blocks_per_thread*16;
//             int last_chan=std::min(nchannels+1, (ith+1)*blocks_per_thread*16);

//             threads.emplace_back(std::thread{fn,
//                         my_src, hits[ith],
//                         nchannels, nsamples,
//                         first_chan, last_chan,
//                         taps, ntaps,
//                         pedsub, filtered});
//         }
//         for(int ith=0; ith<NTH; ++ith){
//             threads[ith].join();
//         }
//         long long runTime=t.getTime();
//         totTime+=runTime;
//         const float speed=dataSizeGB/(1e-6*runTime);
//         // printf("% 6.1fms, % 9.0f GB/s ", 1e-3*totTime/nrun, speed);
//         // printf("\n");
//         ret.push_back(std::make_pair(runTime, speed));
//         free(my_src);
//     }
//     for(int i=0; i<NTH; ++i) delete[] hits[i];

//     return ret;
// }

// //======================================================================
// void timefnNThreads(const char* name,
//                     processing_fn fn,
//                     TPCData* data)
// {
//     int nsize=nchannels*nsamples;
//     const float msData=1e3*nsamples/sampling_rate;
//     const float APAmsData=msData*nchannels/ncollection_per_apa;
//     const float dataSizeGB=float(nsize)*sizeof(SAMPLE_TYPE)/(1024*1024*1024);

//     std::vector<std::pair<float, float> > times;
//     std::vector<float> minTimes, maxTimes;
//     printf("% 20s ", name);

//     std::ofstream fout(std::string(name)+"-times");
//     for(int i=0; i<nnthreads; ++i){
//         fout << nthreads[i] << " ";
//         std::vector<std::pair<float, float> > nthread_times=timefnThreads(fn, src,
//                                                                           nchannels, nsamples,
//                                                                           begin_chan, end_chan,
//                                                                           taps, ntaps,
//                                                                           pedsub, filtered, nthreads[i]);
//         float totTime=0;
//         float minTime=std::numeric_limits<float>::max();
//         float maxTime=0;
//         for(auto const& p: nthread_times){
//             totTime+=p.first;
//             fout << p.first << " ";
//             minTime=std::min(minTime, p.first);
//             maxTime=std::max(maxTime, p.first);
//         }
//         fout << std::endl;
//         const unsigned int nrun=nthread_times.size();
//         times.push_back(std::make_pair(1e-3*totTime/nrun, dataSizeGB/(1e-6*totTime/nrun)));
//         minTimes.push_back(1e-3*minTime);
//         maxTimes.push_back(1e-3*maxTime);
//     }
//     // ---------------------------------------------
//     // Print min times
//     for(int i=0; i<nnthreads; ++i){
//         printf("% 8.1f ", minTimes[i]);
//     }
//     printf(" min ms\n");

//     // ---------------------------------------------
//     // Print average times
//     printf("% 20s ", "");
//     for(int i=0; i<nnthreads; ++i){
//         printf("% 8.1f ", times[i].first);
//     }
//     printf(" avg ms\n");

//     // ---------------------------------------------
//     // Print max times
//     printf("% 20s ", "");
//     for(int i=0; i<nnthreads; ++i){
//         printf("% 8.1f ", maxTimes[i]);
//     }
//     printf(" max ms\n");

//     // ---------------------------------------------
//     // Print "APA/server"
//     printf("% 20s ", "");
//     for(int i=0; i<nnthreads; ++i){
//         printf("% 8.1f ", APAmsData/times[i].first);
//     }
//     printf(" APA/server\n");

//     // ---------------------------------------------
//     // Print average speeds
//     printf("% 20s ", "");
//     for(int i=0; i<nnthreads; ++i){
//         printf("% 8.1f ", times[i].second);
//     }
//     printf(" GB/s\n");

//     printf("\n");
// }

//======================================================================
void do_processing_naive(TPCData* data)
{
    unsigned short* output_loc=hits;
    int nhits=0;

    for(int ichan=data->begin_chan; ichan<data->end_chan; ichan++){

        // Variables for pedestal finding
#ifdef TRANSPOSE_MEMORY
        SAMPLE_TYPE median=data->src[ichan];
#else
        SAMPLE_TYPE median=data->src[ichan*nsamples];
#endif
        SAMPLE_TYPE accum=0;

        // Variables for filtering
        SAMPLE_TYPE prev_samp[NTAPS]={0};

        // Variables for hit finding
        bool prev_was_over=false; // was the previous sample over threshold?
        unsigned short hit_start=0; // start time of the hit
        unsigned short hit_charge=0;
        unsigned short hit_tover=0; // time over threshold

        for(int isample=0; isample<data->nsamples; ++isample){
#ifdef TRANSPOSE_MEMORY
            const size_t index=isample*data->nchannels+ichan;
#else
            const size_t index=ichan*data->nsamples+isample;
#endif
            // --------------------------------------------------------------
            // Pedestal finding/coherent noise removal
            // --------------------------------------------------------------

#ifdef COHERENT_PEDSUB
            // Use frugal streaming to estimate the median value of
            // the next 16 channels at this tick, and subtract it from
            // all of them. Sadly this just doesn't work: if the
            // starting value is on a big hit, then the median
            // estimation never recovers, and our estimate of the
            // median ends up way too large

            // SAMPLE_TYPE coh_median=src[index];
            // if(ichan%16==0){
            //     for(int iblock=0; iblock<16; ++iblock){
            //         int blockindex=(ichan+iblock)*nsamples+isample;
            //         if(src[blockindex]>coh_median) ++coh_median;
            //         if(src[blockindex]<coh_median) --coh_median;
            //     }
            //     for(int iblock=0; iblock<16; ++iblock){
            //         int blockindex=(ichan+iblock)*nsamples+isample;
            //         src[blockindex]=src[blockindex]-coh_median;
            //     }
            // }
            // SAMPLE_TYPE sample=src[index];

            // Instead, let's use a sorting network to find the exact median



#define SWAP(ii, jj) { \
                SAMPLE_TYPE min1=std::min(my_src[ii], my_src[jj]); \
                SAMPLE_TYPE max1=std::max(my_src[ii], my_src[jj]); \
                my_src[ii]=min1; \
                my_src[jj]=max1; \
            }
            if(ichan%16==0){
                // Sorting network works in-place, so we need a copy to work on
                SAMPLE_TYPE my_src[16];
                for(int i=0; i<16; ++i) my_src[i]=src[(ichan+i)*data->nsamples+isample];
                // Swaps for "best" 16-element sorting network, from:
                // http://pages.ripco.net/~jgamble/nw.html
                SWAP(0, 1);
                SWAP(2, 3);
                SWAP(4, 5);
                SWAP(6, 7);
                SWAP(8, 9);
                SWAP(10, 11);
                SWAP(12, 13);
                SWAP(14, 15);
                SWAP(0, 2);
                SWAP(4, 6);
                SWAP(8, 10);
                SWAP(12, 14);
                SWAP(1, 3);
                SWAP(5, 7);
                SWAP(9, 11);
                SWAP(13, 15);
                SWAP(0, 4);
                SWAP(8, 12);
                SWAP(1, 5);
                SWAP(9, 13);
                SWAP(2, 6);
                SWAP(10, 14);
                SWAP(3, 7);
                SWAP(11, 15);
                SWAP(0, 8);
                SWAP(1, 9);
                SWAP(2, 10);
                SWAP(3, 11);
                SWAP(4, 12);
                SWAP(5, 13);
                SWAP(6, 14);
                SWAP(7, 15);
                SWAP(5, 10);
                SWAP(6, 9);
                SWAP(3, 12);
                SWAP(13, 14);
                SWAP(7, 11);
                SWAP(1, 2);
                SWAP(4, 8);
                SWAP(1, 4);
                SWAP(7, 13);
                SWAP(2, 8);
                SWAP(11, 14);
                SWAP(2, 4);
                SWAP(5, 6);
                SWAP(9, 10);
                SWAP(11, 13);
                SWAP(3, 8);
                SWAP(7, 12);
                SWAP(6, 8);
                SWAP(10, 12);
                SWAP(3, 5);
                SWAP(7, 9);
                SWAP(3, 4);
                SWAP(5, 6);
                SWAP(7, 8);
                SWAP(9, 10);
                SWAP(11, 12);
                SWAP(6, 7);
                SWAP(8, 9);

                // Subtract the median
                for(int i=0; i<16; ++i) src[(ichan+i)*data->nsamples+isample] -= my_src[7];
            } // end if(ichan%16==0)

            SAMPLE_TYPE sample=data->src[index];

#else
            if(data->src[index]>median) ++accum;
            if(data->src[index]<median) --accum;
            if(accum>10){
                ++median;
                accum=0;
            }
            if(accum<-10){
                --median;
                accum=0;
            }

            SAMPLE_TYPE sample=data->src[index]-median;
#endif

#ifdef STORE_INTERMEDIATE
            data->pedsub[index]=sample;
#endif
            // --------------------------------------------------------------
            // Filtering
            // --------------------------------------------------------------

            SAMPLE_TYPE filt=0;
            for(int j=0; j<NTAPS; ++j){
                filt+=data->taps[j]*prev_samp[(j+isample)%NTAPS];
            }
            prev_samp[isample%NTAPS]=sample;

#ifdef STORE_INTERMEDIATE
            data->filtered[index]=filt;
#endif

            // --------------------------------------------------------------
            // Hit finding
            // --------------------------------------------------------------
            bool is_over=filt > threshold;
            if(is_over && !prev_was_over) hit_start=isample;
            if(is_over){
                hit_charge+=filt;
                hit_tover++;
                prev_was_over=true;
            }
            if(prev_was_over && !is_over){
                // We reached the end of the hit: write it out
                (*output_loc++) = ichan;
                (*output_loc++) = hit_start;
                (*output_loc++) = hit_charge;
                (*output_loc++) = hit_tover;

                hit_start=0;
                hit_charge=0;
                hit_tover=0;

                ++nhits;
                prev_was_over=false;
            } // end if left hit
            // dst[i]=median;
        } // end loop over samples
    } // end loop over channels
    // printf("Found %d hits\n", nhits);
    // Write a magic "end-of-hits" value into the list of hits
    for(int i=0; i<4; ++i) (*output_loc++) = MAGIC;
}

//======================================================================
void do_processing(TPCData* data)
{
    // Make AVX registers containing the values of the filter taps,
    // which we'll need later
    __m256i tap_256[NTAPS];
    for(int i=0; i<NTAPS; ++i) tap_256[i]= _mm256_set1_epi16(data->taps[i]);

    // Pointer to keep track of where we'll write the next output hit
    __m256i* output_loc=(__m256i*)data->hits;

    // Loop over channels. We go 16 at a time because that's the
    // number of (short int) data points that fit in a 256-bit
    // register
    for(int ichan=data->begin_chan; ichan<data->end_chan-15; ichan+=16){

        // ------------------------------------
        // Variables for pedestal subtraction

        // The current estimate of the pedestal in each channel
        __m256i median;
        // The accumulator that we increase/decrease when the current
        // sample is greater/less than the median
        __m256i accum=_mm256_setzero_si256();

        // ------------------------------------
        // Variables for filtering

        // The (unfiltered) samples `n` places before the current one
        __m256i prev_samp[NTAPS];
        for(int j=0; j<NTAPS; ++j) prev_samp[j]= _mm256_setzero_si256();

        // ------------------------------------
        // Variables for hit finding

        // Was the previous step over threshold?
        __m256i prev_was_over=_mm256_setzero_si256();
        // The time value at which the current hit started
        __m256i hit_start=_mm256_setzero_si256();
        // The integrated charge (so far) of the current hit
        __m256i hit_charge=_mm256_setzero_si256();
        // The time-over-threshold (so far) of the current hit
        __m256i hit_tover=_mm256_setzero_si256();

        // The channel numbers in each of the slots in the register
        __m256i channels=_mm256_set_epi16(ichan + 15,
                                          ichan + 14,
                                          ichan + 13,
                                          ichan + 12,
                                          ichan + 11,
                                          ichan + 10,
                                          ichan +  9,
                                          ichan +  8,
                                          ichan +  7,
                                          ichan +  6,
                                          ichan +  5,
                                          ichan +  4,
                                          ichan +  3,
                                          ichan +  2,
                                          ichan +  1,
                                          ichan +  0);

        // Loop over samples in this block of channels
        for(int isample=0; isample<nsamples; ++isample){

            // --------------------------------------------------------------
            // Pedestal finding
            // --------------------------------------------------------------
#ifdef COHERENT_PEDSUB
            // Try using frugal streaming to estimate the median
            // across channels.  Sadly this just doesn't work: if the
            // starting value is on a big hit, then the median
            // estimation never recovers, and our estimate of the
            // median ends up way too large

            // if(isample%16==0){
            //     __m256i coh_median=_mm256_loadu_si256((__m256i*)(src+ichan*nsamples+isample));
            //     for(int iblock=0; iblock<16; ++iblock){
            //         // This is a loop over the 16 channels that are
            //         // done in one step in the outer loop, with each
            //         // __m256i register holding the next 16 ticks
            //         __m256i s=_mm256_loadu_si256((__m256i*)(src+(ichan+iblock)*nsamples+isample));
            //         // Do the frugal streaming here
            //         __m256i is_gt=_mm256_cmpgt_epi16(s, coh_median);
            //         __m256i is_eq=_mm256_cmpeq_epi16(s, coh_median);

            //         __m256i to_add = _mm256_set1_epi16(-1);
            //         // Really want an epi16 version of this, but the cmpgt and
            //         // cmplt functions set their epi16 parts to 0xff or 0x0,
            //         // so treating everything as epi8 works the same
            //         to_add = _mm256_blendv_epi8(to_add, _mm256_set1_epi16(1), is_gt);
            //         to_add = _mm256_blendv_epi8(to_add, _mm256_set1_epi16(0), is_eq);
            //         coh_median = _mm256_add_epi16(coh_median, to_add);
            //     }
            //     for(int iblock=0; iblock<16; ++iblock){
            //         // This is a loop over the 16 channels that are
            //         // done in one step in the outer loop, with each
            //         // __m256i register holding the next 16 ticks
            //         __m256i s2 = _mm256_loadu_si256((__m256i*)(src+(ichan+iblock)*nsamples+isample));
            //         s2 = _mm256_sub_epi16(s2, coh_median);
            //         _mm256_storeu_si256((__m256i*)(src+(ichan+iblock)*nsamples+isample), s2);
            //     }
            // }

            // // TODO: Make this line less horrible, and maybe also make
            // // its memory access pattern less horrible with cache
            // // blocking?
            // __m256i s=_mm256_set_epi16(src[(ichan + 15)*nsamples+isample],
            //                            src[(ichan + 14)*nsamples+isample],
            //                            src[(ichan + 13)*nsamples+isample],
            //                            src[(ichan + 12)*nsamples+isample],
            //                            src[(ichan + 11)*nsamples+isample],
            //                            src[(ichan + 10)*nsamples+isample],
            //                            src[(ichan +  9)*nsamples+isample],
            //                            src[(ichan +  8)*nsamples+isample],
            //                            src[(ichan +  7)*nsamples+isample],
            //                            src[(ichan +  6)*nsamples+isample],
            //                            src[(ichan +  5)*nsamples+isample],
            //                            src[(ichan +  4)*nsamples+isample],
            //                            src[(ichan +  3)*nsamples+isample],
            //                            src[(ichan +  2)*nsamples+isample],
            //                            src[(ichan +  1)*nsamples+isample],
            //                            src[(ichan +  0)*nsamples+isample]);

#define SWAP256(ii, jj) { \
                __m256i min1=_mm256_min_epi16(my_src[ii], my_src[jj]); \
                __m256i max1=_mm256_max_epi16(my_src[ii], my_src[jj]); \
                my_src[ii]=min1; \
                my_src[jj]=max1; \
            }
            if(isample%16==0){
                __m256i my_src[16];
                for(int i=0; i<16; ++i) my_src[i]=_mm256_loadu_si256((__m256i*)(src+(ichan+i)*data->nsamples+isample));

                // Swaps for "best" 16-element sorting network, from:
                // http://pages.ripco.net/~jgamble/nw.html
                SWAP256(0, 1);
                SWAP256(2, 3);
                SWAP256(4, 5);
                SWAP256(6, 7);
                SWAP256(8, 9);
                SWAP256(10, 11);
                SWAP256(12, 13);
                SWAP256(14, 15);
                SWAP256(0, 2);
                SWAP256(4, 6);
                SWAP256(8, 10);
                SWAP256(12, 14);
                SWAP256(1, 3);
                SWAP256(5, 7);
                SWAP256(9, 11);
                SWAP256(13, 15);
                SWAP256(0, 4);
                SWAP256(8, 12);
                SWAP256(1, 5);
                SWAP256(9, 13);
                SWAP256(2, 6);
                SWAP256(10, 14);
                SWAP256(3, 7);
                SWAP256(11, 15);
                SWAP256(0, 8);
                SWAP256(1, 9);
                SWAP256(2, 10);
                SWAP256(3, 11);
                SWAP256(4, 12);
                SWAP256(5, 13);
                SWAP256(6, 14);
                SWAP256(7, 15);
                SWAP256(5, 10);
                SWAP256(6, 9);
                SWAP256(3, 12);
                SWAP256(13, 14);
                SWAP256(7, 11);
                SWAP256(1, 2);
                SWAP256(4, 8);
                SWAP256(1, 4);
                SWAP256(7, 13);
                SWAP256(2, 8);
                SWAP256(11, 14);
                SWAP256(2, 4);
                SWAP256(5, 6);
                SWAP256(9, 10);
                SWAP256(11, 13);
                SWAP256(3, 8);
                SWAP256(7, 12);
                SWAP256(6, 8);
                SWAP256(10, 12);
                SWAP256(3, 5);
                SWAP256(7, 9);
                SWAP256(3, 4);
                SWAP256(5, 6);
                SWAP256(7, 8);
                SWAP256(9, 10);
                SWAP256(11, 12);
                SWAP256(6, 7);
                SWAP256(8, 9);

                // Subtract the median
                for(int i=0; i<16; ++i){
                    // TODO: This overwrites the input, but we have a
                    // pedsub buffer that's meant for just this
                    // output. Store the output there and retrieve it
                    // below
                    __m256i* addr=(__m256i*)(src+(ichan+i)*data->nsamples+isample);
                    __m256i orig=_mm256_loadu_si256(addr);
                    _mm256_storeu_si256(addr, _mm256_sub_epi16(orig, my_src[7]));
                }
            } // end if(isample%16==0)

            // TODO: Make this line less horrible, and maybe also make
            // its memory access pattern less horrible with cache
            // blocking?
            const int nsamples=data->nsamples;
            __m256i s=_mm256_set_epi16(src[(ichan + 15)*nsamples+isample],
                                       src[(ichan + 14)*nsamples+isample],
                                       src[(ichan + 13)*nsamples+isample],
                                       src[(ichan + 12)*nsamples+isample],
                                       src[(ichan + 11)*nsamples+isample],
                                       src[(ichan + 10)*nsamples+isample],
                                       src[(ichan +  9)*nsamples+isample],
                                       src[(ichan +  8)*nsamples+isample],
                                       src[(ichan +  7)*nsamples+isample],
                                       src[(ichan +  6)*nsamples+isample],
                                       src[(ichan +  5)*nsamples+isample],
                                       src[(ichan +  4)*nsamples+isample],
                                       src[(ichan +  3)*nsamples+isample],
                                       src[(ichan +  2)*nsamples+isample],
                                       src[(ichan +  1)*nsamples+isample],
                                       src[(ichan +  0)*nsamples+isample]);
#else

#ifdef TRANSPOSE_MEMORY
            __m256i s=_mm256_loadu_si256((__m256i*)(src+isample*data->nchannels+ichan));
#else
            // TODO: Make this line less horrible, and maybe also make
            // its memory access pattern less horrible with cache
            // blocking?
            const int nsamples=data->nsamples;
            __m256i s=_mm256_set_epi16(src[(ichan + 15)*nsamples+isample],
                                       src[(ichan + 14)*nsamples+isample],
                                       src[(ichan + 13)*nsamples+isample],
                                       src[(ichan + 12)*nsamples+isample],
                                       src[(ichan + 11)*nsamples+isample],
                                       src[(ichan + 10)*nsamples+isample],
                                       src[(ichan +  9)*nsamples+isample],
                                       src[(ichan +  8)*nsamples+isample],
                                       src[(ichan +  7)*nsamples+isample],
                                       src[(ichan +  6)*nsamples+isample],
                                       src[(ichan +  5)*nsamples+isample],
                                       src[(ichan +  4)*nsamples+isample],
                                       src[(ichan +  3)*nsamples+isample],
                                       src[(ichan +  2)*nsamples+isample],
                                       src[(ichan +  1)*nsamples+isample],
                                       src[(ichan +  0)*nsamples+isample]);
#endif

            if(isample==0) median=s;

            // if the sample is greater than the median, add one to the accumulator
            // if the sample is less than the median, subtract one from the accumulator.

            // For reasons that I don't understand, there's no cmplt
            // for "compare less-than", so we have to compare greater,
            // compare equal, and take everything else to be compared
            // less-then
            __m256i is_gt=_mm256_cmpgt_epi16(s, median);
            __m256i is_eq=_mm256_cmpeq_epi16(s, median);

            __m256i to_add = _mm256_set1_epi16(-1);
            // Really want an epi16 version of this, but the cmpgt and
            // cmplt functions set their epi16 parts to 0xff or 0x0,
            // so treating everything as epi8 works the same
            to_add = _mm256_blendv_epi8(to_add, _mm256_set1_epi16(1), is_gt);
            to_add = _mm256_blendv_epi8(to_add, _mm256_set1_epi16(0), is_eq);

            accum = _mm256_add_epi16(accum, to_add);

            // if(isample<20){
            //     printf("Step 1\n");
            //     printf("vals:  "); print256(s); printf("\n");
            //     printf("med:   "); print256(median); printf("\n");
            //     printf("is_gt: "); print256(is_gt); printf("\n");
            //     printf("is_eq: "); print256(is_eq); printf("\n");
            //     printf("to_add:"); print256(to_add); printf("\n");
            //     printf("acc:   "); print256(accum); printf("\n");
            //     printf("\n");
            // }

            // if the accumulator is >10, add one to the median and
            // set the accumulator to zero. if the accumulator is
            // <-10, subtract one from the median and set the
            // accumulator to zero
            is_gt=_mm256_cmpgt_epi16(accum, _mm256_set1_epi16(10));
            __m256i is_lt=_mm256_cmpgt_epi16(_mm256_sign_epi16(accum, _mm256_set1_epi16(-10)), _mm256_set1_epi16(10));

            to_add = _mm256_setzero_si256();
            to_add = _mm256_blendv_epi8(to_add, _mm256_set1_epi16(1), is_gt);
            to_add = _mm256_blendv_epi8(to_add, _mm256_set1_epi16(-1), is_lt);

            median = _mm256_add_epi16(median, to_add);

            // Reset the channels that were >10 or <-10 to zero, leaving the others unchanged
            accum = _mm256_blendv_epi8(accum, _mm256_setzero_si256(), _mm256_or_si256(is_lt, is_gt));

            // if(isample<20){
            //     printf("Step 2\n");
            //     printf("vals:  "); print256(s); printf("\n");
            //     printf("med:   "); print256(median); printf("\n");
            //     printf("is_gt: "); print256(is_gt); printf("\n");
            //     printf("is_lt: "); print256(is_lt); printf("\n");
            //     printf("to_add:"); print256(to_add); printf("\n");
            //     printf("acc:   "); print256(accum); printf("\n");
            //     printf("med:   "); print256(median); printf("\n");
            //     printf("----------------------------------------------------\n\n");
            // }

            // Actually subtract the pedestal
            s = _mm256_sub_epi16(s, median);

#endif

#ifdef STORE_INTERMEDIATE
            // This stores transposed. Not sure how to fix it without transposing the entire input structure too
            _mm256_storeu_si256((__m256i*)(pedsub+data->nchannels*isample+ichan), s);
#endif

            // --------------------------------------------------------------
            // Filtering
            // --------------------------------------------------------------

            // NB we're doing integer multiplication of
            // (short)*(short) here, so the result can be larger than
            // will fit in a short. `mullo` gives us back the low half
            // of the result, so it will be bogus if the result didn't
            // fit in a short. I'm relying on the fact that we started
            // with a 12-bit number which included the pedestal, and
            // we've subtracted the pedestal and put the result in a
            // 16-bit space to get us never to overflow. Another
            // approach might be to make the filter coeffs so large
            // that we _always_ overflow, and use `mulhi` instead

            // Try pipelining this with multiple accumulators, but it doesn't really make any difference
            __m256i filt0 = _mm256_setzero_si256();
            __m256i filt1 = _mm256_setzero_si256();
            __m256i filt2 = _mm256_setzero_si256();
            __m256i filt3 = _mm256_setzero_si256();

            // % would be slow, but we're making sure that NTAPS is a power of two
#define ADDTAP0(x) filt0 = _mm256_add_epi16(filt0, _mm256_mullo_epi16(tap_256[x], prev_samp[(x+isample)%NTAPS]));
#define ADDTAP1(x) filt1 = _mm256_add_epi16(filt1, _mm256_mullo_epi16(tap_256[x], prev_samp[(x+isample)%NTAPS]));
#define ADDTAP2(x) filt2 = _mm256_add_epi16(filt2, _mm256_mullo_epi16(tap_256[x], prev_samp[(x+isample)%NTAPS]));
#define ADDTAP3(x) filt3 = _mm256_add_epi16(filt3, _mm256_mullo_epi16(tap_256[x], prev_samp[(x+isample)%NTAPS]));

            ADDTAP0(0);
            ADDTAP1(1);
            ADDTAP2(2);
            ADDTAP3(3);
            ADDTAP0(4);
            ADDTAP1(5);
            ADDTAP2(6);
#if NTAPS==32
            ADDTAP3(7);
            ADDTAP0(8);
            ADDTAP1(9);
            ADDTAP2(10);
            ADDTAP3(11);
            ADDTAP0(12);
            ADDTAP1(13);
            ADDTAP2(14);
            ADDTAP3(15);
            ADDTAP0(16);
            ADDTAP1(17);
            ADDTAP2(18);
            ADDTAP3(19);
            ADDTAP0(20);
            ADDTAP1(21);
            ADDTAP2(22);
            ADDTAP3(23);
            ADDTAP0(24);
            ADDTAP1(25);
            ADDTAP2(26);
            ADDTAP3(27);
            ADDTAP0(28);
            ADDTAP1(29);
            ADDTAP2(30);
            ADDTAP3(31);
#endif
            __m256i filt = _mm256_add_epi16(_mm256_add_epi16(filt0, filt1), _mm256_add_epi16(filt2, filt3));
            prev_samp[isample%NTAPS]=s;

#ifdef STORE_INTERMEDIATE
            // This stores transposed. Not sure how to fix it without transposing the entire input structure too
            _mm256_storeu_si256((__m256i*)(filtered+data->nchannels*isample+ichan), filt);
#endif

            // --------------------------------------------------------------
            // Hit finding
            // --------------------------------------------------------------
            // Mask for channels that are over the threshold in this step
            __m256i is_over=_mm256_cmpgt_epi16(filt, _mm256_set1_epi16(threshold));
            // Mask for channels that entered "over threshold" state
            // this step. andnot(a,b)=(NOT a AND b)
            __m256i entered=_mm256_andnot_si256(prev_was_over, is_over);
            // Mask for channels that left "over threshold" state this step
            __m256i left=_mm256_andnot_si256(is_over, prev_was_over);

            //-----------------------------------------
            // Update hit start times for the channels where a hit started
            hit_start=_mm256_blendv_epi8(hit_start, _mm256_set1_epi16(isample), entered);

            //-----------------------------------------
            // Accumulate charge and time-over-threshold in the is_over channels

            // Really want an epi16 version of this, but the cmpgt and
            // cmplt functions set their epi16 parts to 0xff or 0x0,
            // so treating everything as epi8 works the same
            __m256i to_add_charge=_mm256_blendv_epi8(_mm256_set1_epi16(0), filt, is_over);
            hit_charge=_mm256_add_epi16(hit_charge, to_add_charge);

            __m256i to_add_tover=_mm256_blendv_epi8(_mm256_set1_epi16(0), _mm256_set1_epi16(1), is_over);
            hit_tover=_mm256_add_epi16(hit_tover, to_add_tover);


            // if(isample>2700 && isample<3713){
            //     printf("Sample %d\n", isample);
            //     printf("vals:      "); print256(s); printf("\n");
            //     printf("is_over:   "); print256(is_over); printf("\n");
            //     printf("entered:   "); print256(entered); printf("\n");
            //     printf("left:      "); print256(left); printf("\n");
            //     printf("to_add_ch: "); print256(to_add_charge); printf("\n");
            //     printf("to_add_t:  "); print256(to_add_tover); printf("\n");
            //     printf("hit_charge:"); print256(hit_charge); printf("\n");
            //     printf("hit_tover: "); print256(hit_tover); printf("\n");

            //     printf("\n");
            // }

            // Only store the values if there are >0 hits ending on
            // this sample. We have to save the entire 16-channel
            // register, which is inefficient, but whatever

            // Testing whether a whole register is zeroes turns out to be tricky. Here's a way:
            //
            // https://stackoverflow.com/questions/22674205/is-there-an-or-equivalent-to-ptest-in-x64-assembly
            //
            // In x64 assembly, PTEST %XMM0 -> %XMM1 sets the
            // zero-flag if none of the same bits are set in %XMM0 and
            // %XMM1, and sets the carry-flag if everything that is
            // set in %XMM0 is also set in %XMM1:
            int no_hits_to_store=_mm256_testc_si256(_mm256_setzero_si256(), left);
            // printf("no_hits: %d\n", no_hits_to_store);
            if(!no_hits_to_store){

                // We have to save the whole register, including the
                // lanes that don't have anything interesting, but
                // we'll mask them to zero so they're easy to remove
                // in a later processing step. (TODO: Maybe we should
                // do that processing step in this function?)
#define STORE_MASK(x) _mm256_storeu_si256(output_loc++, _mm256_blendv_epi8(_mm256_set1_epi16(0), x, left));
                STORE_MASK(channels);
                STORE_MASK(hit_start);
                STORE_MASK(hit_charge);
                STORE_MASK(hit_tover);

                // reset hit_start, hit_charge and hit_tover in the channels we saved
                hit_start=_mm256_blendv_epi8(hit_start, _mm256_set1_epi16(0), left);
                hit_charge=_mm256_blendv_epi8(hit_charge, _mm256_set1_epi16(0), left);
                hit_tover=_mm256_blendv_epi8(hit_tover, _mm256_set1_epi16(0), left);
            }

            prev_was_over=is_over;
        }
    }
    for(int i=0; i<4; ++i) _mm256_storeu_si256(output_loc++, _mm256_set1_epi16(MAGIC));
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

    std::string inputfile=vm["input"].as<std::string>();
    const int max_channels=-1;
    // How many times to repeat the data in memory so as to make sure
    // it's bigger than the cache
    const int nrepeat=16;

    const bool transpose=false;
    TPCData data(inputfile, nrepeat, transpose);
    data->setTaps(NTAPS, 100);
    //----------------------------------------------------------------
    printf("Using %d samples, %d collection channels (%d unique) = %.1f APA*ms with type size %ld bytes. Total size %.1f MB\n",
           data.nsamples, data.nchannels, data.nchannels_uniq, data.APAmsData(), sizeof(SAMPLE_TYPE),
           data.dataSizeMB());

    const int first_chan=0;
    const int last_chan=data.nchannels;

    // ---------------------------------------------------------
    // Run for benchmarking
    // ---------------------------------------------------------

    // printf("% 20s     Threads\n", "");
    // printf("% 20s ", "");
    // for(int i=0; i<nnthreads; ++i){
    //     printf("% 8d ", nthreads[i]);
    // }
    // printf("\n");

    // timefnNThreads("intrin", do_processing,
    //                src,
    //                nchannels, nsamples,
    //                first_chan, last_chan,
    //                taps, ntaps,
    //                pedsub, filtered);

    // timefnNThreads("naive", do_processing_naive,
    //                src,
    //                nchannels, nsamples,
    //                first_chan, last_chan,
    //                taps, ntaps,
    //                pedsub, filtered);

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
    // // The processing overwrites the input, so give each run its own copy
    // SAMPLE_TYPE* my_src      = (SAMPLE_TYPE*)malloc(nsize*sizeof(SAMPLE_TYPE));
    // memcpy(my_src, src, nsize*sizeof(SAMPLE_TYPE));

    do_processing_naive(&data);

    // saveNaiveHitsToFile(hits, "hits-" TAG "-naive");
#ifdef STORE_INTERMEDIATE
    // saveOutputToFile(pedsub, nchannels_uniq, nsamples, "pedsub-" TAG "-naive", false);
    // saveOutputToFile(filtered, nchannels_uniq, nsamples, "filtered-" TAG "-naive", false);
#endif

    // memcpy(my_src, src, nsize*sizeof(SAMPLE_TYPE));
    do_processing(&data);
    // saveIntrinHitsToFile(hits, "hits-" TAG "-intrin");
#ifdef STORE_INTERMEDIATE
    // saveOutputToFile(pedsub, nchannels_uniq, nsamples, "pedsub-" TAG "-intrin", true);
    // saveOutputToFile(filtered, nchannels_uniq, nsamples, "filtered-" TAG "-intrin", true);
#endif
}
