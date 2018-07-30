#ifndef PROCESS_SAMPLES_H
#define PROCESS_SAMPLES_H

#include "../waveform-tools/read_samples.h"
#include "design_fir.h"

#include <cstring> // for memset

const int nchannels_per_apa=2560;
const int ncollection_per_apa=960;
const float sampling_rate=2e6; // 2 MHz sampling rate of each channe

struct TPCData
{
    TPCData(const char* inputfile, int nrepeat, bool transpose)
    {
        Waveforms<SAMPLE_TYPE> v=read_samples_text<SAMPLE_TYPE>(inputfile, -1);
        // std::vector<int>& channels=v.channels;
        std::vector<std::vector<SAMPLE_TYPE> >& samples_vector=v.samples;
        
        nchannels_uniq=samples_vector.size();
        // Round down to the nearest multiple of 16 so we don't have to
        // worry about tail effects
        nchannels=((nchannels_uniq*nrepeat)/16)*16;
        nsamples=samples_vector[0].size(); // read_samples checks that all the channels have the same #samples
        nsize_uniq = nsamples*nchannels_uniq;
        nsize    = nsize_uniq*nrepeat;
        
        // The source data
        src      = (SAMPLE_TYPE*)malloc(nsize*sizeof(SAMPLE_TYPE));
        // Intermediate steps' data
        pedsub   = (SAMPLE_TYPE*)malloc(nsize*sizeof(SAMPLE_TYPE));
        filtered = (SAMPLE_TYPE*)malloc(nsize*sizeof(SAMPLE_TYPE));
        
        // The list of output hits
        hits=(unsigned short*)malloc(nsize*sizeof(unsigned short));

        if(transpose){
            // TODO this needs fixing for nrepeat!=1. Copy from avx512 branch
            for(int isample=0; isample<nsamples; ++isample){
                for(int ichan=0; ichan<nchannels; ++ichan){
                    const int index=isample*nchannels+ichan;
                    src[index]=samples_vector[ichan][isample];
                }
            }
        }
        else{
            for(int irep=0; irep<nrepeat; ++irep){
                for(int ichan=0; ichan<nchannels_uniq; ++ichan){
                    for(int isample=0; isample<nsamples; ++isample){
                        const int index=irep*nsize_uniq + ichan*nsamples + isample;
                        src[index]=samples_vector[ichan][isample];
                    }
                }
            }
        }

    }

    void setTaps(int ntaps, int multiplier)
    {
        //----------------------------------------------------------------
        // Create the filter taps

        // Some evil going on here: I want NTAPS to be a preprocessor
        // variable so I can build versions of the code with different
        // numbers of taps. I also want it to be a power of two for modulo
        // reasons in do_processing. But (I think) the number of taps for
        // my lowpass filter has to be odd or I get funny behaviour. So
        // make a FIR filter with one less tap than NTAPS, and append a
        // zero to make it the NTAPS long again
        ntaps=NTAPS;
        std::vector<double> coeffs_double(firwin(ntaps-1, 0.1));
        coeffs_double.push_back(0);

        // The coefficients of the FIR filter are floating point #s <1, so
        // if we want the filtering to do anything sensible in 'short'
        // mode, we have to multiply them up by something
        multiplier=100;

        printf("FIR coeffs: ");
        for(int i=0; i<ntaps; ++i){
            printf("%.3f ", coeffs_double[i]);
            taps[i]=(SAMPLE_TYPE)multiplier*coeffs_double[i];
        }
        printf("\n");

    }

    ~TPCData()
    {
        free(src);
        free(pedsub);
        free(filtered);
        free(hits);
    }

    void zeroOutput()
    {
        memset(pedsub, 0, nsize*sizeof(SAMPLE_TYPE));
        memset(filtered, 0, nsize*sizeof(SAMPLE_TYPE));
        memset(hits, 0, nsize*sizeof(unsigned short));
    }

    float msData() { return nsamples*sampling_rate*1000; }
    float APAmsData() { return msData()*nchannels/ncollection_per_apa; }
    float dataSizeMB() { return float(nchannels*nsamples)*sizeof(SAMPLE_TYPE)/(1024*1024); }

    int nchannels;
    int nsamples;

    int nchannels_uniq;
    int nsamples_uniq;

    int nsize;
    int nsize_uniq;

    int begin_chan;
    int end_chan;

    SAMPLE_TYPE * __restrict__ taps;
    int ntaps;

    SAMPLE_TYPE* __restrict__ src;
    unsigned short*  __restrict__  hits;
    SAMPLE_TYPE * __restrict__ pedsub;
    SAMPLE_TYPE * __restrict__ filtered;
};

typedef void (*processing_fn)(TPCData*);


#endif
