#ifndef RUN_ALGS_NAIVE_H
#define RUN_ALGS_NAIVE_H

#include "constants.h"

//======================================================================
template<bool TRANSPOSE, bool COHERENT_PEDSUB, bool STORE_INTERMEDIATE, int NTAPS>
void run_algs_naive(TPCData* data, int begin_chan, int end_chan, int hit_offset)
{
    unsigned short* output_loc=data->hits+hit_offset;
    int nhits=0;
    
    for(int ichan=begin_chan; ichan < end_chan; ++ichan){
        SAMPLE_TYPE median=TRANSPOSE ? data->src[ichan] : data->src[ichan*data->nsamples];
        SAMPLE_TYPE accum=0;
        
        // Variables for filtering
        SAMPLE_TYPE prev_samp[NTAPS]={0};

        // Variables for hit finding
        bool prev_was_over=false; // was the previous sample over threshold?
        unsigned short hit_start=0; // start time of the hit
        unsigned short hit_charge=0;
        unsigned short hit_tover=0; // time over threshold
        
        for(int isample=0; isample<data->nsamples; ++isample){
            
            const size_t index=TRANSPOSE ? isample*data->nchannels+ichan : ichan*data->nsamples+isample;
            
            // --------------------------------------------------------------
            // Pedestal finding/coherent noise removal
            // --------------------------------------------------------------
            SAMPLE_TYPE sample;
            
            if(COHERENT_PEDSUB){
                // Use a sorting network to find the exact median
#define SWAP(ii, jj) {                                                  \
                    SAMPLE_TYPE min1=std::min(my_src[ii], my_src[jj]);  \
                    SAMPLE_TYPE max1=std::max(my_src[ii], my_src[jj]);  \
                    my_src[ii]=min1;                                    \
                    my_src[jj]=max1;                                    \
                }
                if(ichan%16==0){
                    // Sorting network works in-place, so we need a copy to work on
                    SAMPLE_TYPE my_src[16];
                    for(int i=0; i<16; ++i) my_src[i]=data->src[(ichan+i)*data->nsamples+isample];
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
                    for(int i=0; i<16; ++i) data->src[(ichan+i)*data->nsamples+isample] -= my_src[7];
                } // end if(ichan%16==0)
                
                sample=data->src[index];
            } // end if coherent_pedsub
            else{
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
                
                sample=data->src[index]-median;
            }
            
            if(STORE_INTERMEDIATE) data->pedsub[index]=sample;
            
            // --------------------------------------------------------------
            // Filtering
            // --------------------------------------------------------------
            
            SAMPLE_TYPE filt=0;
            for(int j=0; j<NTAPS; ++j){
                filt+=data->taps[j]*prev_samp[(j+isample)%NTAPS];
            }
            prev_samp[isample%NTAPS]=sample;
            
            if(STORE_INTERMEDIATE) data->filtered[index]=filt;
            
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

// Make some necessary template specializations
//                                  TRANSP   COH  STORE TAPS
template<> void run_algs_naive<false, false, true,    8>(TPCData* data, int begin_chan, int end_chan, int hit_offset);
template<> void run_algs_naive<true,  false, true,    8>(TPCData* data, int begin_chan, int end_chan, int hit_offset);
template<> void run_algs_naive<false, false, false,  32>(TPCData* data, int begin_chan, int end_chan, int hit_offset);

#endif
