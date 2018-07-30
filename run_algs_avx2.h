#ifndef RUN_ALGS_INTRINSIC_H
#define RUN_ALGS_INTRINSIC_H

#include "constants.h"
#include <immintrin.h>

//======================================================================
template<bool TRANSPOSE, bool COHERENT_PEDSUB, bool STORE_INTERMEDIATE, int NTAPS>
void run_algs_avx2(TPCData* data, int begin_chan, int end_chan, int hit_offset)
{
    // Make AVX registers containing the values of the filter taps,
    // which we'll need later
    __m256i tap_256[NTAPS];
    for(int i=0; i<NTAPS; ++i) tap_256[i]= _mm256_set1_epi16(data->taps[i]);
        
    // Pointer to keep track of where we'll write the next output hit
    __m256i* output_loc=(__m256i*)(data->hits+hit_offset);
        
    // Loop over channels. We go 16 at a time because that's the
    // number of (short int) data points that fit in a 256-bit
    // register
    for(int ichan=begin_chan; ichan<end_chan-15; ichan+=16){

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
        for(int isample=0; isample<data->nsamples; ++isample){
            __m256i s;
            // --------------------------------------------------------------
            // Pedestal finding
            // --------------------------------------------------------------
            if(COHERENT_PEDSUB){
#define SWAP256(ii, jj) {                                               \
                    __m256i min1=_mm256_min_epi16(my_src[ii], my_src[jj]); \
                    __m256i max1=_mm256_max_epi16(my_src[ii], my_src[jj]); \
                    my_src[ii]=min1;                                    \
                    my_src[jj]=max1;                                    \
                }
                if(isample%16==0){
                    __m256i my_src[16];
                    for(int i=0; i<16; ++i) my_src[i]=_mm256_loadu_si256((__m256i*)(data->src+(ichan+i)*data->nsamples+isample));

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
                        __m256i* addr=(__m256i*)(data->pedsub+(ichan+i)*data->nsamples+isample);
                        __m256i orig=_mm256_loadu_si256(addr);
                        _mm256_storeu_si256(addr, _mm256_sub_epi16(orig, my_src[7]));
                    }
                } // end if(isample%16==0)

                // TODO: Make this line less horrible, and maybe also make
                // its memory access pattern less horrible with cache
                // blocking?
                const int nsamples=data->nsamples;
                const SAMPLE_TYPE* pedsub=data->pedsub;
                s=_mm256_set_epi16(pedsub[(ichan + 15)*nsamples+isample],
                                   pedsub[(ichan + 14)*nsamples+isample],
                                   pedsub[(ichan + 13)*nsamples+isample],
                                   pedsub[(ichan + 12)*nsamples+isample],
                                   pedsub[(ichan + 11)*nsamples+isample],
                                   pedsub[(ichan + 10)*nsamples+isample],
                                   pedsub[(ichan +  9)*nsamples+isample],
                                   pedsub[(ichan +  8)*nsamples+isample],
                                   pedsub[(ichan +  7)*nsamples+isample],
                                   pedsub[(ichan +  6)*nsamples+isample],
                                   pedsub[(ichan +  5)*nsamples+isample],
                                   pedsub[(ichan +  4)*nsamples+isample],
                                   pedsub[(ichan +  3)*nsamples+isample],
                                   pedsub[(ichan +  2)*nsamples+isample],
                                   pedsub[(ichan +  1)*nsamples+isample],
                                   pedsub[(ichan +  0)*nsamples+isample]);
            } // end if(COHERENT_PEDSUB)
            else{
                if(TRANSPOSE){
                    s=_mm256_loadu_si256((__m256i*)(data->src+isample*data->nchannels+ichan));
                }
                else{
                    // TODO: Make this line less horrible, and maybe also make
                    // its memory access pattern less horrible with cache
                    // blocking?
                    const int nsamples=data->nsamples;
                    const SAMPLE_TYPE* src=data->src;
                    s=_mm256_set_epi16(src[(ichan + 15)*nsamples+isample],
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
                }
                
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

            } // end if(coherent_pedsub) else block

            if(STORE_INTERMEDIATE){
                // This stores transposed. Not sure how to fix it without transposing the entire input structure too
                _mm256_storeu_si256((__m256i*)(data->pedsub+data->nchannels*isample+ichan), s);
            }

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
            if(NTAPS==32){
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
            }

            __m256i filt = _mm256_add_epi16(_mm256_add_epi16(filt0, filt1), _mm256_add_epi16(filt2, filt3));
            prev_samp[isample%NTAPS]=s;
            
            if(STORE_INTERMEDIATE){
                // This stores transposed. Not sure how to fix it without transposing the entire input structure too
                _mm256_storeu_si256((__m256i*)(data->filtered+data->nchannels*isample+ichan), filt);
            }
            
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

// Make some necessary template specializations
//                                  TRANSP   COH  STORE TAPS
template<> void run_algs_avx2<false, false, true,    8>(TPCData* data, int begin_chan, int end_chan, int hit_offset);
template<> void run_algs_avx2<true,  false, true,    8>(TPCData* data, int begin_chan, int end_chan, int hit_offset);
template<> void run_algs_avx2<false, false, false,  32>(TPCData* data, int begin_chan, int end_chan, int hit_offset);

#endif
