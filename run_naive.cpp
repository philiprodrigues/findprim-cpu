#include <map>
#include <set>
#include <vector>
#include <algorithm> // For set_union

#include "ExpandedCollectionFile.hh"

#include "pthread.h"

#include "constants.hh"
#include "design_fir.h"
#include "process_naive.h"
#include "process_avx2.h"

struct Hit
{
    uint16_t chan, hit_end, hit_charge, hit_tover;
    bool operator<(const Hit& rhs) const {
        if(chan!=rhs.chan)             return chan<rhs.chan;
        if(hit_end!=rhs.hit_end)   return hit_end<rhs.hit_end;
        if(hit_charge!=rhs.hit_charge) return hit_charge<rhs.hit_charge;
        return hit_tover<rhs.hit_tover;
    }

    bool operator==(const Hit& rhs) const {
        return chan==rhs.chan           &&
            hit_end==rhs.hit_end    &&
            hit_charge==rhs.hit_charge  &&
            hit_tover==rhs.hit_tover;
    }
};

//======================================================================
std::set<Hit> get_naive_hits(const uint16_t* input_loc)
{
    std::set<Hit> ret;

    while(*input_loc!=MAGIC){
        uint16_t chan=*input_loc++;
        uint16_t hit_end=*input_loc++;
        uint16_t hit_charge=*input_loc++;
        uint16_t hit_tover=*input_loc++;
        ret.insert(Hit{chan, hit_end, hit_charge, hit_tover});
    }
    return ret;
}

//======================================================================
int main(int, char**)
{
    pthread_setname_np(pthread_self(), "main");

    // mf::setStandAloneMessageThreshold({"INFO"});

    ExpandedCollectionFile f("/data/lar/dunedaq/rodrigues/felix-long-readout/dump-link-raw-2019-06-24/expanded-frame-dump-link0-4000000-ticks-timestamp-0x1155bcddf21d100.dat");
    // char* fragment=reinterpret_cast<char*>(f.fragment(0));
    // size_t length=FrameFile::frames_per_fragment*sizeof(dune::FelixFrame);
    size_t n_messages=f.num_messages();

    const uint8_t tap_exponent=6;
    const int multiplier=1<<tap_exponent; // 64
    std::vector<int16_t> taps=firwin_int(7, 0.1, multiplier);
    taps.push_back(0); // Make it 8 long so it's a power of two

    // TODO: ProcessingInfo points to but doesn't copy the taps it's
    // passed, so we make a long-lived copy here. We ought to not leak
    // this (and also do it properly)
    int16_t* taps_p=new int16_t[taps.size()];
    for(size_t i=0; i<taps.size(); ++i) taps_p[i]=taps[i];

    // Temporary place to stash the hits
    uint16_t* primfind_dest_naive=new uint16_t[100000];
    
    ProcessingInfo pi_naive(nullptr,
                            FRAMES_PER_MSG, // We'll just process one message
                            0, // First register
                            REGISTERS_PER_FRAME, // Last register
                            primfind_dest_naive,
                            taps_p, (uint8_t)taps.size(),
                            tap_exponent,
                            0, // nhits
                            0); // absTimeModNTAPS

    size_t nhits_naive=0;
    for(size_t imessage=0; imessage<std::min(n_messages, 500UL); ++imessage){
        // printf("imessage=%d, info.absTimeModNTAPS=%d\n", imessage, pi_naive.absTimeModNTAPS);
        // SUPERCHUNK_CHAR_STRUCT* scs=reinterpret_cast<SUPERCHUNK_CHAR_STRUCT*>(fragment+imessage*NETIO_MSG_SIZE);
        // RegisterArray<REGISTERS_PER_FRAME*FRAMES_PER_MSG> expanded=expand_message_adcs(*scs);
        MessageCollectionADCs* mcadc=f.message(imessage);
        if(imessage==0){
            pi_naive.setState(mcadc);
        }
        pi_naive.input=mcadc;
        process_window_naive(pi_naive);
        std::set<Hit> naive_hits=get_naive_hits(primfind_dest_naive);
        nhits_naive+=naive_hits.size();
        for(auto const& hit: naive_hits){
            printf("% 5d % 5d % 5d % 5d\n", hit.chan, imessage*FRAMES_PER_MSG+hit.hit_end-hit.hit_tover, hit.hit_charge, hit.hit_tover);
        }
    }

    delete[] primfind_dest_naive;

}

/* Local Variables:  */
/* mode: c++         */
/* c-basic-offset: 4 */
/* End:              */
