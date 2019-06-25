#pragma once

#include <limits>
#include <cstdint>

#include "immintrin.h"

const unsigned short MAGIC = std::numeric_limits<unsigned short>::max();
const int16_t threshold=2000;

// How many frames are concatenated in one netio message
static constexpr size_t FRAMES_PER_MSG=12;
// How many collection-wire AVX2 registers are returned per
// frame.
static constexpr size_t REGISTERS_PER_FRAME=6;
// How many bytes are in an AVX2 register
static constexpr size_t BYTES_PER_REGISTER=32;
// How many samples are in a register
static constexpr size_t SAMPLES_PER_REGISTER=16;

static const size_t FELIX_FRAME_SIZE=464;
static const size_t NETIO_MSG_SIZE=FELIX_FRAME_SIZE*FRAMES_PER_MSG;

// One netio message's worth of collection channel ADCs after
// expansion: 12 frames per message times 8 registers per frame times
// 32 bytes (256 bits) per register
static const size_t collection_adcs_size=BYTES_PER_REGISTER*REGISTERS_PER_FRAME*FRAMES_PER_MSG;
struct MessageCollectionADCs {
    char fragments[collection_adcs_size];
};

struct NETIO_MSG_STRUCT {
    char fragments[NETIO_MSG_SIZE];
};

struct WindowCollectionADCs {
    WindowCollectionADCs(size_t numMessages_, MessageCollectionADCs* fragments_)
        : numMessages(numMessages_),
          fragments(fragments_)
    {}

    // Get a pointer to register `ireg` at time `itime`, as an AVX2 int register
    const __m256i* get256(size_t ireg, size_t itime) const
    {
            const size_t msg_index=itime/12;
            const size_t msg_time_offset=itime%12;
            const size_t index=msg_index*REGISTERS_PER_FRAME*FRAMES_PER_MSG+FRAMES_PER_MSG*ireg+msg_time_offset;
            const __m256i* rawp=reinterpret_cast<const __m256i*>(fragments)+index;
            return rawp;
    }

    uint16_t get16(size_t ichan, size_t itime) const
    {
        const size_t register_index=ichan/SAMPLES_PER_REGISTER;
        const size_t register_offset=ichan%SAMPLES_PER_REGISTER;
        const size_t register_t0_start=register_index*SAMPLES_PER_REGISTER*FRAMES_PER_MSG;
        const size_t msg_index=itime/12;
        const size_t msg_time_offset=itime%12;
        // The index in uint16_t of the start of the message we want
        const size_t msg_start_index=msg_index*sizeof(MessageCollectionADCs)/sizeof(uint16_t);
        const size_t offset_within_msg=register_t0_start+SAMPLES_PER_REGISTER*msg_time_offset+register_offset;
        const size_t index=msg_start_index+offset_within_msg;
        return *(reinterpret_cast<uint16_t*>(fragments)+index);
    }

    size_t numMessages;
    MessageCollectionADCs* __restrict__ fragments;
};

// A little wrapper around an array of 256-bit registers, so that we
// can explicitly access it as an array of 256-bit registers or as an
// array of uint16_t
template<int N>
class RegisterArray
{
public:
    // Get the value at the ith position as a 256-bit register
    __m256i ymm(size_t i) { return _mm256_lddqu_si256(reinterpret_cast<__m256i*>(m_array)+i); }
    void set_ymm(size_t i, __m256i val) { _mm256_storeu_si256(reinterpret_cast<__m256i*>(m_array)+i, val); }

    uint16_t uint16(size_t i) { return m_array[i]; }
    void set_uint16(size_t i, uint16_t val) { m_array[i]=val; }

    // Access the jth entry in the ith register
    uint16_t uint16(size_t i, size_t j) { return m_array[16*i+j]; }
    void set_uint16(size_t i, size_t j, uint16_t val) { m_array[16*i+j]=val; }

    uint16_t* data() { return m_array; }
    const uint16_t* data() const { return m_array; }
private:
    alignas(32) uint16_t __restrict__ m_array[N*16];
};
