#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <limits>

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

#endif
