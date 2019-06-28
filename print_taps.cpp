#include "design_fir.h"
#include <iostream>

int main(int, char**)
{
    const uint8_t tap_exponent=6;
    const int multiplier=1<<tap_exponent; // 64
    std::vector<int16_t> taps=firwin_int(7, 0.1, multiplier);
    taps.push_back(0); // Make it 8 long so it's a power of two
    for(auto tap: taps) std::cout << tap << " ";
    std::cout << std::endl;
}
