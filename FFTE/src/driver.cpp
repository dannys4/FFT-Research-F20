#include "test.hpp"
#include "Engine.hpp"
#include <iostream>

/**
 * Code Author: Danny Sharp
 * This file is part of FFTE (Fast Fourier Transform Engine)
 */
using namespace FFTE;
// Main function. Used for calling different parts of the testing code
int main() {
    const int n = 4*9*7;
    auto input  = FFTE::complex_alloc<double,4>(n);
    auto output = FFTE::complex_alloc<double,4>(n);
    for(int i = 0; i < n; i++) {
        std::complex<double> tmp1{1.*i, 2.*i};
        std::complex<double> tmp2{3.*i, 4.*i};
        std::complex<double> tmp[2] {tmp1, tmp2};
        input[i] = FFTE::Complex<double,4>(tmp);
    }
    std::cout << "Output:\n";
    FFTE::fft<double, 4, n>(input, output, FFTE::Direction::inverse);
    for(int i = 0; i < n; i++) {
        std::cout << output[i]/((double) n) << "\n";
    }
    _mm_free(input);
    _mm_free(output);
}