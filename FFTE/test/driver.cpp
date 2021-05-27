#include "test.hpp"
#include "Engine.hpp"
#include <iostream>

/**
 * Code Author: Danny Sharp
 * This file is part of the implementation for a stock FFT algorithm intended for HeFFTe
 */
using namespace stock_fft;
using F = float;
// Main function. Used for calling different parts of the testing code
int main() {
    const int P = 5;
    const int Q = 5;
    auto input = std::array<std::array<std::complex<double>,Q>,P>();
    for(int p = 0; p < P; p++) {
        for(int q = 0; q < Q; q++) {
            input[p][q] = std::complex<double>(p, q);
            auto eol = (q < (Q-1)) ? ", " : ";\n";
            std::cout << input[p][q] << eol;
        }
    }

    auto output = STOCK_FFT::dim2engine<P,Q>::fft2(input, STOCK_FFT::Direction::forward, STOCK_FFT::Major::row);
    for(int p = 0; p < P; p++) {
        for(int q = 0; q < Q; q++) {
            auto eol = (q < (Q-1)) ? ", " : ";\n";
            std::cout << output[p][q] << eol;
        }
    }
}