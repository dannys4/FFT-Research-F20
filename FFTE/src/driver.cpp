#include "test.hpp"
#include "Engine.hpp"
#include <iostream>

/**
 * Code Author: Danny Sharp
 * This file is part of FFTE (Fast Fourier Transform Engine)
 */
using namespace FFTE;
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

    auto output = FFTE::fft2(input, FFTE::Direction::forward);
    for(int p = 0; p < P; p++) {
        for(int q = 0; q < Q; q++) {
            auto eol = (q < (Q-1)) ? ", " : ";\n";
            std::cout << output[p][q] << eol;
        }
    }
}
