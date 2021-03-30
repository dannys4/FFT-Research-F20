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
    time_batch_fft();
}

    /* const int n = 2*3*5;
    const int L = 8;
    FFTE::complex_vector<F, L> input {};
    for(int i = 0; i < n; i++) {
        input.push_back(Complex<F, L> {1.f*i, 2.f*i, 3.f*i, 4.f*i, 5.f*i, 6.f*i, 7.f*i, 8.f*i});
    }
    std::cout << "Output:\n";
    auto output = engine<F, L>::fft(input, FFTE::Direction::inverse);
    for(int i = 0; i < n; i++) {
        std::cout  << output[i]/((F) n) << " \n";
    } */