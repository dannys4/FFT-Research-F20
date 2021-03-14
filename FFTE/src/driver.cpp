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
    const int n = 2*3*5;
    FFTE::complex_vector<F, 4> input {};
    for(int i = 0; i < n; i++) {
        input.push_back(Complex<F, 4> {1.f*i, 2.f*i, 3.f*i, 4.f*i});
    }
    std::cout << "Output:\n";
    auto output = engine<F, 4>::fft(input, FFTE::Direction::inverse);
    for(int i = 0; i < n; i++) {
        std::cout  << output[i]/((F) n) << "\n";
    }
    // check_complex_ops();
}