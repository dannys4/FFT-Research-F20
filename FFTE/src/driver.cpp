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
    const int n = 2*3*5;
    FFTE::complex_vector<double, 4> input {};
    for(int i = 0; i < n; i++) {
        input.push_back(Complex<double, 4> {1.*i, 2.*i, 3.*i, 4.*i});
    }
    std::cout << "Output:\n";
    auto output = engine<double, 4>::fft(input, FFTE::Direction::inverse);
    for(int i = 0; i < n; i++) {
        std::cout << input[i] << "\t\t" << output[i]/((double) n) << "\n";
    }
}