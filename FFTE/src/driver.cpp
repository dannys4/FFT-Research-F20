#include "test.hpp"
#include "Engine.hpp"
#include <iostream>
#include <vector>

/**
 * Code Author: Danny Sharp
 * This file is part of FFTE (Fast Fourier Transform Engine)
 */
using namespace FFTE;
// Main function. Used for calling different parts of the testing code
int main() {
    const int n = 4*9*7;
    FFTE::complex_vector<double, 4> input ();
    // FFTE::Complex<double, 4>* input2 = new FFTE::Complex<double, 4>[20];
    for(int i = 0; i < n; i++) {
        Complex<double, 4> tmp {1.*i, 2.*i, 3.*i, 4.*i};
        input.push_back(tmp);
        std::cout << "tmp = " << tmp << ", input[i] = " << input[i] << "\n";
    }
    std::cout << "Output:\n";
    auto output = fft<double, 4>(input, FFTE::Direction::forward);
    // for(int i = 0; i < n; i++) {
    //     std::cout << output[i] << "\n";
    // }
}