#include "test.hpp"
#include "new_complex.hpp"

/**
 * Code Author: Danny Sharp
 * This file is part of FFTE (Fast Fourier Transform Engine)
 */

// Main function. Used for calling different parts of the testing code
int main() {
    // check_fft(FFTE::Direction::forward);
    double d4[] {1., 2., 3., 4.};
    double d2[] {1., 2.};
    float  s4[] {1., 2., 3., 4.};

    FFTE::NewComplex<double, 4> xd4 (d4);
    FFTE::NewComplex<double, 2> xd2 (d2);
    FFTE::NewComplex<float, 4>  xs4 (s4);

    // [1+2i, 3+4i] .* [1+2i, 3+4i]
    auto nn = xd4*xd4;
    nn.getComplex(d4);
    std::cout << "d4 = " << d4[0] << " + " << d4[1] << "i, " << d4[2] << " + " << d4[3] << "i\n";
}