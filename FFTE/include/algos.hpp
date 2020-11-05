#ifndef ALGOS_HPP
#define ALGOS_HPP
/**
 * Code Author: Danny Sharp
 * This file is part of FFTE (Fast Fourier Transform Engine)
 */

#include <cstdint>
#include <cmath>
#include "complex.hpp"
#include "omega.hpp"
#include <iostream>

// Need forward declaration for the using directive
class constBiFuncNode;

/* This is the standard fft function call:
 * N, x, y, s_in, s_out.
 * N     = input signal length
 * x     = input signal
 * y     = output signal
 * s_in  = the distance (stride) between each input
 * s_out = the distance (stride) between each output
 */
using fft_fptr = void (*)(Complex*, Complex*, uint64_t, uint64_t, constBiFuncNode*, Omega&);

/* A class to map out what the call-graph of the FFT will look like, then
 * hold it in memory for the FFTs. Theoretically, it's compile-time ready.
 * However, due to the way that compilers work, this may or may not happen.
 */
class constBiFuncNode {
    public:
        fft_fptr fptr = nullptr; // FFT for this call
        uint64_t sz = 0;         // Size of FFT
        uint64_t left = 0;       // Offset in array until left child
        uint64_t right = 0;      // Offset in array until right child
        constexpr constBiFuncNode() = default; // Create default constructor
};

// We use this because it's a variable used in the factor function, so by defining it beforehand,
// going out of bounds on memory is much harder
#define FACTORS_LEN 15 

/* Statically allocated array of factors that are known at compile-time. These
 * are not necessarily prime, just ordered in the way that we prioritize.
 */
static constexpr uint64_t factors[FACTORS_LEN] {4, 2, 3, 5, 7, 11, 13, 16, 17, 19, 23, 29, 31, 37, 41};

// Function to find the smallest usable factor of f at compile-time.
constexpr uint64_t factor(const uint64_t f) {
    uint64_t k = 0;
    // Prioritize factors in the factors array
    for(; k < FACTORS_LEN; k++) {
        if( f % factors[k] == 0) return factors[k];
    }

    // If none of those work, turn to odd numbers that are greater than the last index
    for(k = factors[k - 1]; k*k < f; k+=2) {
        if( f % k == 0 ) return k;
    }

    // return f if no factor was found
    return f;
}

// Functions in the algos implementation file facing externally
Complex omega(uint power, uint N, Direction dir);
void pow2_FFT(Complex* x, Complex* y, uint64_t s_in, uint64_t s_out, constBiFuncNode* sRoot, Omega& w);

void DFT(Complex* x, Complex* y, uint64_t s_in, uint64_t s_out, constBiFuncNode* sLeaf, Omega& w);
void reference_DFT(uint64_t N, Complex* x, Complex* y, Direction dir);

void composite_FFT(Complex* x, Complex* y, uint64_t s_in, uint64_t s_out, constBiFuncNode* sRoot, Omega& w);
void reference_composite_FFT(uint64_t N, Complex* x, Complex* y, Direction dir);

void pow3_FFT(Complex* x, Complex* y, uint64_t s_in, uint64_t s_out, constBiFuncNode* sRoot, Omega& w);


#endif