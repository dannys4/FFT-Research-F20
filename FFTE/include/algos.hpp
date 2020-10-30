#ifndef ALGOS_HPP
#define ALGOS_HPP

#include <cstdint>
#include <cmath>
#include "complex.hpp"
#include "omega.hpp"
#include <iostream>

class constBiFuncNode;

// This is the standard fft function call:
// N, x, y, s_in, s_out.
// N     = input signal length
// x     = input signal
// y     = output signal
// s_in  = the distance (stride) between each input
// s_out = the distance (stride) between each output
using fft_fptr = void (*)(Complex*, Complex*, uint64_t, uint64_t, constBiFuncNode*, Omega&);

template<typename T>
class constBiNode {
    private:
    public:
        T elem;
        uint64_t left = 0;
        uint64_t right = 0;
        constexpr constBiNode() = default;
        constexpr constBiNode(const T e) {elem = e;}
};

class constBiFuncNode {
    public:
        fft_fptr fptr = nullptr;
        uint64_t sz = 0;
        uint64_t left = 0;
        uint64_t right = 0;
        constexpr constBiFuncNode() = default;
};



#define FACTORS_LEN 15
static constexpr uint64_t factors[FACTORS_LEN] {4, 2, 3, 5, 7, 11, 13, 16, 17, 19, 23, 29, 31, 37, 41};

constexpr uint64_t factor(const uint64_t f) {
    uint64_t k = 0;
    for(; k < FACTORS_LEN; k++) {
        if( f % factors[k] == 0) return factors[k];
    }
    for(k = factors[k - 1]; k*k < f; k+=2) {
        if( f % k == 0 ) return k;
    }
    return f;
}

Complex omega(uint power, uint N);
void pow2_FFT(Complex* x, Complex* y, uint64_t s_in, uint64_t s_out, constBiFuncNode* sRoot, Omega& w);

void DFT(Complex* x, Complex* y, uint64_t s_in, uint64_t s_out, constBiFuncNode* sLeaf, Omega& w);
void reference_DFT(uint64_t N, Complex* x, Complex* y);

void composite_FFT(Complex* x, Complex* y, uint64_t s_in, uint64_t s_out, constBiFuncNode* sRoot, Omega& w);
void reference_composite_FFT(uint64_t N, Complex* x, Complex* y, uint64_t s_in, uint64_t s_out);

void pow3_FFT(Complex* x, Complex* y, uint64_t s_in, uint64_t s_out, constBiFuncNode* sRoot, Omega& w);


#endif