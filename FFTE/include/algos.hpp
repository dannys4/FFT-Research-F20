#ifndef ALGOS_HPP
#define ALGOS_HPP
/**
 * Code Author: Danny Sharp
 * This file is part of FFTE (Fast Fourier Transform Engine)
 */
#include <cmath>
#include "complex.hpp"
#include "omega.hpp"
#include <iostream>

#define MODERN_CPP (defined(__cplusplus) && ((__cplusplus / 100) >= 2014))

namespace FFTE {
    // Need forward declaration for the using directive
    class biFuncNode;

    /* This is the standard fft function call:
    * N, x, y, s_in, s_out.
    * N     = input signal length
    * x     = input signal
    * y     = output signal
    * s_in  = the distance (stride) between each input
    * s_out = the distance (stride) between each output
    */
    using fft_fptr = void (*)(Complex*, Complex*, size_t, size_t, biFuncNode*, Omega&);

    // Functions in the algos implementation file facing externally
    Complex omega(uint power, uint N, Direction dir);
    void pow2_FFT(Complex* x, Complex* y, size_t s_in, size_t s_out, biFuncNode* sRoot, Omega& w);

    void DFT(Complex* x, Complex* y, size_t s_in, size_t s_out, biFuncNode* sLeaf, Omega& w);
    void reference_DFT(size_t N, Complex* x, Complex* y, Direction dir);

    void composite_FFT(Complex* x, Complex* y, size_t s_in, size_t s_out, biFuncNode* sRoot, Omega& w);
    void reference_composite_FFT(size_t N, Complex* x, Complex* y, Direction dir);

    void pow3_FFT(Complex* x, Complex* y, size_t s_in, size_t s_out, biFuncNode* sRoot, Omega& w);

    enum fft_type {pow2, pow3, pow4, composite, discrete};

#if MODERN_CPP
    // Functor class for performing these transforms
    class Fourier_Transform {
        private:
            fft_type type;
        public:
            constexpr Fourier_Transform(fft_type fft): type(fft) {}
            constexpr Fourier_Transform(Fourier_Transform& fft): type(fft.type) {}
            void operator()(Complex* x, Complex* y, size_t s_in, size_t s_out, biFuncNode* sRoot, Omega& w) {
                switch(type) {
                    case fft_type::pow2: pow2_FFT(x, y, s_in, s_out, sRoot, w); break;
                    case fft_type::pow3: pow3_FFT(x, y, s_in, s_out, sRoot, w); break;
                    case fft_type::pow4: pow2_FFT(x, y, s_in, s_out, sRoot, w); break;
                    case fft_type::composite: composite_FFT(x, y, s_in, s_out, sRoot, w); break;
                    case fft_type::discrete: DFT(x, y, s_in, s_out, sRoot, w); break;
                    default: std::cerr << "NO\n"; exit(-1);
                }
            }
    };

    /* A class to map out what the call-graph of the FFT will look like, then
    * hold it in memory for the FFTs. Theoretically, it's compile-time ready.
    * However, due to the way that compilers work, this may or may not happen.
    */
    class biFuncNode {
        public:
            Fourier_Transform fptr; // FFT for this call
            size_t sz = 0;         // Size of FFT
            size_t left = 0;       // Offset in array until left child
            size_t right = 0;      // Offset in array until right child
            constexpr biFuncNode(): fptr(fft_type::discrete) {};
            constexpr biFuncNode(fft_type type): fptr(type) {}; // Create default constructor
            constexpr biFuncNode operator=(const biFuncNode& o) {fptr = o.fptr; return *this;}
    };

#else
    // Functor class for performing these transforms
    class Fourier_Transform {
        private:
            fft_type type;
        public:
            Fourier_Transform(fft_type fft): type(fft) {}
            void operator()(Complex* x, Complex* y, size_t s_in, size_t s_out, biFuncNode* sRoot, Omega& w) {
                switch(type) {
                    case fft_type::pow2: pow2_FFT(x, y, s_in, s_out, sRoot, w); break;
                    case fft_type::pow3: pow3_FFT(x, y, s_in, s_out, sRoot, w); break;
                    case fft_type::pow4: pow2_FFT(x, y, s_in, s_out, sRoot, w); break;
                    case fft_type::composite: composite_FFT(x, y, s_in, s_out, sRoot, w); break;
                    case fft_type::discrete: DFT(x, y, s_in, s_out, sRoot, w); break;
                    default: std::cerr << "NO\n"; exit(-1);
                }
            }
    };

    /* A class to map out what the call-graph of the FFT will look like, then
    * hold it in memory for the FFTs. Theoretically, it's compile-time ready.
    * However, due to the way that compilers work, this may or may not happen.
    */
    class biFuncNode {
        public:
            Fourier_Transform fptr; // FFT for this call
            size_t sz = 0;         // Size of FFT
            size_t left = 0;       // Offset in array until left child
            size_t right = 0;      // Offset in array until right child
            biFuncNode(): fptr(fft_type::discrete) {};
            biFuncNode(fft_type type): fptr(type) {}; // Create default constructor
            biFuncNode operator=(const biFuncNode o) {fptr = o.fptr; return *this;}
    };
#endif
}
#endif