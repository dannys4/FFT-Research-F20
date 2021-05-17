#ifndef FFTE_ALGOS_HPP
#define FFTE_ALGOS_HPP
/**
 * Code Author: Danny Sharp
 * This file is part of FFTE (Fast Fourier Transform Engine)
 */
#include "complex.hpp"
#include "direction.hpp"
#include <cmath>
#include <vector>
#include <iostream>

// Check if we can use modern C++ features
#if (defined(__cplusplus) && ((__cplusplus / 100) >= 2014))
#define FFTE_MODERN_CPP 1
#else
#define FFTE_MODERN_CPP 0
#endif

namespace FFTE {
    // Need forward declaration for the using directive
    template<typename F, int L>
    class biFuncNode;

    // Functions in the algos implementation file facing externally
    
    template<typename F, int L>
    struct omega {
        static inline Complex<F, L> get(size_t power, size_t N, Direction dir) {
            F a = 2.*M_PI*((F) power)/((F) N);
            return Complex<F,L>(cos(a), static_cast<int>(dir)*sin(a));
        }
    };

    template<typename F, int L>
    inline void pow2_FFT(Complex<F,L>* x, Complex<F,L>* y, size_t s_in, size_t s_out, biFuncNode<F,L>* sRoot, Direction dir);

    template<typename F, int L>
    inline void DFT(Complex<F,L>* x, Complex<F,L>* y, size_t s_in, size_t s_out, biFuncNode<F,L>* sLeaf, Direction dir);
    
    template<typename F, int L>
    inline void reference_DFT(size_t N, Complex<F,L>* x, Complex<F,L>* y, Direction dir);


    template<typename F, int L>
    inline void composite_FFT(Complex<F,L>* x, Complex<F,L>* y, size_t s_in, size_t s_out, biFuncNode<F,L>* sRoot, Direction dir);

    template<typename F, int L>
    inline void reference_composite_FFT(size_t N, Complex<F,L>* x, Complex<F,L>* y, Direction dir);


    template<typename F, int L>
    inline void pow3_FFT(Complex<F,L>* x, Complex<F,L>* y, size_t s_in, size_t s_out, biFuncNode<F,L>* sRoot, Direction dir);

    template<typename F, int L>
    inline void rader_FFT(Complex<F,L>* x, Complex<F,L>* y, size_t s_in, size_t s_out, biFuncNode<F,L>* sRoot, Direction dir, size_t a, size_t ainv);

    enum fft_type {pow2, pow3, pow4, composite, discrete, rader};

#if FFTE_MODERN_CPP
    // Functor class for performing these transforms
    template<typename F, int L>
    class Fourier_Transform {
        private:
            fft_type type;
        public:
            constexpr Fourier_Transform(fft_type fft): type(fft) {}
            constexpr Fourier_Transform(Fourier_Transform& fft): type(fft.type) {}
            void operator()(Complex<F,L>* x, Complex<F,L>* y, size_t s_in, size_t s_out, biFuncNode<F,L>* sRoot, Direction dir) {
                switch(type) {
                    case fft_type::pow2: pow2_FFT(x, y, s_in, s_out, sRoot, dir); break;
                    case fft_type::pow3: pow3_FFT(x, y, s_in, s_out, sRoot, dir); break;
                    case fft_type::pow4: pow2_FFT(x, y, s_in, s_out, sRoot, dir); break;
                    case fft_type::composite: composite_FFT(x, y, s_in, s_out, sRoot, dir); break;
                    case fft_type::discrete: DFT(x, y, s_in, s_out, sRoot, dir); break;
                    default: std::cerr << "This is not supported or implemented yet\n"; exit(-1);
                }
            }
    };

    /* A class to map out what the call-graph of the FFT will look like, then
    * hold it in memory for the FFTs. Theoretically, it's compile-time ready.
    * However, due to the way that compilers work, this may or may not happen.
    */
    template<typename F, int L>
    class biFuncNode {
        public:
            Fourier_Transform<F,L> fptr; // FFT for this call
            size_t sz = 0;         // Size of FFT
            size_t left = 0;       // Offset in array until left child
            size_t right = 0;      // Offset in array until right child
            constexpr biFuncNode(): fptr(fft_type::discrete) {};
            constexpr biFuncNode(fft_type type): fptr(type) {}; // Create default constructor
            constexpr biFuncNode operator=(const biFuncNode& o) {fptr = o.fptr; return *this;}
    };

#else
    // Functor class for performing these transforms
    template<typename F, int L>
    class Fourier_Transform {
        protected:
            fft_type type;
            size_t root,root_inv;
        public:
            explicit Fourier_Transform() = default;
            explicit Fourier_Transform(fft_type fft): type(fft) {}
            explicit Fourier_Transform(size_t a, size_t ainv) {
                root = a; root_inv = ainv; type = fft_type::rader;
            }
            virtual void operator()(Complex<F,L>* x, Complex<F,L>* y, size_t s_in, size_t s_out, biFuncNode<F,L>* sRoot, Direction dir) {
                switch(type) {
                    case fft_type::pow2: pow2_FFT(x, y, s_in, s_out, sRoot, dir); break;
                    case fft_type::pow3: pow3_FFT(x, y, s_in, s_out, sRoot, dir); break;
                    case fft_type::pow4: pow2_FFT(x, y, s_in, s_out, sRoot, dir); break;
                    case fft_type::composite: composite_FFT(x, y, s_in, s_out, sRoot, dir); break;
                    case fft_type::discrete: DFT(x, y, s_in, s_out, sRoot, dir); break;
                    case fft_type::rader: rader_FFT(x, y, s_in, s_out, sRoot, dir, root, root_inv); break;
                    default: std::cerr << "Problem creating Fourier Transform functor\n"; exit(-1);
                }
            }
    };

    /* A class to map out what the call-graph of the FFT will look like, then
    * hold it in memory for the FFTs. Theoretically, it's compile-time ready.
    * However, due to the way that compilers work, this may or may not happen.
    */
   template<typename F, int L>
    class biFuncNode {
        public:
            Fourier_Transform<F,L> fptr; // FFT for this call
            size_t sz = 0;               // Size of FFT
            size_t left = 0;             // Offset in array until left child
            size_t right = 0;            // Offset in array until right child
            biFuncNode(): fptr(fft_type::discrete) {};
            biFuncNode(fft_type type): fptr(type) {}; // Create default constructor
            biFuncNode(size_t a, size_t ainv): fptr(a,ainv) {};
            biFuncNode<F,L> operator=(const biFuncNode<F,L> o) {fptr = o.fptr; return *this;}
    };
#endif // END DECISIONS BASED ON WHETHER C++14+
}

#include "algos_imp.hpp"
#endif // END FFTE_ALGOS_HPP