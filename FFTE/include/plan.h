/**
 * Code Author: Danny Sharp
 * This file is part of the implementation for a stock FFT algorithm intended for HeFFTe
 */

#ifndef STOCK_FFT_PLAN_H
#define STOCK_FFT_PLAN_H

#include "tree.hpp"

namespace stock_fft {

    // Shortcut to allocate a complex array. This allocates everything aligned in the
    // correct manner.
    template<typename F, int L>
    struct complex_alloc {
        static inline Complex<F,L>* alloc(size_t N) {
            return (Complex<F,L>*) aligned_alloc(alignof(Complex<F,L>), N*sizeof(Complex<F,L>));
        }
    };

    template<typename F> struct stock_fft_plan_dft {};

    struct stock_fft_plan_dft<float> {
        size_t N, P, stride, dist; Direction dir;
        biFuncNode<float, 8>* root;
        stock_fft_plan_dft(size_t N, int P, int stride,
                        int dist, Direction dir):
                        N(N), P(P), stride(stride), dist(dist), dir(dir) {
            int numNodes = getNumNodes(N);
            root = new biFuncNode<float, 8>[numNodes];
            init_fft_tree(root, N);
        }

        ~stock_fft_plan_dft() { delete[] root; }
    };

    struct stock_fft_plan_dft<double> {
        size_t N, P, stride, dist; Direction dir;
        biFuncNode<double, 4>* root;
        stock_fft_plan_dft(size_t N, int P, int stride,
                        int dist, Direction dir):
                        N(N), P(P), stride(stride), dist(dist), dir(dir) {
            int numNodes = getNumNodes(N);
            root = new biFuncNode<double, 4>[numNodes];
            init_fft_tree(root, N);
        }

        void execute(std::complex<double> data[]) {
            // TODO: Incorporate OMP
            auto inp = complex_alloc<double,4>::alloc(N);
            auto out = complex_alloc<double,4>::alloc(N);
            for(int p = 0; p < P; p += 2) {
                
            }
        }

        ~stock_fft_plan_dft() { delete[] root; }
    };

}

#endif // STOCK_FFT_PLAN_H