/**
 * Code Author: Danny Sharp
 * This file is part of the implementation for a stock FFT algorithm intended for HeFFTe
 */

#ifndef STOCK_FFT_PLAN_H
#define STOCK_FFT_PLAN_H

#include "tree.hpp"

namespace stock_fft {

template<typename F> struct stock_fft_plan_dft {};

struct stock_fft_plan_dft<float> {
    int N, P, stride, dist; Direction dir;
    biFuncNode<float, 8>* tree;
    stock_fft_plan_dft(int N, int P, int stride,
                       int dist, Direction dir):
                       N(N), P(P), stride(stride), dist(dist), dir(dir) {}

};

}

#endif // STOCK_FFT_PLAN_H