/**
 * Code Author: Danny Sharp
 * This file is part of FFTE (Fast Fourier Transform Engine)
 */

#ifndef FFTE_ENGINE_HPP
#define FFTE_ENGINE_HPP

#include "algos.hpp"
#include "allocator.hpp"

namespace FFTE {
    
    // Shortcut to allocate a complex array. This allocates everything aligned in the
    // correct manner.
    template<typename F, int L>
    Complex<F,L>* complex_alloc(int N) {
        return (Complex<F,L>*) aligned_alloc(alignof(Complex<F,L>), N*sizeof(Complex<F,L>));
    }

    // Interface for interacting with the backend of the whole ordeal.
    template<typename F, int L, int sig_length>
    void fft(Complex<F,L>* input, Complex<F,L>* output, Direction dir) {
        int num_nodes = getNumNodes(sig_length);
        biFuncNode<F, L> root[num_nodes];
        init_fft_tree(root, sig_length);
        root->fptr(input, output, 1, 1, root, dir);
    }

    template<typename F, int L>
    complex_vector<F,L> fft(complex_vector<F,L> input, Direction dir) {
        auto N = input.size();
        auto out = input.get_allocator().allocate(N);
        int num_nodes = getNumNodes(N);
        biFuncNode<F,L> root[num_nodes];
        init_fft_tree(root, N);
        root->fptr(input.data(), out, 1, 1, root, dir);
        auto ret = complex_vector<F,L> (out, out+N);
        return ret;
    }
}

#endif // FFTE_ENGINE_HPP