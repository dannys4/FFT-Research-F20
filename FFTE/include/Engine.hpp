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

    // Interface for interacting with the backend of the engine
    template<typename F, int L>
    struct engine {
        // Most basic way of interaction with backend, normal arrays
        static inline void fft(Complex<F,L>* input, Complex<F,L>* output, int sig_length, Direction dir) {
            int num_nodes = getNumNodes(sig_length);
            biFuncNode<F, L> root[num_nodes];
            init_fft_tree(root, sig_length);
            root->fptr(input, output, 1, 1, root, dir);
        }

        // Interacting using vectors of FFTE complex numbers
        static inline complex_vector<F,L> fft(complex_vector<F,L> input, Direction dir) {
            auto N = input.size();
            auto out = input.get_allocator().allocate(N);
            int num_nodes = getNumNodes(N);
            biFuncNode<F,L> root[num_nodes];
            init_fft_tree(root, N);
            root->fptr(input.data(), out, 1, 1, root, dir);
            auto ret = complex_vector<F,L> (out, out+N);
            return ret;
        }

        // Interacting using vectors of std::complex numbers
        static inline std::vector<std::complex<F>*> fft(std::vector<std::complex<F>*> input, Direction dir) {
            auto N = input.size();
            auto in  = (Complex<F,L>*) aligned_alloc(alignof(Complex<F,L>), N*sizeof(Complex<F,L>));
            auto out = (Complex<F,L>*) aligned_alloc(alignof(Complex<F,L>), N*sizeof(Complex<F,L>));
            int num_nodes = getNumNodes(N);
            biFuncNode<F,L> root[num_nodes];
            init_fft_tree(root, N);
            for(int i = 0; i < N; i++) in[i] = Complex<F,L>(out[i]);
            root->fptr(in, out, 1, 1, root, dir);
            auto ret = std::vector<std::complex<F>*> {};
            for(int i = 0; i < N; i++) {
                std::complex<F>* tmp = new std::complex<F>[L];
                mm_store<F,L>::store(reinterpret_cast<F*>(tmp), out[i].get());
                ret.push_back(tmp);
            }
            return ret;
        }

    };
}

#endif // FFTE_ENGINE_HPP