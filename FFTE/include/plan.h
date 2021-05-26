/**
 * Code Author: Danny Sharp
 * This file is part of the implementation for a stock FFT algorithm intended for HeFFTe
 */

#ifndef STOCK_FFT_PLAN_H
#define STOCK_FFT_PLAN_H

#include <exception>

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

    // Templated fft plan depending on precision
    template<typename F> struct stock_fft_plan_dft {};

    // Create fft plan for single precision
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

        void execute(std::complex<float> data[]) {
            // TODO: Incorporate OMP

            // Allocate input and output temporary arrays
            auto inp = complex_alloc<float,8>::alloc(N);
            auto out = complex_alloc<float,8>::alloc(N);

            // Perform batch transform on everything save for the remainder
            for(int p = 0; p < P-3; p += 4) {
                // Convert types
                for(int i = 0; i < N; i++) {
                    auto idx = p*dist + i*stride;
                    inp[i] = Complex<float, 8> {data[idx+0*dist], data[idx+1*dist],
                                                data[idx+2*dist], data[idx+3*dist]};
                }
                // Perform fft
                root->fptr(inp, out, 1, 1, root, dir);
                // Convert type back
                for(int i = 0; i < N; i++) {
                    auto idx = p*dist + i*stride;
                    auto arr = reinterpret_cast<std::complex<float>*>(&out[i]);
                    data[idx+0*dist] = arr[0]; data[idx+1*dist] = arr[1];
                    data[idx+2*dist] = arr[2]; data[idx+3*dist] = arr[3];
                }
            }

            // Handle remainder
            if(P%4 > 0) {
                auto K = P % 4;
                // Init p for ease of use
                auto p = P - K;
                for(int i = 0; i < N; i++) {
                    auto idx = p*dist + i*stride;
                    // remainder columns are all zeros
                    auto z = std::complex<float> {};
                    switch(K) {
                        case 1: inp[i] = Complex<float,8> {data[idx+0*dist],
                                                           z,
                                                           z,
                                                           z};
                                break;
                        case 2: inp[i] = Complex<float,8> {data[idx+0*dist],
                                                           data[idx+1*dist],
                                                           z,
                                                           z};
                                break;
                        case 3: inp[i] = Complex<float,8> {data[idx+0*dist],
                                                           data[idx+1*dist],
                                                           data[idx+2*dist],
                                                           z};
                                break;
                        default: throw std::runtime_error("Something went wrong in stock fft!\n");
                    }
                }
                root->fptr(inp, out, 1, 1, root, dir);
                for(int i = 0; i < N; i++) {
                    auto idx = p*dist + i*stride;
                    auto arr = reinterpret_cast<std::complex<double>*>(&out[i]);
                    // Only need first k columns
                    for(int k = 0; k < K; k++) {
                        data[idx + k*dist] = arr[k];
                    }
                }
            }
        }

        // Destructor
        ~stock_fft_plan_dft() { delete[] root; }
    };

    // Create fft plan struct for double precision
    struct stock_fft_plan_dft<double> {
        size_t N, P, stride, dist; Direction dir;
        biFuncNode<double, 4>* root;

        // Constructor for plan, initializes root
        stock_fft_plan_dft(size_t N, int P, int stride,
                        int dist, Direction dir):
                        N(N), P(P), stride(stride), dist(dist), dir(dir) {
            int numNodes = getNumNodes(N);
            root = new biFuncNode<double, 4>[numNodes];
            init_fft_tree(root, N);
        }

        // Executes the plan on given data
        void execute(std::complex<double> data[]) {
            // TODO: Incorporate OMP

            // Allocate input and output temporary arrays
            auto inp = complex_alloc<double,4>::alloc(N);
            auto out = complex_alloc<double,4>::alloc(N);

            // Perform batch transform on everything save for the remainder
            for(int p = 0; p < P-1; p += 2) {
                // Convert types
                for(int i = 0; i < N; i++) {
                    auto idx = p*dist + i*stride;
                    inp[i] = Complex<double, 4> {data[idx], data[idx+dist]};
                }
                // Perform fft
                root->fptr(inp, out, 1, 1, root, dir);
                // Convert type back
                for(int i = 0; i < N; i++) {
                    auto idx = p*dist + i*stride;
                    auto arr = reinterpret_cast<std::complex<double>*>(&out[i]);
                    data[idx] = arr[0]; data[idx+dist] = arr[1];
                }
            }

            // Handle remainder
            if(P%2 == 1) {
                // Init p for ease of use
                auto p = P-1;
                for(int i = 0; i < N; i++) {
                    auto idx = p*dist + i*stride;
                    // Second column is all zeros
                    inp[i] = Complex<double, 4> {data[idx], std::complex<double> {}};
                }
                root->fptr(inp, out, 1, 1, root, dir);
                for(int i = 0; i < N; i++) {
                    auto idx = p*dist + i*stride;
                    auto arr = reinterpret_cast<std::complex<double>*>(&out[i]);
                    // Only need first column
                    data[idx] = arr[0];
                }
            }
        }

        ~stock_fft_plan_dft() { delete[] root; }
    };

}

#endif // STOCK_FFT_PLAN_H