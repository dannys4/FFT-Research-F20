/**
 * Code Author: Danny Sharp
 * This file is part of the implementation for a stock FFT algorithm intended for HeFFTe
 */

#ifndef STOCK_FFT_PLAN_H
#define STOCK_FFT_PLAN_H

#include <exception>

#include "tree.h"

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
    template<>
    struct stock_fft_plan_dft<float> {
        size_t N, P, stride, idist, odist; Direction dir;
        biFuncNode<float, 8>* root;
        stock_fft_plan_dft(size_t N, int P, int stride,
                        int idist, int odist, Direction dir):
                        N(N), P(P), stride(stride), idist(idist), odist(odist), dir(dir) {
            int numNodes = getNumNodes(N);
            root = new biFuncNode<float, 8>[numNodes];
            init_fft_tree(root, N);
        }
        stock_fft_plan_dft(const stock_fft_plan_dft& plan): N(plan.N), P(plan.P), stride(plan.stride), idist(plan.idist), odist(plan.odist), dir(plan.dir) {}

        void execute(std::complex<float> data[]) {
            // TODO: Incorporate OMP

            // Allocate input and output temporary arrays
            auto inp = complex_alloc<float,8>::alloc(N);
            auto out = complex_alloc<float,8>::alloc(N);

            // Perform batch transform on everything save for the remainder
            for(int p = 0; p < P-3; p += 4) {
                // Convert types
                for(int i = 0; i < N; i++) {
                    auto idx = p*idist + i*stride;
                    inp[i] = Complex<float, 8> {data[idx+0*idist], data[idx+1*idist],
                                                data[idx+2*idist], data[idx+3*idist]};
                }
                // Perform fft
                root->fptr(inp, out, 1, 1, root, dir);
                // Convert type back
                for(int i = 0; i < N; i++) {
                    auto idx = p*odist + i*stride;
                    auto arr = reinterpret_cast<std::complex<float>*>(&out[i]);
                    data[idx+0*odist] = arr[0]; data[idx+1*odist] = arr[1];
                    data[idx+2*odist] = arr[2]; data[idx+3*odist] = arr[3];
                }
            }

            // Handle remainder
            if(P%4 > 0) {
                auto K = P % 4;
                // Init p for ease of use
                auto p = P - K;
                for(int i = 0; i < N; i++) {
                    auto idx = p*idist + i*stride;
                    // remainder columns are all zeros
                    auto z = std::complex<float> {};
                    switch(K) {
                        case 1: inp[i] = Complex<float,8> {data[idx+0*idist],
                                                           z,
                                                           z,
                                                           z};
                                break;
                        case 2: inp[i] = Complex<float,8> {data[idx+0*idist],
                                                           data[idx+1*idist],
                                                           z,
                                                           z};
                                break;
                        case 3: inp[i] = Complex<float,8> {data[idx+0*idist],
                                                           data[idx+1*idist],
                                                           data[idx+2*idist],
                                                           z};
                                break;
                        default: throw std::runtime_error("Something went wrong in stock fft!\n");
                    }
                }
                root->fptr(inp, out, 1, 1, root, dir);
                for(int i = 0; i < N; i++) {
                    auto idx = p*odist + i*stride;
                    auto arr = reinterpret_cast<std::complex<double>*>(&out[i]);
                    // Only need first k columns
                    for(int k = 0; k < K; k++) {
                        data[idx + k*odist] = arr[k];
                    }
                }
            }
        }
        void execute(std::complex<float> const idata[], float odata[]) {
            // TODO: Incorporate OMP

            // Allocate input and output temporary arrays
            auto inp = complex_alloc<float,8>::alloc(N);
            auto out = complex_alloc<float,8>::alloc(N);

            // Perform batch transform on everything save for the remainder
            for(int p = 0; p < P-3; p += 4) {
                // Convert types
                for(int i = 0; i < N; i++) {
                    auto idx = p*idist + i*stride;
                    inp[i] = Complex<float, 8> {idata[idx+0*idist], idata[idx+1*idist],
                                                idata[idx+2*idist], idata[idx+3*idist]};
                }
                // Perform fft
                root->fptr(inp, out, 1, 1, root, dir);
                // Convert type back
                for(int i = 0; i < N; i++) {
                    auto idx = p*odist + i*stride;
                    auto arr = reinterpret_cast<std::complex<float>*>(&out[i]);
                    odata[idx+0*odist] = arr[0].real(); odata[idx+1*odist] = arr[1].real();
                    odata[idx+2*odist] = arr[2].real(); odata[idx+3*odist] = arr[3].real();
                }
            }

            // Handle remainder
            if(P%4 > 0) {
                auto K = P % 4;
                // Init p for ease of use
                auto p = P - K;
                for(int i = 0; i < N; i++) {
                    auto idx = p*idist + i*stride;
                    // remainder columns are all zeros
                    auto z = std::complex<float> {};
                    switch(K) {
                        case 1: inp[i] = Complex<float,8> {idata[idx+0*idist],
                                                           z,
                                                           z,
                                                           z};
                                break;
                        case 2: inp[i] = Complex<float,8> {idata[idx+0*idist],
                                                           idata[idx+1*idist],
                                                           z,
                                                           z};
                                break;
                        case 3: inp[i] = Complex<float,8> {idata[idx+0*idist],
                                                           idata[idx+1*idist],
                                                           idata[idx+2*idist],
                                                           z};
                                break;
                        default: throw std::runtime_error("Something went wrong in stock fft!\n");
                    }
                }
                root->fptr(inp, out, 1, 1, root, dir);
                for(int i = 0; i < N; i++) {
                    auto idx = p*odist + i*stride;
                    auto arr = reinterpret_cast<std::complex<double>*>(&out[i]);
                    // Only need first k columns
                    for(int k = 0; k < K; k++) {
                        odata[idx + k*odist] = arr[k].real();
                    }
                }
            }
        }
        void execute(float const idata[], std::complex<float> odata[]) {
            // TODO: Incorporate OMP

            // Allocate input and output temporary arrays
            auto inp = complex_alloc<float,8>::alloc(N);
            auto out = complex_alloc<float,8>::alloc(N);

            // Perform batch transform on everything save for the remainder
            for(int p = 0; p < P-3; p += 4) {
                // Convert types
                for(int i = 0; i < N; i++) {
                    auto idx = p*idist + i*stride;
                    inp[i] = Complex<float, 8> {idata[idx+0*idist], 0, idata[idx+1*idist], 0,
                                                idata[idx+2*idist], 0, idata[idx+3*idist], 0};
                }
                // Perform fft
                root->fptr(inp, out, 1, 1, root, dir);
                // Convert type back
                for(int i = 0; i < N; i++) {
                    auto idx = p*odist + i*stride;
                    auto arr = reinterpret_cast<std::complex<float>*>(&out[i]);
                    odata[idx+0*odist] = arr[0]; odata[idx+1*odist] = arr[1];
                    odata[idx+2*odist] = arr[2]; odata[idx+3*odist] = arr[3];
                }
            }

            // Handle remainder
            if(P%4 > 0) {
                auto K = P % 4;
                // Init p for ease of use
                auto p = P - K;
                for(int i = 0; i < N; i++) {
                    auto idx = p*idist + i*stride;
                    // remainder columns are all zeros
                    switch(K) {
                        case 1: inp[i] = Complex<float,8> {idata[idx+0*idist], 0,
                                                           0                 , 0,
                                                           0                 , 0,
                                                           0                 , 0};
                                break;
                        case 2: inp[i] = Complex<float,8> {idata[idx+0*idist], 0,
                                                           idata[idx+1*idist], 0,
                                                           0                 , 0,
                                                           0                 , 0};
                                break;
                        case 3: inp[i] = Complex<float,8> {idata[idx+0*idist], 0,
                                                           idata[idx+1*idist], 0,
                                                           idata[idx+2*idist], 0,
                                                           0                 , 0};
                                break;
                        default: throw std::runtime_error("Something went wrong in stock fft!\n");
                    }
                }
                root->fptr(inp, out, 1, 1, root, dir);
                for(int i = 0; i < N; i++) {
                    auto idx = p*odist + i*stride;
                    auto arr = reinterpret_cast<std::complex<double>*>(&out[i]);
                    // Only need first k columns
                    for(int k = 0; k < K; k++) {
                        odata[idx + k*odist] = arr[k];
                    }
                }
            }
        }
        // Destructor
        ~stock_fft_plan_dft() { delete[] root; }
    };

    // Create fft plan struct for double precision
    template<>
    struct stock_fft_plan_dft<double> {
        size_t N, P, stride, idist, odist; Direction dir;
        biFuncNode<double, 4>* root;

        // Constructor for plan, initializes root
        stock_fft_plan_dft(size_t N, int P, int stride,
                        int idist, int odist, Direction dir):
                        N(N), P(P), stride(stride), idist(idist), odist(odist), dir(dir) {
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
                    auto idx = p*idist + i*stride;
                    inp[i] = Complex<double, 4> {data[idx], data[idx+idist]};
                }
                // Perform fft
                root->fptr(inp, out, 1, 1, root, dir);
                // Convert type back
                for(int i = 0; i < N; i++) {
                    auto idx = p*odist + i*stride;
                    auto arr = reinterpret_cast<std::complex<double>*>(&out[i]);
                    data[idx] = arr[0]; data[idx+odist] = arr[1];
                }
            }

            // Handle remainder
            if(P%2 == 1) {
                // Init p for ease of use
                auto p = P-1;
                for(int i = 0; i < N; i++) {
                    auto idx = p*idist + i*stride;
                    // Second column is all zeros
                    inp[i] = Complex<double, 4> {data[idx], std::complex<double> {}};
                }
                root->fptr(inp, out, 1, 1, root, dir);
                for(int i = 0; i < N; i++) {
                    auto idx = p*odist + i*stride;
                    auto arr = reinterpret_cast<std::complex<double>*>(&out[i]);
                    // Only need first column
                    data[idx] = arr[0];
                }
            }
        }
        // Executes the plan on given data
        void execute(std::complex<double> const idata[], double odata[]) {
            // TODO: Incorporate OMP

            // Allocate input and output temporary arrays
            auto inp = complex_alloc<double,4>::alloc(N);
            auto out = complex_alloc<double,4>::alloc(N);

            // Perform batch transform on everything save for the remainder
            for(int p = 0; p < P-1; p += 2) {
                // Convert types
                for(int i = 0; i < N; i++) {
                    auto idx = p*idist + i*stride;
                    inp[i] = Complex<double, 4> {idata[idx], idata[idx+idist]};
                }
                // Perform fft
                root->fptr(inp, out, 1, 1, root, dir);
                // Convert type back
                for(int i = 0; i < N; i++) {
                    auto idx = p*odist + i*stride;
                    auto arr = reinterpret_cast<std::complex<double>*>(&out[i]);
                    odata[idx] = arr[0].real(); odata[idx+odist] = arr[1].real();
                }
            }

            // Handle remainder
            if(P%2 == 1) {
                // Init p for ease of use
                auto p = P-1;
                for(int i = 0; i < N; i++) {
                    auto idx = p*idist + i*stride;
                    // Second column is all zeros
                    inp[i] = Complex<double, 4> {idata[idx], std::complex<double> {}};
                }
                root->fptr(inp, out, 1, 1, root, dir);
                for(int i = 0; i < N; i++) {
                    auto idx = p*odist + i*stride;
                    auto arr = reinterpret_cast<std::complex<double>*>(&out[i]);
                    // Only need first column
                    odata[idx] = arr[0].real();
                }
            }
        }
        // Executes the plan on given data
        void execute(double const idata[], std::complex<double> odata[]) {
            // TODO: Incorporate OMP

            // Allocate input and output temporary arrays
            auto inp = complex_alloc<double,4>::alloc(N);
            auto out = complex_alloc<double,4>::alloc(N);

            // Perform batch transform on everything save for the remainder
            for(int p = 0; p < P-1; p += 2) {
                // Convert types
                for(int i = 0; i < N; i++) {
                    auto idx = p*idist + i*stride;
                    inp[i] = Complex<double, 4> {idata[idx], 0, idata[idx+idist], 0};
                }
                // Perform fft
                root->fptr(inp, out, 1, 1, root, dir);
                // Convert type back
                for(int i = 0; i < N; i++) {
                    auto idx = p*odist + i*stride;
                    auto arr = reinterpret_cast<std::complex<double>*>(&out[i]);
                    odata[idx] = arr[0]; odata[idx+odist] = arr[1];
                }
            }

            // Handle remainder
            if(P%2 == 1) {
                // Init p for ease of use
                auto p = P-1;
                for(int i = 0; i < N; i++) {
                    auto idx = p*idist + i*stride;
                    // Second column is all zeros
                    inp[i] = Complex<double, 4> {idata[idx], 0, 0, 0};
                }
                root->fptr(inp, out, 1, 1, root, dir);
                for(int i = 0; i < N; i++) {
                    auto idx = p*odist + i*stride;
                    auto arr = reinterpret_cast<std::complex<double>*>(&out[i]);
                    // Only need first column
                    odata[idx] = arr[0];
                }
            }
        }

        ~stock_fft_plan_dft() { delete[] root; }
    };

}

#endif // STOCK_FFT_PLAN_H