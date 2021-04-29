/**
 * Code Author: Danny Sharp
 * This file is part of FFTE (Fast Fourier Transform Engine)
 */

#ifndef FFTE_ENGINE_HPP
#define FFTE_ENGINE_HPP

#define FFTE_IN_PARALLEL 0

#include "algos.hpp"
#include "allocator.hpp"
#if FFTE_MODERN_CPP
#include "constTree.hpp"
#else
#include "tree.hpp"
#endif

namespace FFTE {
    template<typename F, int L>
    using std_arrvec = std::array<std::vector<std::complex<F>>,L>;

    enum Major {column, row};

    // Shortcut to allocate a complex array. This allocates everything aligned in the
    // correct manner.
    template<typename F, int L>
    struct complex_alloc {
        static inline Complex<F,L>* alloc(size_t N) {
            return (Complex<F,L>*) aligned_alloc(alignof(Complex<F,L>), N*sizeof(Complex<F,L>));
        }
    };

    // Interface for interacting with the backend of the engine
    template<typename F, int L>
    struct engine {
        // Most basic way of interaction with backend, normal arrays
        static inline void fft(Complex<F,L>* input, Complex<F,L>* output, size_t sig_length, Direction dir, int stride_in, int stride_out) {
            int num_nodes = getNumNodes(sig_length);
            biFuncNode<F, L> root[num_nodes];
            init_fft_tree(root, sig_length);
            root->fptr(input, output, stride_in, stride_out, root, dir);
        }

        static inline void fft(Complex<F,L>* input, Complex<F,L>* output, size_t sig_length, Direction dir) {
            engine::fft(input, output, sig_length, dir, 1, 1);
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
        static inline std::vector<std::complex<F>*> fft(std::vector<std::complex<F>*>& input, Direction dir) {
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
                std::complex<F>* tmp = new std::complex<F>[L/2];
                mm_store<F,L>::store(reinterpret_cast<F*>(tmp), out[i].get());
                ret.push_back(tmp);
            }
            return ret;
        }

        // Perform an FFT on each of L/2 vectors of std::complex<F>
        static inline std::vector<std::complex<F>>* fft(std::vector<std::complex<F>>* input, Direction dir) {
            auto N = input[0].size();
            auto in = (Complex<F,L>*) aligned_alloc(alignof(Complex<F,L>), N*sizeof(Complex<F,L>));
            auto out = (Complex<F,L>*) aligned_alloc(alignof(Complex<F,L>), N*sizeof(Complex<F,L>));
            int num_nodes = getNumNodes(N);
            biFuncNode<F,L> root[num_nodes];
            init_fft_tree(root, N);
            auto tmp = new std::complex<F>[L/2];
            for(int i = 0; i < N; i++) {
                for(int j = 0; j < L/2; j++) tmp[j] = input[j][i];
                in[i] = Complex<F,L>(tmp);
            }
            root->fptr(in, out, 1, 1, root, dir);
            auto ret = new std::vector<std::complex<F>>[L/2];
            
            for(int i = 0; i < N; i++) {
                mm_store<F,L>::store(reinterpret_cast<F*>(tmp), out[i].get());
                for(int j = 0; j < L/2; j++) ret[j].push_back(tmp[j]);
            }
            delete[] tmp;
            return ret;
        }

        // Perform fft on one vector of std::complex numbers in Complex<F,L> types (with undefined last L/2-1 numbers when L > 2)
        static inline std::vector<std::complex<F>> fft(std::vector<std::complex<F>>& input, Direction dir) {
            auto N = input.size();
            auto in = (Complex<F,L>*) aligned_alloc(alignof(Complex<F,L>), N*sizeof(Complex<F,L>));
            auto out = (Complex<F,L>*) aligned_alloc(alignof(Complex<F,L>), N*sizeof(Complex<F,L>));
            int num_nodes = getNumNodes(N);
            biFuncNode<F,L> root[num_nodes];
            init_fft_tree(root, N);
            for(int i = 0; i < N; i++) {
                in[i] = Complex<F,L>(input[i]);
            }
            root->fptr(in, out, 1, 1, root, dir);
            auto ret = std::vector<std::complex<F>>{};
            
            std::complex<F> tmp[L/2];
            for(int i = 0; i < N; i++) {
                mm_store<F,L>::store(reinterpret_cast<F*>(tmp), out[i].get());
                ret.push_back(tmp[0]);
            }
            return ret;
        }

        // Perform fft on one input and put in predetmined location.
        static inline void fft(std::vector<std::complex<F>>& input, std::vector<std::complex<F>>& output, Direction dir) {
            size_t N = input.size();
            auto in = (Complex<F,L>*) aligned_alloc(alignof(Complex<F,L>), N*sizeof(Complex<F,L>));
            auto out = (Complex<F,L>*) aligned_alloc(alignof(Complex<F,L>), N*sizeof(Complex<F,L>));
            int num_nodes = getNumNodes(N);
            biFuncNode<F,L> root[num_nodes];
            init_fft_tree(root, N);
            for(size_t i = 0; i < N; i++) {
                in[i] = Complex<F,L>(input[i]);
            }
            root->fptr(in, out, 1, 1, root, dir);
            
            std::complex<F> tmp[L/2];
            for(size_t i = 0; i < N; i++) {
                mm_store<F,L>::store(reinterpret_cast<F*>(tmp), out[i].get());
                output.push_back(tmp[0]);
            }
        }

        // Perform fft on a set of L/2 inputs and put it in a set of L/2 outputs
        static inline void fft(std::vector<std::complex<F>>* input, std::vector<std::complex<F>>* output, Direction dir) {
            auto N = input[0].size();
            auto in = (Complex<F,L>*) aligned_alloc(alignof(Complex<F,L>), N*sizeof(Complex<F,L>));
            auto out = (Complex<F,L>*) aligned_alloc(alignof(Complex<F,L>), N*sizeof(Complex<F,L>));
            int num_nodes = getNumNodes(N);
            biFuncNode<F,L> root[num_nodes];
            init_fft_tree(root, N);
            auto tmp = new std::complex<F>[L/2];
            for(size_t i = 0; i < N; i++) {
                for(int j = 0; j < L/2; j++) tmp[j] = input[j][i];
                in[i] = Complex<F,L>(tmp);
            }
            root->fptr(in, out, 1, 1, root, dir);

            for(size_t i = 0; i < N; i++) {
                mm_store<F,L>::store(reinterpret_cast<F*>(tmp), out[i].get());
                for(int j = 0; j < L/2; j++) output[j].push_back(tmp[j]);
            }
            delete[] tmp;
        }
    };

    // Perform an arbitrary number of 1-D ffts on floats
    template<int P>
    std_arrvec<float, P> batch_fft(std_arrvec<float, P>& input, Direction dir) {
        std_arrvec<float, P> ret {};
        auto in_ptr = input.data();
        auto out_ptr = ret.data();
        #if FFTE_IN_PARALLEL
        #pragma omp parallel for
        #endif
        for(int p = 0; p <= P-4; p += 4) {
            engine<float, 8>::fft(in_ptr+p, out_ptr+p, dir);
        }
        if(P % 4 >= 2) {
            engine<float,4>::fft(in_ptr+P-(P%4), out_ptr+P-(P%4), dir);
        }
        if(P % 2 >= 1) {
            engine<float, 4>::fft(*(in_ptr+P-1), *(out_ptr+P-1), dir);
        }
        return ret;
    }

    template<size_t P, size_t Q>
    struct dim2engine {
        static inline std::array<std::array<std::complex<double>,Q>,P> fft2(std::array<std::array<std::complex<double>,Q>,P>& input, Direction dir, Major maj) {
            auto out1_ptr = (Complex<double, 2>*) aligned_alloc(alignof(Complex<double,2>), sizeof(Complex<double,2>)*P*Q);
            auto out2_ptr = (Complex<double, 2>*) aligned_alloc(alignof(Complex<double,2>), sizeof(Complex<double,2>)*P*Q);

            for(size_t p = 0; p < P; p++) {
                for(size_t q = 0; q < Q; q++) {
                    out1_ptr[Q*p + q] = Complex<double, 2>(&input[p][q]);
                }
            }
            
            switch(maj) {
                case row: 
                    #if FFTE_IN_PARALLEL
                    #pragma omp parallel for num_threads(4)
                    #endif
                    for(size_t p = 0; p < P; p++) {
                        engine<double,2>::fft(&out1_ptr[p*Q], &out2_ptr[p], Q, dir, 1, Q);
                    }

                    #if FFTE_IN_PARALLEL
                    #pragma omp parallel for num_threads(4)
                    #endif
                    for(size_t p = 0; p < P; p++) {
                        engine<double,2>::fft(&out2_ptr[p*Q], &out1_ptr[p], Q, dir, 1, Q);
                    }
                    break;
                case column:
                    #if FFTE_IN_PARALLEL
                    #pragma omp parallel for num_threads(4)
                    #endif
                    for(size_t p = 0; p < P; p++) {
                        engine<double,2>::fft(&out1_ptr[p], &out2_ptr[p*Q], Q, dir, Q, 1);
                    }
                    break;

                    #if FFTE_IN_PARALLEL
                    #pragma omp parallel for num_threads(4)
                    #endif
                    for(size_t p = 0; p < P; p++) {
                        engine<double,2>::fft(&out2_ptr[p], &out1_ptr[p*Q], Q, dir, Q, 1);
                    }
            }

            free(out2_ptr);
            auto output = std::array<std::array<std::complex<double>,Q>,P>();
            for(size_t p = 0; p < P; p++) {
                for(size_t q = 0; q < Q; q++) {
                    auto tmp = out1_ptr[Q*p + q].get();
                    output[p][q] = std::complex<double>(tmp[0], tmp[1]);
                }
            }
            free(out1_ptr);
            return output;
        }

        static inline Complex<double,2>* fft2(Complex<double,2>* input, Direction dir, Major maj) {
            auto out1_ptr = (Complex<double, 2>*) aligned_alloc(alignof(Complex<double,2>), sizeof(Complex<double,2>)*P*Q);
            auto out2_ptr = (Complex<double, 2>*) aligned_alloc(alignof(Complex<double,2>), sizeof(Complex<double,2>)*P*Q);
            
            
            
            switch(maj) {
                case row: 
                    #if FFTE_IN_PARALLEL
                    #pragma omp parallel for num_threads(4)
                    #endif
                    for(size_t p = 0; p < P; p++) {
                        engine<double,2>::fft(&input[p*Q], &out2_ptr[p], Q, dir, 1, Q);
                    }

                    #if FFTE_IN_PARALLEL
                    #pragma omp parallel for num_threads(4)
                    #endif
                    for(size_t p = 0; p < P; p++) {
                        engine<double,2>::fft(&out2_ptr[p*Q], &out1_ptr[p], Q, dir, 1, Q);
                    }
                    break;
                case column:
                    #if FFTE_IN_PARALLEL
                    #pragma omp parallel for num_threads(4)
                    #endif
                    for(size_t p = 0; p < P; p++) {
                        engine<double,2>::fft(&input[p], &out2_ptr[p*Q], Q, dir, Q, 1);
                    }
                    break;

                    #if FFTE_IN_PARALLEL
                    #pragma omp parallel for num_threads(4)
                    #endif
                    for(size_t p = 0; p < P; p++) {
                        engine<double,2>::fft(&out2_ptr[p], &out1_ptr[p*Q], Q, dir, Q, 1);
                    }
            }

            free(out2_ptr);
            return out1_ptr;
        }
    };

    // Perform an arbitrary number of 1-D ffts on doubles
    template<int P>
    std_arrvec<double, P> batch_fft(std_arrvec<double, P> input, Direction dir) {
        int p = 0;
        std_arrvec<double, P> ret {};
        auto sz = input[0].size();
        auto in_ptr = input.data();
        auto out_ptr = ret.data();
        #if FFTE_IN_PARALLEL
        #pragma omp parallel for num_threads(4)
        #endif
        for(int p = 0; p <= P-2; p += 2) {
            engine<double,4>::fft(in_ptr+p, out_ptr+p, sz, dir);
        }
        if( P % 2 >= 1) {
            engine<double,2>::fft(in_ptr+P-1, out_ptr+P-1, sz, dir);
        }
        return ret;
    }
}

#endif // FFTE_ENGINE_HPP