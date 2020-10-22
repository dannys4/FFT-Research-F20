#ifndef ALGOS_HPP
#define ALGOS_HPP

#include <cstdint>
#include <cmath>
#include "utils.hpp"
#include <iostream>

void pow2_FFT(Complex* sig_in, uint64_t stride, uint64_t size, Complex* sig_out);
void DFT(uint64_t size, Complex* sig_in, Complex* sig_out, uint64_t s_in, uint64_t s_out);
void composite_FFT(uint64_t N, Complex* x, Complex* y, uint64_t s_in, uint64_t s_out);
void pow3_FFT(uint64_t N, Complex* x, Complex* y, uint64_t s_in);
void pow2_FFT_it(Complex* sig_in, uint64_t size, Complex* sig_out);

#endif