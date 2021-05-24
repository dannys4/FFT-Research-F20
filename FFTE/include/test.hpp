#ifndef TEST_HPP
#define TEST_HPP

#include <stdlib.h>
#include "algos.hpp"
#include "Engine.hpp"
#include <random>
#include <algorithm>
#include <functional>
#include <chrono>
#include <complex>
#include <cassert>
#include <fstream>

#if STOCK_FFT_MODERN_CPP
#include "constTree.hpp"
#else
#include "tree.hpp"
#endif


// Check the error in the FFT in each direction
void check_fft(STOCK_FFT::Direction dir);

// Check the error in the FFT for multiple simultaneous packed transforms
void check_fft_multidim(STOCK_FFT::Direction dir);

// Check the way that the FFT branches
void check_fft_tree();

// Check the error in the omega values
void check_omega();

// Check all arithmetic operations in complex class
void check_complex_ops();

// Check examples of batch ffts
void check_batch_fft();

// Time an FFT
void time_fft();

// Time complex multiplication
std::pair<long long int, long long int> time_complex_mult(size_t, bool);

// Compare timing of const tree to runtime tree
// C++14+ only!
void time_const_tree();

// Time the batch FFT
void time_batch_fft();

// Time the FFTs for many numbers and input the times into a csv
void test_FFT_into_csv(std::string filename, int maxSize);
#endif