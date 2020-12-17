#ifndef TEST_HPP
#define TEST_HPP

#include <stdlib.h>
#include "algos.hpp"
#include <random>
#include <algorithm>
#include <functional>
#include <chrono>
#include <complex>
#include "omega.hpp"
#include <cassert>
#include <fstream>

#if MODERN_CPP
#include "constTree.hpp"
#else
#include "tree.hpp"
#endif


// Check the error in the FFT in each direction
void check_fft(FFTE::Direction dir);

// Check the way that the FFT branches
void check_fft_tree();

// Check the error in the omega values
void check_omega();

// Time an FFT
void time_fft();

// Time complex multiplication
void time_complex_mult();

// Compare timing of const tree to runtime tree
// C++14+ only!
void time_const_tree();

// Time using the omega class
void time_omega();

// Time the FFTs for many numbers and input the times into a csv
void test_FFT_into_csv(std::string filename, int maxSize);
#endif