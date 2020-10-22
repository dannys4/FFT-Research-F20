#ifndef TEST_HPP
#define TEST_HPP

#include <stdlib.h>
#include "algos.hpp"
#include <random>
#include <algorithm>
#include <functional>
#include <chrono>
#include <complex>
#include "constTree.hpp"

void check_fft();
volatile void check_fft_tree();

void time_fft();
void time_complex_mult();
void time_const_tree();

#endif