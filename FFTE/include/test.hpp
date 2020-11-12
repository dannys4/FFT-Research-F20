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

#define MODERN_CPP (defined(__cplusplus) && ((__cplusplus / 100) >= 2014))
#if MODERN_CPP
#include "constTree.hpp"
#else
#include "tree.hpp"
#endif



void check_fft(Direction dir);
void check_fft_tree();
void check_omega();

void time_fft();
void time_complex_mult();
void time_const_tree();
void time_omega();

#endif