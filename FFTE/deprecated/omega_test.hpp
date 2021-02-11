/**
 * Code Author: Danny Sharp
 * This file is part of FFTE (Fast Fourier Transform Engine)
 */

#ifndef OMEGA_TEST_HPP
#define OMEGA_TEST_HPP

#include "omega.hpp"
#include "algos.hpp"
#include "complex.hpp"
#include <random>
#include <algorithm>
#include <functional>
#include <chrono>
#include <complex>
#include <cassert>
#include <fstream>


// Check the error in the omega values
void check_omega();

// Time using the omega class
void time_omega();

#endif // OMEGA_TEST_HPP