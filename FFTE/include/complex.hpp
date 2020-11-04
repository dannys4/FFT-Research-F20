#ifndef COMPLEX_HPP
#define COMPLEX_HPP
/**
 * Code Author: Danny Sharp
 * This file is part of FFTE (Fast Fourier Transform Engine)
 */

#include <immintrin.h>
#include <iostream>
#include <cmath>

class Complex {
    
    // A pack of two doubles that will hold the complex number
    __m128d var;
    
    public:

        // Default complex number constructor
        Complex() {
            var = _mm_setzero_pd();
        }

        // Creates a complex number from two doubles
        Complex(double a, double b) {
            var = _mm_set_pd(b, a);
        }

        // Creates a complex number from a pack
        Complex(__m128d& v) { var = v; }

        // Calculates the sum of two complex numbers
        Complex operator+(const Complex& y) {
            __m128d ret = _mm_add_pd(var, y.var);
            return Complex(ret);
        }

        // Calculates the difference of two complex numbers
        Complex operator-(const Complex& y) {
            __m128d ret = _mm_sub_pd(var, y.var);
            return Complex(ret);
        }

        /* 
         * Performs multiplication for two complex numbers. Equivalent code below:
         * auto cc = _mm_permute_pd(y.var, 0);
         * auto ba = _mm_permute_pd(this->var, 1);
         * auto dd = _mm_permute_pd(y.var, 3);
         * auto dba = _mm_mul_pd(ba, dd);
         * auto mult = _mm_fmaddsub_pd(this->var, cc, dba);
         */
        Complex operator*(const Complex& y) {
            auto mult = _mm_fmaddsub_pd(var, _mm_permute_pd(y.var, 0), _mm_mul_pd(_mm_permute_pd(var, 1), _mm_permute_pd(y.var, 3)));
            return Complex(mult);
        }

        // Returns the var member of this class
        __m128d getVar() const {
            return var;
        }

        /* 
         * Gets the modulus of the complex number
         * NOTE: a dangerous cast is performed here which depends on how AVX2 is
         * formatted in memory. Could be different per machine, but I doubt it.
         */
        double modulus() {
            auto dp = _mm_dp_pd(var, var, 0xF1);
            return *((double*) &dp);
        }
};

// Determines how my complex number should be printed to an ostream
inline std::ostream& operator<<(std::ostream& os, const Complex& dt){
    double* tmp = (double*) aligned_alloc(16, 2*sizeof(double));
    _mm_store_pd(tmp, dt.getVar());
    if(tmp[1] > 0) os << tmp[0] << " + " << tmp[1] << "i";
    else os << tmp[0] << " - " << -tmp[1] << "i";
    free(tmp);
    return os;
}
#endif