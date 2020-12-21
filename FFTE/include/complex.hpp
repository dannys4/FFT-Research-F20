#ifndef FFTE_COMPLEX_HPP
#define FFTE_COMPLEX_HPP
/**
 * Code Author: Danny Sharp
 * This file is part of FFTE (Fast Fourier Transform Engine)
 */

#include <immintrin.h>
#include <iostream>
#include <cmath>

namespace FFTE {
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
            Complex(__m128d& v) : var(v) {}

            // Creates a complex number from a pack
            Complex(__m128d v) : var(v) {}

            // Calculates the sum of two complex numbers
            Complex operator+(const Complex& y) {
                return Complex(_mm_add_pd(var, y.var));
            }

            // Calculates the difference of two complex numbers
            Complex operator-(const Complex& y) {
                return Complex(_mm_sub_pd(var, y.var));
            }

            // Divides a complex number by a double
            Complex operator/(const double y) {
                return Complex(_mm_div_pd(var, _mm_set_pd1(y)));
            }

            // Divides a complex number by a double and resets this value
            Complex operator/=(const double y) {
                var = _mm_div_pd(var, _mm_set_pd1(y));
                return *this;
            }

            // Return the complex conjugate of this double
            Complex conjugate() {
                return Complex(_mm_addsub_pd(_mm_setzero_pd(), var));
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
                return Complex(_mm_fmaddsub_pd(var, _mm_permute_pd(y.var, 0), _mm_mul_pd(_mm_permute_pd(var, 1), _mm_permute_pd(y.var, 3))));
            }

            // Performs multiplication between this and y then resets this
            Complex operator*=(const Complex& y) {
                var = _mm_fmaddsub_pd(var, _mm_permute_pd(y.var, 0), _mm_mul_pd(_mm_permute_pd(var, 1), _mm_permute_pd(y.var, 3)));
                return *this;
            }

            // Performs multiplication between a double and Complex number
            Complex operator*(double y) {
                return Complex(_mm_mul_pd(_mm_set_pd1(y), var));
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
        double tmp[2];
        _mm_storeu_pd(tmp, dt.getVar());
        if(tmp[1] > 0) os << tmp[0] << " + " << tmp[1] << "i";
        else os << tmp[0] << " - " << -tmp[1] << "i";
        return os;
    }

    // Define multiplication between double and Complex number
    inline Complex operator*(double d, Complex c) {
        return c*d;
    }

    // Define division of a double by a complex number
    inline Complex operator/(double d, Complex c) {
        Complex ret = d * c.conjugate() / c.modulus();
        return ret;
    }
}
#endif // END FFTE_COMPLEX_HPP