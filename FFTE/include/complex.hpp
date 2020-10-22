#ifndef COMPLEX_HPP
#define COMPLEX_HPP

#include <immintrin.h>
#include <iostream>
#include <cmath>

class Complex {
    __m128d var;
    public:
        Complex(double a, double b) {
            var = _mm_set_pd(b, a);
        }

        Complex(__m128d& v) { var = v; }

        Complex operator+(const Complex& y) {
            __m128d ret = _mm_add_pd(var, y.var);
            return Complex(ret);
        }

        Complex operator-(const Complex& y) {
            __m128d ret = _mm_sub_pd(var, y.var);
            return Complex(ret);
        }

        Complex operator*(const Complex& y) {
            // auto cc = _mm_permute_pd(y.var, 0);
            // auto ba = _mm_permute_pd(this->var, 1);
            // auto dd = _mm_permute_pd(y.var, 3);
            // auto dba = _mm_mul_pd(ba, dd);
            // auto mult = _mm_fmaddsub_pd(this->var, cc, dba);
            auto mult = _mm_fmaddsub_pd(var, _mm_permute_pd(y.var, 0), _mm_mul_pd(_mm_permute_pd(var, 1), _mm_permute_pd(y.var, 3)));
            return Complex(mult);
        }

        Complex twiddle(int j, int L) {
            double frac = (double) j / (double) L;
            return Complex(cos(2*frac*M_PI), -sin(2*frac*M_PI));
        }

        __m128d getVar() const {
            return var;
        }

        // ATTENTION: DANGEROUS CAST!
        // It's really fast though, so ¯\_(ツ)_/¯
        double modulus() {
            auto dp = _mm_dp_pd(var, var, 0xF1);
            return *((double*) &dp);
        }
};

inline std::ostream& operator<<(std::ostream& os, const Complex& dt){
    double* tmp = (double*) aligned_alloc(16, 2*sizeof(double));
    _mm_store_pd(tmp, dt.getVar());
    if(tmp[1] > 0) os << tmp[0] << " + " << tmp[1] << "i";
    else os << tmp[0] << " - " << -tmp[1] << "i";
    free(tmp);
    return os;
}
#endif