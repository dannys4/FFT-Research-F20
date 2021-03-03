/**
 * Code Author: Danny Sharp
 * This file is part of FFTE (Fast Fourier Transform Engine)
 */

#ifndef FFTE_COMPLEX_HPP
#define FFTE_COMPLEX_HPP
#include <type_traits>
#include <iostream>
#include "vec_types.hpp"

/**
 * This file holds the tools of using Complex numbers in a templated fashion
 * using vectorized intrinsics.
 */
namespace FFTE {
    template<typename F, int L>
    class alignas(L*sizeof(F)) Complex {
        public:
            // One 64-bit Complex-- 2 doubles-- pack<double, 2>::type == _m128d
            // Two 64-bit Complex-- 4 doubles-- pack<double, 4>::type == _m256d
            // Two 64-bit Complex-- 4 floats -- pack<float, 4>::type == _m128
            // Four 64-bit Complex-- 8 floats -- pack<float, 8>::type == _m256
            explicit Complex(F* const f): var(mm_load<F,L>::load(f)) {}

            explicit Complex(std::initializer_list<F> il): var(mm_load<F,L>::load(il.begin())) {};

            explicit Complex(typename pack<F,L>::type v): var(v) {}
            
            explicit Complex(F x, F y): var(mm_pair_set<F,L>::set(x, y)) {}

            explicit Complex(std::complex<F>& c): var(mm_pair_set<F,L>::set(c.real, c.imag)) {}

            explicit Complex(std::complex<F>* c): var(mm_complex_load<F,L>::load(c)) {}

            explicit Complex(std::initializer_list<std::complex<F>> il): var(mm_complex_load<F,L>::load(il.begin())) {};

            explicit Complex(): var(mm_zero<F,L>::get()) {}

            ///////////////////////////////////////////////////////////
            /* Basic operations with another pack of complex numbers */
            ///////////////////////////////////////////////////////////

            // Add with another pack of complex number
            Complex<F,L> operator+(Complex<F,L> const &o) {
                return Complex(mm_add(var, o.var));
            }

            // Subtract another pack of complex number
            Complex<F,L> operator-(Complex<F,L> const &o) {
                return Complex(mm_sub(var, o.var));
            }

            // Multiply by another pack of complex number
            Complex<F,L> operator*(Complex<F,L> const & o) {
                return Complex(mm_complex_mul(var, o.var));
            }

            // Divide by another pack of complex number
            Complex<F,L> operator/(Complex<F,L> const & o) {
                return Complex(mm_complex_div(var, o.var));
            }

            // Add with another complex number
            Complex<F,L> operator+=(Complex<F,L> const &o) {
                var = mm_add(var, o.var);
                return *this;
            }

            // Subtract another complex number from this
            Complex<F,L> operator-=(Complex<F,L> const &o) {
                var = mm_sub(var, o.var);
                return *this;
            }

            // Multiply by another complex number
            Complex<F,L> operator*=(Complex<F,L> const &o) {
                var = mm_complex_mul(var, o.var);
                return *this;
            }

            // Divide by another complex number
            Complex<F,L> operator/=(Complex<F,L> const &o) {
                var = mm_complex_div(var, o.var);
                return *this;
            }

            ///////////////////////////////////////////////////////////////
            /* Basic operations with a single real floating point number */
            ///////////////////////////////////////////////////////////////

            // Add with a floating point number
            Complex<F,L> operator+(F o) {
                return Complex(mm_add(var, mm_set1<F,L>::set(o)));
            }

            // Subtract a floating point number
            Complex<F,L> operator-(F o) {
                return Complex(mm_sub(var, mm_set1<F,L>::set(o)));
            }

            // Multiply by a floating point number
            Complex<F,L> operator*(F o) {
                return Complex(mm_mul(var, mm_set1<F,L>::set(o)));
            }

            // Divide by a floating point number
            Complex<F,L> operator/(F o) {
                return Complex(mm_div(var, mm_set1<F,L>::set(o)));
            }

            // Add with a floating point number
            Complex<F,L> operator+=(F o) {
                var = mm_add(var, mm_set1<F,L>::set(o));
                return *this;
            }

            // Subtract a floating point number
            Complex<F,L> operator-=(F o) {
                var = mm_sub(var, mm_set1<F,L>::set(o));
                return *this;
            }

            // Multiply by a floating point number
            Complex<F,L> operator*=(F o) {
                var = mm_mul(var, mm_set1<F,L>::set(o));
                return *this;
            }

            // Divide by a floating point number
            Complex<F,L> operator/=(F o) {
                var = mm_div(var, mm_set1<F,L>::set(o));
                return *this;
            }

            ///////////////////
            /* Other methods */
            ///////////////////

            // Store the modulus of the complex number in an array of size L/2
            void modulus(F* dest) {
                auto res = mm_complex_mod(var);
                for(int i = 0; i < L/2; i++) {
                    dest[i] = res[i*2];
                }
            }

            // Return the modulus of the complex number in a vector pack
            typename pack<F,L>::type modulus() {
                return mm_complex_mod(var);
            }

            // Conjugate the current complex number
            Complex<F,L> conjugate() {
                return Complex(mm_complex_conj(var));
            }

            // Store the Complex number in an array of length L
            void get(F* dest) {
                mm_store<F,L>::store(dest, var);
            }

            // Return a vector pack representation of this number
            typename pack<F,L>::type get() const {
                return var;
            }

        private:
            typename pack<F,L>::type var {};
    };

    // Determines how my complex number should be printed to an ostream
    template<typename F, int L>
    inline std::ostream& operator<<(std::ostream& os, const Complex<F,L>& dt){
        auto var = dt.get();
        os << "( ";
        for(int i = 0; i < L; i+=2) {
            os << var[i];
            if(var[i+1] < 0) os << " - " << -var[i+1] << "i";
            else os << " + " << var[i+1] << "i";
            if(i+2 < L) os << ", ";
        }
        os << " )";
        return os;
    }
}

#endif // FFTE_COMPLEX_HPP