#ifndef FFTE_COMPLEX_HPP
#define FFTE_COMPLEX_HPP
#include <type_traits>
#include <iostream>
#include "vec_types.hpp"

namespace FFTE {
    template<typename F, class = typename std::enable_if<std::is_floating_point<F>::value, F>::type>
    class Twiddle {
        F re; F im;
    };

    template<typename F, int L>
    class Complex {
        public:
            // One 64-bit Complex-- 2 doubles-- pack<double, 2>::type == _m128d
            // Two 64-bit Complex-- 4 doubles-- pack<double, 4>::type == _m256d
            // Two 64-bit Complex-- 4 floats -- pack<float, 4>::type == _m128
            // Four 64-bit Complex-- 8 floats -- pack<float, 8>::type == _m256
            explicit Complex(F* const f): var(mm_load<F,L>::load(f)) {}

            explicit Complex(typename pack<F,L>::type v): var(v) {}
            
            explicit Complex(F x, F y): var(mm_pair_set<F,L>::set(x, y)) {}

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

            ///////////////////////////////////////////////////
            /* Basic operations with a single complex number */
            ///////////////////////////////////////////////////

            // Add with another pack of complex number
            Complex<F,L> operator+(Twiddle<F> const &o) {
                return Complex(mm_add(var, mm_pair_set<F, L>::set(o.re, o.im)));
            }

            // Subtract another pack of complex number
            Complex<F,L> operator-(Twiddle<F> const &o) {
                return Complex(mm_sub(var, mm_pair_set<F, L>::set(o.re, o.im)));
            }

            // Multiply by another pack of complex number
            Complex<F,L> operator*(Twiddle<F> const & o) {
                return Complex(mm_complex_mul(var, mm_pair_set<F, L>::set(o.re, o.im)));
            }

            // Divide by another pack of complex number
            Complex<F,L> operator/(Twiddle<F> const & o) {
                return Complex(mm_complex_div(var, mm_pair_set<F, L>::set(o.re, o.im)));
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
            if(var[i+1] < 0) os << " - " << -var[i+1];
            else os << " + " << var[i+1];
            if(i+2 < L) os << ", ";
        }
        os << " )";
        return os;
    }
}

#endif // FFTE_COMPLEX_HPP