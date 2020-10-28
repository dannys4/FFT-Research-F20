#ifndef OMEGA_HPP
#define OMEGA_HPP

#include "complex.hpp"

enum Direction {forward, inverse};

class Omega {
    private:
        uint64_t N;
        Complex* data;

        Complex initForward() {
            double a = 2.*M_PI*(1.)/((double) N);
            Complex w0 {cos(a), -sin(a)};
            return w0;
        }
        Complex initInverse() {
            double a = 2.*M_PI*(1.)/((double) N);
            Complex w0 {cos(a), sin(a)};
            return w0;
        }
    public:
        Omega(uint64_t length, Direction d) {
            N = length;
            Complex w0;
            switch(d) {
                case forward: w0 = initForward(); break;
                case inverse: w0 = initInverse(); break;
            }
            Complex w = w0;
            data = (Complex*) malloc(N*sizeof(Complex));
            for(int k = 0; k < N; k++) {
                data[k] = w;
                w = w*w0;
            }
        }

        ~Omega() {
            free(data);
        }

        Complex operator()(uint64_t num, uint64_t denom) {
            if(N % denom) {
                throw "Cannot use this object for this denominator!\n";
            }
            uint64_t s = N / denom;
            return data[num*s];
        }
};

#endif // OMEGA_HPP