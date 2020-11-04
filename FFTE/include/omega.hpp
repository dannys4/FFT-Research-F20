#ifndef OMEGA_HPP
#define OMEGA_HPP
/**
 * Code Author: Danny Sharp
 * This file is part of FFTE (Fast Fourier Transform Engine)
 */

#include "complex.hpp"
#include <memory>

// An enum for knowing whether the transform is forward or inverse
enum class Direction {forward, inverse};

// A class to hold the exponential factors of forward or inverse Fourier transforms
class Omega { // exp(2pi i j / N)
    // For automatic freeing when we're done using this object
    using ComplexArr = std::shared_ptr<Complex[]>;

    private:
        uint64_t N; // Size of array
        ComplexArr data = nullptr; // Where data is held

        // Initialize N and data (assumes dir is predetermined)
        void init(uint64_t length) {
            N = length;

            // Base the initialization on an exponential factor with
            // minimal increment
            Complex w0;
            switch(dir) {
                case Direction::forward: w0 = initForward(); break;
                case Direction::inverse: w0 = initInverse(); break;
            }

            // Allocate and initialize data
            Complex w = w0;
            data = ComplexArr(new Complex[N]);
            data[0] = Complex(1., 0.);
            for(uint64_t k = 1; k < N; k++) {
                data[k] = w;
                w = w*w0;
            }
        }

        // Base initialization for a forward FFT
        Complex initForward() {
            double a = 2.*M_PI*(1.)/((double) N);
            Complex w0 {cos(a), -sin(a)};
            return w0;
        }

        // Base initialization for an inverse FFT
        Complex initInverse() {
            double a = 2.*M_PI*(1.)/((double) N);
            Complex w0 {cos(a), sin(a)};
            return w0;
        }

    public:
        // Public facing member that ultimately determines whether
        // the FFT is forward or inverse
        const Direction dir;

        Omega() = delete;

        // Basic constructor
        explicit Omega(Direction d): dir(d) {}

        // Constructor that also initializes the object
        explicit Omega(uint64_t length, Direction d): dir(d) {
            init(length);
        }

        // Constructor from other Omega
        explicit Omega(const Omega& o): N(o.N), data(o.data), dir(o.dir) {};

        // Copy assignment
        Omega& operator=(const Omega& o) {
            return *this;
        }

        // Initialize the Omega object post-construction
        void operator()(uint64_t length) {
            init(length);
        }

        // Return if this Omega object has any data
        bool operator()() {return data != nullptr;}

        // Return the appropriate twiddle, if possible
        Complex operator()(uint64_t num, uint64_t denom) {
            if(N % denom || num > denom) {
                throw "Invalid arguments for this Omega!\n";
            }
            uint64_t s = N / denom;
            return data[num*s];
        }

        #define MIN_INT(X, Y) (((X) < (Y)) ? (X) : (Y))
        #define MAX_INT(X, Y) (((X) > (Y)) ? (X) : (Y))

        Complex operator()(uint64_t num1, uint64_t num2, uint64_t denom) {
            if(num1*num2 < denom) return this->operator()(num1*num2, denom);
            
            uint64_t min = MIN_INT(num1, num2);
            uint64_t max = MAX_INT(num1, num2);
            
            Complex w0 = this->operator()(max, denom);
            Complex w = w0;
            for(uint64_t i = 1; i < min; i++) w = w0*w;
            return w;
        }
};

#endif // OMEGA_HPP