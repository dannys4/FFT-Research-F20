#ifndef OMEGA_HPP
#define OMEGA_HPP
/**
 * Code Author: Danny Sharp
 * This file is part of FFTE (Fast Fourier Transform Engine)
 */

#include "complex.hpp"
#include <vector>

// An enum for knowing whether the transform is forward or inverse
enum class Direction {forward = -1, inverse = 1};

// A class to hold the exponential factors of forward or inverse Fourier transforms
class Omega { // exp(2pi i j / N)
    // For automatic freeing when we're done using this object
    using ComplexArr = std::vector<Complex>;

    private:
        size_t N; // Size of array
        ComplexArr data; // Where data is held

        // Initialize N and data (assumes dir is predetermined)
        void init(size_t length) {
            N = length;

            // Base the initialization on an exponential factor with
            // minimal increment
            double a = 2.*M_PI*(1.)/((double) N);
            Complex w0 {cos(a), static_cast<int>(dir) * sin(a)};

            // Allocate and initialize data
            Complex w = w0;
            data = ComplexArr(N);
            data[0] = Complex(1., 0.);
            for(size_t k = 1; k < N; k++) {
                data[k] = w;
                w = w*w0;
            }
        }

    public:
        // Public facing member that ultimately determines whether
        // the FFT is forward or inverse
        const Direction dir;

        Omega() = delete;

        // Basic constructor
        explicit Omega(Direction d): dir(d) {}

        // Constructor that also initializes the object
        explicit Omega(size_t length, Direction d): dir(d) {
            init(length);
        }

        // Constructor from other Omega
        explicit Omega(const Omega& o): N(o.N), data(o.data), dir(o.dir) {};

        // Copy assignment
        Omega& operator=(const Omega& o) {
            return *this;
        }

        // Initialize the Omega object post-construction
        void operator()(size_t length) {
            init(length);
        }

        // Return if this Omega object has any data
        bool operator()() {return data.size() != 0;}

        // Return the appropriate twiddle, if possible
        Complex operator()(size_t num, size_t denom) {
            if(N % denom || num > denom) {
                throw "Invalid arguments for this Omega!\n";
            }
            size_t s = N / denom;
            return data[num*s];
        }

        // Constructs appropriate exponential factor if number > denom
        Complex operator()(size_t num1, size_t num2, size_t denom) {
            if(num1*num2 < denom) return this->operator()(num1*num2, denom);
            
            size_t min = std::min<size_t>(num1, num2);
            size_t max = std::max<size_t>(num1, num2);
            
            Complex w0 = this->operator()(max, denom);
            Complex w = w0;
            for(size_t i = 1; i < min; i++) w = w0*w;
            return w;
        }
};

#endif // OMEGA_HPP