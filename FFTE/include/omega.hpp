#ifndef FFTE_OMEGA_HPP
#define FFTE_OMEGA_HPP
/**
 * Code Author: Danny Sharp
 * This file is part of FFTE (Fast Fourier Transform Engine)
 */

#include "complex.hpp"
#include <vector>

/* The ability to cache exponential factors.
 * This isn't very well supported and is turned off,
 * as it actively hampers performance
 */
#define FFTE_CACHE_OMEGA 0

namespace FFTE {
#if FFTE_CACHE_OMEGA
    // Omega caching abilities
    class OmegaCache {
        private:
            size_t num;
            size_t denom;
            Complex inc_data;
            Complex curr_data;
        public:
            explicit OmegaCache() {
                init(0, Complex());
            }

            explicit OmegaCache(size_t denominator, Complex increment) {
                init(denominator, increment);
            }

            bool operator()(size_t numerator, size_t denominator) {
                if(denominator == denom && numerator == num+1) {
                    curr_data *= inc_data;
                    num++;
                    return true;
                }
                return false;
            }

            void init(size_t denominator, Complex increment) {
                num = 1;
                denom = denominator;
                inc_data = increment;
                curr_data = increment;
            }

            Complex getData() {
                return curr_data;
            }
    };
#endif
    // An enum for knowing whether the transform is forward or inverse
    enum class Direction {forward = -1, inverse = 1};

    // A class to hold the exponential factors of forward or inverse Fourier transforms
    class Omega { // exp(2pi i j / N)
        // For automatic freeing when we're done using this object
        using ComplexArr = std::vector<Complex>;

        private:

            size_t N; // Size of array
            ComplexArr data; // Where data is held

#if CACHE_OMEGA
            OmegaCache cache; // To hold the cache
#endif // END CACHE CHECKING
            // Log of N base 2
            int log2(size_t N) {
                int ret = 0;
                while((N>>ret) != 0) ret++;
                return ret;
            }

            // Initialize N and data (assumes dir is predetermined)
            void init(size_t length) {
                N = length;

                int arr_len = 2*log2(length);

                // Base the initialization on an exponential factor with
                // minimal increment
                double a = 2.*M_PI*(1.)/((double) N);
                Complex w0 {cos(a), static_cast<int>(dir) * sin(a)};

                // Allocate and initialize data
                // TODO: fix this accuracy stuff
                Complex w = w0;
                data = ComplexArr(arr_len);
                data[0] = Complex(1., 0.);

                for(size_t k = 1; k < N*N; k*=2) {
                    data[k] = Complex(cos(2*M_PI*((double) k)/N), static_cast<int>(dir) * sin(2*M_PI*((double) k)/N));
                }
            }

            /* Private function to get the appropriate exponential factor
             * for the given twiddle factor
             */
            Complex getElement(size_t s) {
                int k = 1;
                Complex ret = data[0];
                while(s != 0) {
                    if((s & 0x1)) {
                        ret *= data[k];
                    }
                    k++;
                    s >>= 1;
                }
                return ret;
            }

        public:
            // Public facing member that ultimately determines whether
            // the FFT is forward or inverse
            const Direction dir;

            // Don't create a default constructor.
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
                if(N % denom || num > denom*denom) {
                    std::cerr << "Invalid arguments for this Omega!\n";
                    exit(-1);
                }
#if FFTE_CACHE_OMEGA
                if(cache(num, denom)) {
                    return cache.getData();
                }
                else {
                    cache.init(denom, getElement(N/denom));
                }
#endif // END CACHE CHECKING
                size_t s = num * N / denom;
                
                return getElement(s);
            }

            // Constructs appropriate exponential factor if number > denom
            Complex operator()(size_t num1, size_t num2, size_t denom) {
                return this->operator()(num1*num2, denom);
            }

            // "inverts" the data
            Direction inv() {
                return (Direction) (-1*((int) dir));
            }

            // Gets the actual direction of the data
            Direction direction() {
                return dir;
            }
    };
}
#endif // FFTE_OMEGA_HPP