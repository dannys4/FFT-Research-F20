#include <cstddef>
#include <memory>
#include <vector>
#include <random>
#include <functional>
#include <iostream>
#include <chrono>
#include <stdlib.h>

#include <immintrin.h>
#define SPACE_NAME sig_space
#define REG_WIDTH 8

namespace SPACE_NAME {

    using double_avx = __m256d;

    template<typename T>
    class Aligned_Allocator{
        private:
            std::size_t alignment = 64/sizeof(T);
        public:
            using value_type            = T;
            using pointer               = T*;
            using const_pointer         = T const*;
            using void_pointer          = void*;
            using const_void_pointer    = void const*;
            using size_type             = std::size_t;
            using difference_type       = std::ptrdiff_t;

            template<typename U>
            struct rebind {
                typedef Aligned_Allocator<U> other;
            };

            Aligned_Allocator() noexcept {};
            Aligned_Allocator(const Aligned_Allocator&) noexcept {};
            bool operator==(const Aligned_Allocator&) const noexcept { return true; };
            bool operator!=(const Aligned_Allocator&) const noexcept { return false; };
            pointer allocate(size_type n, const_void_pointer=0) const {
                static_assert(64 % sizeof(T) == 0);
                return reinterpret_cast<pointer>(aligned_alloc(alignment*sizeof(value_type), (n+alignment - n % alignment)*sizeof(value_type)));
            }

            void deallocate(pointer p, size_type) const noexcept { free(p); }
    };

    

    class Signal
    {
    private:
        using avx_vec = std::vector<SPACE_NAME::double_avx, SPACE_NAME::Aligned_Allocator<SPACE_NAME::double_avx> >;
        avx_vec sig;
        uint8_t remainder;
        Signal(avx_vec arr, uint8_t rem) {
            sig = arr;
            remainder = rem;
        }
    public:
        Signal(const std::vector<double> input) {
            auto sz = input.size();
            sig = avx_vec();
            for(uint i = 0; i < input.size()/4; i++) {
                sig.push_back(_mm256_setr_pd(input[4*i], input[4*i+1], input[4*i+2], input[4*i+3]));
            }
            switch(input.size() % 4) {
                case 0:                                                                          remainder = 4; break;
                case 1: sig.push_back(_mm256_setr_pd(input[sz-1], 0          , 0          , 0)); remainder = 1; break;
                case 2: sig.push_back(_mm256_setr_pd(input[sz-2], input[sz-1], 0          , 0)); remainder = 2; break;
                case 3: sig.push_back(_mm256_setr_pd(input[sz-3], input[sz-2], input[sz-1], 0)); remainder = 3; break;
            }
        }


        double& operator[](std::size_t n) {
            return (reinterpret_cast<double*>(&sig[n/4]))[n%4];
        }

        Signal operator+(const Signal& other) {
            if(!(other.sig.size() == sig.size() && other.remainder == remainder)) exit(-1);
            avx_vec out {};
            for(uint i = 0; i < sig.size(); i++) {
                out.push_back(_mm256_add_pd(sig[i], other.sig[i]));
            }
            return Signal(out, remainder);
        }

        uint size() { return 4*(sig.size()-1) + remainder; }

    };
};



int main() {
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;
    typedef std::chrono::high_resolution_clock clock;

    std::vector<double> in1 {0.25, 0.50, 0.75, 1.00};
    std::vector<double> in2 {0.33, 0.67, 1.00, 1.33};
    std::vector<double> in3 {};
    std::vector<double> in4 {};
    auto sig1 = SPACE_NAME::Signal(in1);
    auto sig2 = SPACE_NAME::Signal(in2);

    auto rand = std::bind(std::uniform_real_distribution<>{0.0,10.0}, std::default_random_engine{});
    for(uint64_t i = 0; i < 1+5e7; i++) {
        in3.push_back(rand());
        in4.push_back(rand());
    }
    auto sig3 = SPACE_NAME::Signal(in3);
    auto sig4 = SPACE_NAME::Signal(in4);
    std::vector<double> compare{};
    std::cout << "sig3.size() = " << sig3.size() << '\n';
    auto start = clock::now();
    for(uint64_t i = 0; i < sig3.size(); i++) {
        compare.push_back(sig3[i]+sig4[i]);
    }
    auto end = clock::now();
    std::cout << "Standard addition took " << duration_cast<nanoseconds>(end-start).count() << "ns\n";
    start = clock::now();
    auto sig5 = sig3+sig4;
    end = clock::now();
    std::cout << "Vector addition took " << duration_cast<nanoseconds>(end-start).count() << "ns\n";
    for(uint64_t i = 0; i < sig5.size(); i++) {
        if(compare[i] != sig5[i]) {
            printf("THESE AREN'T EQUAL: compare[%d] = %f, sig5[%d] = %f\n", i, compare[i], i, sig5[i]);
            exit(-1);
        }
    }
}