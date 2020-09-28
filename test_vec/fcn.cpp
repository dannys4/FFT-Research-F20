#include <cstddef>
#include <memory>
#include <vector>
#include <stdlib.h>

#include <immintrin.h>
#define SPACE_NAME sig_space
#define REG_WIDTH 8
using double_avx = __m256d;

namespace SPACE_NAME {


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
    using avx_vec = std::vector<double_avx, SPACE_NAME::Aligned_Allocator<double_avx>>;
    private:
        avx_vec sig;
        uint8_t remainder;
    public:

        Signal(const std::vector<double> input) {
            auto sz = input.size();
            sig = avx_vec();
            for(uint i = 0; i < input.size()/4; i++) {
                sig.push_back(_mm256_setr_pd(input[4*i], input[4*i+1], input[4*i+2], input[4*i+3]));
            }
            switch(input.size() % 4) {
                case 0:                                                                       remainder = 4; break;
                case 1: sig[sz/4] = _mm256_setr_pd(input[sz-1], 0          , 0          , 0); remainder = 1; break;
                case 2: sig[sz/4] = _mm256_setr_pd(input[sz-2], input[sz-1], 0          , 0); remainder = 2; break;
                case 3: sig[sz/4] = _mm256_setr_pd(input[sz-3], input[sz-2], input[sz-1], 0); remainder = 3; break;
            }
        }

        Signal(avx_vec arr) {
            sig = arr;
        }

        double operator[](std::size_t n) {
            return ((double*) &sig[n/4])[n%4];
        }

        Signal operator+(const Signal& other) {
            if(other.sig.size() != sig.size()) exit(-1);
            avx_vec out(sig.size() / 4 + (remainder != 4));
            for(uint i = 0; i < sig.size(); i++) {
                out.push_back(_mm256_add_pd(sig[i], other.sig[i]));
            }
            return Signal(out);
        }

        uint size() { return 4*(sig.size()-1) + remainder; }

    };
};



int main() {
    std::vector<double> in1 {0.25, 0.50, 0.75, 1.00, 1.25, 1.50};
    std::vector<double> in2 {0.33, 0.67, 1.00, 1.33, 1.67, 2.00};
    auto sig1 = sig_space::Signal(in1);

    for(uint i = 0; i < in1.size(); i++) {
        printf("%f\n", sig1[i]);
    }
}