#ifndef STOCK_FFT_ALLOCATOR_HPP
#define STOCK_FFT_ALLOCATOR_HPP

#include <memory>
#include <vector>
#include "complex.hpp"

namespace STOCK_FFT{
    template<typename F>
    class complex_allocator_t {
        public:
            // Taking care of mandatory allocator fields
            typedef F value_type;
            typedef F* pointer;
            typedef const F* const_pointer;
            typedef F& reference;
            typedef const F& const_reference;
            typedef size_t size_type;
            typedef ptrdiff_t difference_type;

            template <class U>
            struct rebind {
                typedef complex_allocator_t<U> other;
            };

            inline complex_allocator_t() noexcept {}
            inline complex_allocator_t(const complex_allocator_t&) noexcept {}

            template <class U>
            inline complex_allocator_t(const complex_allocator_t<U>&) noexcept {}

            inline ~complex_allocator_t() noexcept {}

            inline pointer address(reference r) { return &r; }
            inline const_pointer address(const_reference r) const { return &r; }

            pointer allocate(size_type n, typename std::allocator<void>::const_pointer hint = 0) {
                pointer ret = reinterpret_cast<pointer>(aligned_alloc(alignof(F), n*sizeof(F)));
                return ret;
            }
            inline void deallocate(pointer p, size_type) {
                return free(p);
            };

            inline void construct(pointer p, const_reference value) { new (p) value_type(value); }
            inline void destroy(pointer p) { p->~value_type(); }

            inline size_type max_size() const noexcept { return size_type(-1) / sizeof(F); }

            inline bool operator==(const complex_allocator_t&) { return true; }
            inline bool operator!=(const complex_allocator_t& rhs) { return !operator==(rhs); }
    };

    template<typename F, int L>
    using complex_allocator = complex_allocator_t<Complex<F, L>>;

    template<typename F, int L>
    using complex_vector = std::vector<Complex<F,L>, complex_allocator<F,L>>;
}
#endif // STOCK_FFT_ALLOCATOR_HPP