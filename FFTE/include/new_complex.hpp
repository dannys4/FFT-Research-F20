#include <type_traits>
#include <iostream>
#include "vec_types.hpp"

namespace FFTE {
    template<typename F, int N, class = typename std::enable_if<std::is_floating_point<F>::value, F>::type>
    class NewComplex {
        public:
            explicit NewComplex(F* const f) {
                // 2 doubles-- pack<double, 2>::type == _m128d
                // 4 doubles-- pack<double, 4>::type == _m256d
                // 4 floats -- pack<float, 4>::type == _m128
                typename pack<F, N>::type v = mm_load<F,N>::load(f);
                var = v;
            }

            explicit NewComplex(typename pack<F, N>::type v): var(v) {}

            NewComplex<F,N> operator*(NewComplex<F,N> const &o) {
                return NewComplex(mm_complex_mul(var, o.var));
            }

            void getComplex(F* dest) {
                mm_store<F, N>::store(dest, var);
            }

        private:
            typename pack<F, N>::type var {};
    };
}