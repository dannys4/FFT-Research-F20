/**
 * Intended to be used for FFTE. Collective inspiration and credit
 * to M. Stoyanov's LibHALA.
 * 
 * github.com/LIBHALA/hala/
 */

#ifndef VEC_TYPES_HPP
#define VEC_TYPES_HPP

#include <type_traits>
#include <immintrin.h>
#include <complex>

namespace FFTE {
    // Structs to figure out precise properties of inputs
    template<typename T>
    using is_float = std::is_same<float, typename std::remove_cv<T>::type>;
    template<typename T>
    using is_double = std::is_same<double, typename std::remove_cv<T>::type>;
    template<typename T>
    using is_fcomplex = std::is_same<std::complex<float>, typename std::remove_cv<T>::type>;
    template<typename T>
    using is_dcomplex = std::is_same<std::complex<double>, typename std::remove_cv<T>::type>;

    // Check if a type is real
    template<typename T> struct is_real {
        static constexpr bool value = is_float<T>::value || is_double<T>::value;
    };

    // Check if a type is complex
    template<typename T> struct is_complex {
        static constexpr bool value = is_fcomplex<T>::value || is_dcomplex<T>::value;
    };

    // Used for checking sizes of the packs
    template<int N> struct size_pack_type {};

    // Used for checking if it only uses one complex number
    template<int N, typename I = size_pack_type<1>> struct is_1_pack {};

    template<int N> struct is_1_pack<N, size_pack_type<N>> {
        typedef int type;
    };

    // Used for checking if it's two complex numbers
    template<int N, typename I = size_pack_type<2>> struct is_2_pack {};

    template<int N> struct is_2_pack<N, size_pack_type<N>> {
        typedef int type;
    };


    // Some simple operations that will be useful for vectorized types
    template<typename T, int N, typename = void, typename = void> struct mm_zero {
        static inline T get(){return static_cast<T>(0);}
    };

    template<typename T, int N, typename = void, typename = void> struct mm_load {
        static inline T load(T const *src) { return src[0]; }
    };

    template<typename T, int N, typename = void, typename = void> struct mm_store {
        static inline void store(T const *dest, T const &src) { dest[0] = src; }
    };

    // Real basic arithmetic for the "none" case
    template<typename T> inline T mm_fmadd(T const &a, T const &b, T const &c){ return a * b + c; }
    template<typename T> inline T mm_sqrt(T const &a){ return std::sqrt(a); }
    template<typename T> inline T mm_abs(T const &a){ return std::abs(a); }

    // Complex basic arithmetic for the "none" case
    template<typename T>
    inline typename std::enable_if<is_real<T>::value || is_complex<T>::value, T>::type mm_complex_mul(T const &a, T const &b){ return a * b; }
    template<typename T>
    inline typename std::enable_if<is_real<T>::value || is_complex<T>::value, T>::type mm_complex_fmadd(T const &a, T const &b, T const &c){ return a * b + c; }
    template<typename T>
    inline typename std::enable_if<is_real<T>::value || is_complex<T>::value, T>::type mm_complex_div(T const &a, T const &b){ return a / b; }
    template<typename T>
    inline typename std::enable_if<is_real<T>::value || is_complex<T>::value, T>::type mm_complex_conj(T const &a){ return hala_conj(a); }
    template<typename T>
    inline typename std::enable_if<is_real<T>::value || is_complex<T>::value, T>::type mm_complex_abs(T const &a){ return std::abs(a); }
    template<typename T>
    inline typename std::enable_if<is_real<T>::value || is_complex<T>::value, T>::type mm_complex_abs(std::complex<T> const &a){
        T r = std::abs(a);
        return {r, r};
    }

    // Setting the zero if there are two pairs of single precision complex numbers
    template<typename T, int N>
    struct mm_zero<T, N, typename std::enable_if<is_float<T>::value, T>::type, typename is_2_pack<N>::type> {
        static inline __m128 get(){ return _mm_setzero_ps(); }
    };

    // Loading from a pointer to at least 2 floats into a vectorized type
    template<typename T, int N>
    struct mm_load<T, N, typename std::enable_if<is_float<T>::value, T>::type, typename is_2_pack<N>::type> {
        static inline __m128 load(T const *src) { return __mm_load_ps(src); }
    };

    // Stores two floats from a vectorized type in a pointer of floats
    template<typename T, int N>
    struct mm_store<T, N, typename std::enable_if<is_float<T>::value, T>::type, typename is_2_pack<N>::type>  {
        static inline void store(T const *dest, __m128 const &src) { __mm_store_ps(dest, src); }
    };

    template<typename T>
    inline std::enable_if< mm_fmadd<__m128>(__m128 const &a, __m128 const &b, __m128 const &c){ return a * b + c; }
}

#endif // End VEC_TYPES_HPP

/**
 *  // Setting the zero if there is one pair of double precision complex numbers
 *  template<typename T, int N>
 *  struct mm_zero<T, N, typename std::enable_if<is_float<T>::value, T>::type, typename is_2_pack<N>::type> {
 *      static inline __m128d get(){ return _mm_setzero_pd(); }
 *  };
*/