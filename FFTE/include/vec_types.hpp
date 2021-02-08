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

    // Used for checking if it only uses one number
    template<int N, typename I = size_pack_type<1>> struct is_1_pack {};

    template<int N> struct is_1_pack<N, size_pack_type<N>> {
        typedef int type;
        static constexpr bool value = true;
    };

    // Used for checking if it's two numbers
    template<int N, typename I = size_pack_type<2>> struct is_2_pack {};

    template<int N> struct is_2_pack<N, size_pack_type<N>> {
        typedef int type;
        static constexpr bool value = true;
    };

    // Used for checking if it's four numbers
    template<int N, typename I = size_pack_type<4>> struct is_4_pack {};

    template<int N> struct is_4_pack<N, size_pack_type<N>> {
        typedef int type;
        static constexpr bool value = true;
    };

    // Used for getting vector pack types
    template<typename T, int N> struct pack {};
    template<> struct pack<float, 4> { typedef __m128 type; };
    template<> struct pack<double, 2> { typedef __m128d type; };
    template<> struct pack<double, 4> { typedef __m256d type; };
    template<> struct pack<float, 8> { typedef __m256 type; };

    template<typename T, int N> struct mm_zero {};
    template<typename T, int N> struct mm_load {};
    template<typename T, int N> struct mm_store {};

    // Some simple operations that will be useful for vectorized types
    template<typename T>
    struct mm_zero<T, 1> {
        static inline T get(){return static_cast<T>(0);}
    };

    template<typename T>
    struct mm_load<T, 1> {
        static inline T load(T const *src) { return src[0]; }
    };

    template<typename T>
    struct mm_store<T, 1> {
        static inline void store(T const *dest, T const &src) { dest[0] = src; }
    };

    // Real basic arithmetic for the "none" case
    template<typename T>
    inline T mm_fmadd(T const &a, T const &b, T const &c){ return a * b + c; }
    template<typename T>
    inline T mm_sqrt(T const &a){ return std::sqrt(a); }
    template<typename T>
    inline T mm_abs(T const &a){ return std::abs(a); }

    // Complex basic arithmetic for the "none" case
    template<typename T>
    inline typename std::enable_if<is_real<T>::value || is_complex<T>::value, T>::type mm_complex_mul(T const &a, T const &b){ return a * b; }
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
    template<>
    struct mm_zero<float, 4> {
        static inline __m128 get(){ return _mm_setzero_ps(); }
    };

    // Loading from a pointer to at least 4 floats into a vectorized type
    template<>
    struct mm_load<float, 4> {
        static inline __m128 load(float const *src) { return _mm_loadu_ps(src); }
    };

    // Stores four floats from a vectorized type in a pointer of floats
    template<>
    struct mm_store<float, 4>  {
        static inline void store(float *dest, __m128 const &src) { _mm_storeu_ps(dest, src); }
    };

    // Setting the zero if there is one pair of double precision complex numbers
    template<>
    struct mm_zero<double, 2> {
        static inline __m128d get(){ return _mm_setzero_pd(); }
    };

    // Loading from a pointer to at least 2 doubles into a vectorized type
    template<>
    struct mm_load<double, 2> {
        static inline __m128d load(double const *src) { return _mm_loadu_pd(src); }
    };

    // Stores two doubles from a vectorized type in a pointer of doubles
    template<>
    struct mm_store<double, 2>  {
        static inline void store(double *dest, __m128d const &src) { _mm_storeu_pd(dest, src); }
    };

    // Setting the zero if there is one pair of double precision complex numbers
    template<>
    struct mm_zero<double, 4> {
        static inline __m256d get(){ return _mm256_setzero_pd(); }
    };

    // Loading from a pointer to at least 2 doubles into a vectorized type
    template<>
    struct mm_load<double, 4> {
        static inline __m256d load(double const *src) { return _mm256_loadu_pd(src); }
    };

    // Stores two doubles from a vectorized type in a pointer of doubles
    template<>
    struct mm_store<double, 4>  {
        static inline void store(double *dest, __m256d const &src) { _mm256_storeu_pd(dest, src); }
    };

    // Complex multiplication using vector packs
    
    // Two pairs of floats
    inline pack<float, 4>::type mm_complex_mul(pack<float, 4>::type const &x, pack<float, 4>::type const &y) {
        auto cc = _mm_permute_ps(y, 0b10100000);
        auto ba = _mm_permute_ps(x, 0b10110001);
        auto dd = _mm_permute_ps(y, 0b11110101);
        auto dba = _mm_mul_ps(ba, dd);
        auto mult = _mm_fmaddsub_ps(x, cc, dba);
        return mult;
    }

    // Four pairs of floats
    inline pack<float, 8>::type mm_complex_mul(pack<float, 8>::type const &x, pack<float, 8>::type const &y) {
        auto cc = _mm256_permute_ps(y, 0b10100000);
        auto ba = _mm256_permute_ps(x, 0b10110001);
        auto dd = _mm256_permute_ps(y, 0b11110101);
        auto dba = _mm256_mul_ps(ba, dd);
        auto mult = _mm256_fmaddsub_ps(x, cc, dba);
        return mult;
    }
    
    // Pair of doubles
    inline pack<double, 2>::type mm_complex_mul(pack<double, 2>::type const &x, pack<double, 2>::type const &y) {
        auto cc = _mm_permute_pd(y, 0);
        auto ba = _mm_permute_pd(x, 0b01);
        auto dd = _mm_permute_pd(y, 0b11);
        auto dba = _mm_mul_pd(ba, dd);
        auto mult = _mm_fmaddsub_pd(x, cc, dba);
        return mult;
    }
    
    // Two pairs of doubles
    inline pack<double, 4>::type mm_complex_mul(pack<double, 4>::type const &x, pack<double, 4>::type const &y) {
        auto cc = _mm256_permute_pd(y, 0b0000);
        auto ba = _mm256_permute_pd(x, 0b0101);
        auto dd = _mm256_permute_pd(y, 0b1111);
        auto dba = _mm256_mul_pd(ba, dd);
        auto mult = _mm256_fmaddsub_pd(x, cc, dba);
        return mult;
    }

    // Addition using vector packs

    // Four floats
    inline pack<float, 4>::type mm_add(pack<float, 4>::type const &x,pack<float, 4>::type const &y) {
        return _mm_add_ps(x, y);
    }

    // Eight floats
    inline pack<float, 8>::type mm_add(pack<float, 8>::type const &x, pack<float, 8>::type const &y) {
        return _mm256_add_ps(x, y);
    }

    // Two doubles
    inline pack<double, 2>::type mm_add(pack<double, 2>::type const &x, pack<double, 2>::type const &y) {
        return _mm_add_pd(x, y);
    }

    // Four doubles
    inline pack<double, 4>::type mm_add(pack<double, 4>::type const &x, pack<double, 4>::type const &y) {
        return _mm256_add_pd(x, y);
    }

    // Subtraction using vector packs

    // Four floats
    inline pack<float, 4>::type mm_sub(pack<float, 4>::type const &x,pack<float, 4>::type const &y) {
        return _mm_sub_ps(x, y);
    }

    // Eight floats
    inline pack<float, 8>::type mm_sub(pack<float, 8>::type const &x, pack<float, 8>::type const &y) {
        return _mm256_sub_ps(x, y);
    }

    // Two doubles
    inline pack<double, 2>::type mm_sub(pack<double, 2>::type const &x, pack<double, 2>::type const &y) {
        return _mm_sub_pd(x, y);
    }

    // Four doubles
    inline pack<double, 4>::type mm_sub(pack<double, 4>::type const &x, pack<double, 4>::type const &y) {
        return _mm256_sub_pd(x, y);
    }

}

#endif