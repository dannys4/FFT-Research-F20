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


    // Some simple operations that will be useful for vectorized types.
    // These cannot be in partially specialized template functions because
    // the only thing that changes in the function is the return type.

    template<typename F, int L> struct mm_zero {};
    template<typename F, int L> struct mm_load {};
    template<typename F, int L> struct mm_store {};
    template<typename F, int L> struct mm_pair_set{};
    template<typename F, int L> struct mm_set1{};
    template<typename F, int L> struct mm_complex_load{};

    /* Specializations for singular type */

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

    template<typename T>
    struct mm_set1<T,1> {
        static inline T set(T src) { return src; }
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

    //////////////////////////////////////////
    /* Below are structs for pack<float, 4> */
    //////////////////////////////////////////

    // Setting the zero if there are two pairs of single precision complex numbers
    template<>
    struct mm_zero<float, 4> {
        static inline pack<float, 4>::type get(){ return _mm_setzero_ps(); }
    };

    // Loading from a pointer to at least 4 floats into a vectorized type
    template<>
    struct mm_load<float, 4> {
        static inline pack<float, 4>::type load(float const *src) { return _mm_loadu_ps(src); }
    };

    // Stores four floats from a vectorized type in a pointer of floats
    template<>
    struct mm_store<float, 4>  {
        static inline void store(float *dest, pack<float, 4>::type const &src) { _mm_storeu_ps(dest, src); }
    };

    // Stores pair of floats into vectorized type
    template<>
    struct mm_pair_set<float, 4> {
        static inline pack<float, 4>::type set(float x, float y) { return _mm_setr_ps(x, y, x, y); }
    };

    // Set a vectorized type as a repeated float
    template<>
    struct mm_set1<float,4> {
        static inline pack<float, 4>::type set(float x) { return _mm_set1_ps(x); }
    };

    // Create a vectorized type from pointer to two 32-bit floating point complex numbers
    template<>
    struct mm_complex_load<float, 4> {
        static inline pack<float, 4>::type load(std::complex<float> const *src) {
            return _mm_setr_ps(src[0].real(), src[0].imag(), src[1].real(), src[1].imag());
        }
    };

    //////////////////////////////////////////
    /* Below are structs for pack<float, 8> */
    //////////////////////////////////////////

    // Setting the zero if there are four pairs of single precision complex numbers
    template<>
    struct mm_zero<float, 8> {
        static inline pack<float, 8>::type get(){ return _mm256_setzero_ps(); }
    };

    // Loading from a pointer to at least 8 floats into a vectorized type
    template<>
    struct mm_load<float, 8> {
        static inline pack<float, 8>::type load(float const *src) { return _mm256_loadu_ps(src); }
    };

    // Stores four floats from a vectorized type in a pointer of floats
    template<>
    struct mm_store<float, 8>  {
        static inline void store(float *dest, pack<float, 8>::type const &src) { _mm256_storeu_ps(dest, src); }
    };

    // Stores pair of floats into vectorized type
    template<>
    struct mm_pair_set<float, 8> {
        static inline pack<float, 8>::type set(float x, float y) { return _mm256_setr_ps(x, y, x, y, x, y, x, y); }
    };

    // Set a vectorized type as a repeated float
    template<>
    struct mm_set1<float,8> {
        static inline pack<float, 8>::type set(float x) { return _mm256_set1_ps(x); }
    };

    // Create a vectorized type from pointer to four 32-bit floating point complex numbers
    template<>
    struct mm_complex_load<float, 8> {
        static inline pack<float, 8>::type load(std::complex<float> const *src) {
            return _mm256_setr_ps(src[0].real(), src[0].imag(),
                                  src[1].real(), src[1].imag(),
                                  src[2].real(), src[2].imag(),
                                  src[3].real(), src[3].imag());
        }
    };

    ///////////////////////////////////////////
    /* Below are structs for pack<double, 2> */
    ///////////////////////////////////////////

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

    // Stores pair of doubles into vectorized type
    template<>
    struct mm_pair_set<double, 2> {
        static inline pack<double, 2>::type set(double x, double y) { return _mm_setr_pd(x, y); }
    };

    // Set a vectorized type as a repeated double
    template<>
    struct mm_set1<double, 2> {
        static inline pack<double, 2>::type set(double x) { return _mm_set1_pd(x); }
    };

    // Create a vectorized type from pointer to one 64-bit floating point complex number
    template<>
    struct mm_complex_load<double, 2> {
        static inline pack<double, 2>::type load(std::complex<double> const *src) {
            return _mm_setr_pd(src[0].real(), src[0].imag());
        }
    };

    ///////////////////////////////////////////
    /* Below are structs for pack<double, 4> */
    ///////////////////////////////////////////

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

    // Stores pair of floats into vectorized type
    template<>
    struct mm_pair_set<double, 4> {
        static inline pack<double, 4>::type set(double x, double y) { return _mm256_setr_pd(x, y, x, y); }
    };

    // Set a vectorized type as a repeated float
    template<>
    struct mm_set1<double, 4> {
        static inline pack<double, 4>::type set(double x) { return _mm256_set1_pd(x); }
    };

    // Create a vectorized type from pointer to two 64-bit floating point complex numbers
    template<>
    struct mm_complex_load<double, 4> {
        static inline pack<double, 4>::type load(std::complex<double> const *src) {
            return _mm256_setr_pd(src[0].real(), src[0].imag(), src[1].real(), src[1].imag());
        }
    };

    ///////////////////////////////////////////////////
    /* Elementary binary operations for vector packs */
    ///////////////////////////////////////////////////

    /* Addition */

    // Perform addition on vectorized packs of four floats
    inline pack<float, 4>::type mm_add(pack<float, 4>::type const &x,pack<float, 4>::type const &y) {
        return _mm_add_ps(x, y);
    }

    // Perform addition on vectorized packs of eight floats
    inline pack<float, 8>::type mm_add(pack<float, 8>::type const &x, pack<float, 8>::type const &y) {
        return _mm256_add_ps(x, y);
    }

    // Perform addition on vectorized packs of two doubles
    inline pack<double, 2>::type mm_add(pack<double, 2>::type const &x, pack<double, 2>::type const &y) {
        return _mm_add_pd(x, y);
    }

    // Perform addition on vectorized packs of four doubles
    inline pack<double, 4>::type mm_add(pack<double, 4>::type const &x, pack<double, 4>::type const &y) {
        return _mm256_add_pd(x, y);
    }

    /* Subtraction */

    // Perform subtraction on vectorized packs of four floats
    inline pack<float, 4>::type mm_sub(pack<float, 4>::type const &x,pack<float, 4>::type const &y) {
        return _mm_sub_ps(x, y);
    }

    // Perform subtraction on vectorized packs of eight floats
    inline pack<float, 8>::type mm_sub(pack<float, 8>::type const &x, pack<float, 8>::type const &y) {
        return _mm256_sub_ps(x, y);
    }

    // Perform subtraction on vectorized packs of two doubles
    inline pack<double, 2>::type mm_sub(pack<double, 2>::type const &x, pack<double, 2>::type const &y) {
        return _mm_sub_pd(x, y);
    }

    // Perform subtraction on vectorized packs of four doubles
    inline pack<double, 4>::type mm_sub(pack<double, 4>::type const &x, pack<double, 4>::type const &y) {
        return _mm256_sub_pd(x, y);
    }

    /* Multiplication */

    // Perform multiplication on vectorized packs of four floats
    inline pack<float, 4>::type mm_mul(pack<float, 4>::type const &x, pack<float, 4>::type const &y) {
        return _mm_mul_ps(x, y);
    }

    // Perform multiplication on vectorized packs of eight floats
    inline pack<float, 8>::type mm_mul(pack<float, 8>::type const &x, pack<float, 8>::type const &y) {
        return _mm256_mul_ps(x, y);
    }

    // Perform multiplication on vectorized packs of two doubles
    inline pack<double, 2>::type mm_mul(pack<double, 2>::type const &x, pack<double, 2>::type const &y) {
        return _mm_mul_pd(x, y);
    }

    // Perform multiplication on vectorized packs of four doubles
    inline pack<double, 4>::type mm_mul(pack<double, 4>::type const &x, pack<double, 4>::type const &y) {
        return _mm256_mul_pd(x, y);
    }

    /* Division */

    // Perform division on vectorized packs of four floats
    inline pack<float, 4>::type mm_div(pack<float, 4>::type const &x,pack<float, 4>::type const &y) {
        return _mm_mul_ps(x, y);
    }

    // Perform division on vectorized packs of eight floats
    inline pack<float, 8>::type mm_div(pack<float, 8>::type const &x, pack<float, 8>::type const &y) {
        return _mm256_div_ps(x, y);
    }

    // Perform division on vectorized packs of two doubles
    inline pack<double, 2>::type mm_div(pack<double, 2>::type const &x, pack<double, 2>::type const &y) {
        return _mm_div_pd(x, y);
    }

    // Perform division on vectorized packs of four doubles
    inline pack<double, 4>::type mm_div(pack<double, 4>::type const &x, pack<double, 4>::type const &y) {
        return _mm256_div_pd(x, y);
    }

    
    ///////////////////////////////////////////
    /* Complex operations using vector packs */
    ///////////////////////////////////////////
    
    // Complex Multiplication

    // Complex multiply two pairs of floats
    inline pack<float, 4>::type mm_complex_mul(pack<float, 4>::type const &x, pack<float, 4>::type const &y) {
        auto cc = _mm_permute_ps(y, 0b10100000);
        auto ba = _mm_permute_ps(x, 0b10110001);
        auto dd = _mm_permute_ps(y, 0b11110101);
        auto dba = _mm_mul_ps(ba, dd);
        auto mult = _mm_fmaddsub_ps(x, cc, dba);
        return mult;
    }

    // Complex multiply four pairs of floats
    inline pack<float, 8>::type mm_complex_mul(pack<float, 8>::type const &x, pack<float, 8>::type const &y) {
        auto cc = _mm256_permute_ps(y, 0b10100000);
        auto ba = _mm256_permute_ps(x, 0b10110001);
        auto dd = _mm256_permute_ps(y, 0b11110101);
        auto dba = _mm256_mul_ps(ba, dd);
        auto mult = _mm256_fmaddsub_ps(x, cc, dba);
        return mult;
    }
    
    // Complex multiply one pair of doubles
    inline pack<double, 2>::type mm_complex_mul(pack<double, 2>::type const &x, pack<double, 2>::type const &y) {
        auto cc = _mm_permute_pd(y, 0);
        auto ba = _mm_permute_pd(x, 0b01);
        auto dd = _mm_permute_pd(y, 0b11);
        auto dba = _mm_mul_pd(ba, dd);
        auto mult = _mm_fmaddsub_pd(x, cc, dba);
        return mult;
    }
    
    // Complex multiply two pairs of doubles
    inline pack<double, 4>::type mm_complex_mul(pack<double, 4>::type const &x, pack<double, 4>::type const &y) {
        auto cc = _mm256_permute_pd(y, 0b0000);
        auto ba = _mm256_permute_pd(x, 0b0101);
        auto dd = _mm256_permute_pd(y, 0b1111);
        auto dba = _mm256_mul_pd(ba, dd);
        auto mult = _mm256_fmaddsub_pd(x, cc, dba);
        return mult;
    }
    
    // Squared modulus of the complex numbers in a pack

    // Squared modulus of two 32-bit complex numbers in a pack
    inline pack<float, 4>::type mm_complex_sq_mod(pack<float, 4>::type const &x) {
        return _mm_or_ps(_mm_dp_ps(x, x, 0b11001100), _mm_dp_ps(x, x, 0b00110011));
    }

    // Squared modulus of four 32-bit complex numbers in a pack
    inline pack<float, 8>::type mm_complex_sq_mod(pack<float, 8>::type const &x) {
        return _mm256_or_ps(_mm256_dp_ps(x, x, 0b11001100), _mm256_dp_ps(x, x, 0b00110011));
    }

    // Squared modulus of one 64-bit complex number in a pack
    inline pack<double, 2>::type mm_complex_sq_mod(pack<double, 2>::type const &x) {
        return _mm_dp_pd(x, x, 0b11111111);
    }

    // Squared modulus of two 64-bit complex numbers in a pack
    inline pack<double, 4>::type mm_complex_sq_mod(pack<double, 4>::type const &x) {
        auto a = _mm256_mul_pd(x, x);
        return _mm256_hadd_pd(a, a);
    }

    // Moduli (with square root) of complex numbers

    // Moduli of two 32-bit complex numbers in a pack
    inline pack<float, 4>::type mm_complex_mod(pack<float, 4>::type const &x) {
        return _mm_sqrt_ps(mm_complex_sq_mod(x));
    }

    // Moduli of four 32-bit complex numbers in a pack
    inline pack<float, 8>::type mm_complex_mod(pack<float, 8>::type const &x) {
        return _mm256_sqrt_ps(mm_complex_sq_mod(x));
    }
    
    // Modulus of one 64-bit complex number in a pack
    inline pack<double, 2>::type mm_complex_mod(pack<double, 2>::type const &x) {
        return _mm_sqrt_pd(mm_complex_sq_mod(x));
    }
    
    // Moduli of two 64-bit complex numbers in a pack
    inline pack<double, 4>::type mm_complex_mod(pack<double, 4>::type const &x) {
        return _mm256_sqrt_pd(mm_complex_sq_mod(x));
    }

    // Conjugation of two 32-bit complex numbers
    inline pack<float, 4>::type mm_complex_conj(pack<float, 4>::type const &x) {
        return _mm_blend_ps(x, -x, 0b1010);
    }

    // Conjugation of four 32-bit complex numbers
    inline pack<float, 8>::type mm_complex_conj(pack<float, 8>::type const &x) {
        return _mm256_blend_ps(x, -x, 0b10101010);
    }
    
    // Conjugation of two 32-bit complex numbers
    inline pack<double, 2>::type mm_complex_conj(pack<double, 2>::type const &x) {
        return _mm_blend_pd(x, -x, 0b10);
    }

    // Conjugation of four 32-bit complex numbers
    inline pack<double, 4>::type mm_complex_conj(pack<double, 4>::type const &x) {
        return _mm256_blend_pd(x, -x, 0b1010);
    }
    
    // Complex division

    // Divide x by y, where x and y are each 2 32-bit complex numbers
    inline pack<float, 4>::type mm_complex_div(pack<float, 4>::type const &x, pack<float, 4>::type const &y) {
        return _mm_div_ps(mm_complex_mul(x, mm_complex_conj(y)), mm_complex_sq_mod(y));
    }
    
    // Divide x by y, where x and y are each 4 32-bit complex numbers
    inline pack<float, 8>::type mm_complex_div(pack<float, 8>::type const &x, pack<float, 8>::type const &y) {
        return _mm256_div_ps(mm_complex_mul(x, mm_complex_conj(y)), mm_complex_sq_mod(y));
    }
    
    // Divide x by y, where x and y are each 1 64-bit complex number
    inline pack<double, 2>::type mm_complex_div(pack<double, 2>::type const &x, pack<double, 2>::type const &y) {
        return _mm_div_pd(mm_complex_mul(x, mm_complex_conj(y)), mm_complex_sq_mod(y));
    }
    
    // Divide x by y, where x and y are each 2 64-bit complex numbers
    inline pack<double, 4>::type mm_complex_div(pack<double, 4>::type const &x, pack<double, 4>::type const &y) {
        return _mm256_div_pd(mm_complex_mul(x, mm_complex_conj(y)), mm_complex_sq_mod(y));
    }

}

#endif