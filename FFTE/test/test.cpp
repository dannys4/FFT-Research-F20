/**
 * Code Author: Danny Sharp
 * This file is part of the implementation for a stock FFT algorithm intended for HeFFTe
 */

#include "test.hpp"
#define CHECKSUM_COMP 1
#define ENTRYWISE_COMP 0
#define BATCH_PRINT 1
#define FFT_LENGTH 1024

/*
 * Functions to test veracity of outputs. These check against references
 * or expected values to get their answers, so sometimes (due to floating point
 * inaccuracies), there might be slight discrepencies
 */
using namespace stock_fft;

// Checks how correct the forward Fourier transforms are
void check_fft(Direction dir) {
    switch(dir) {
        case Direction::forward: std::cout << "Checking the correctness of the FFT against the comparison... \n"; break;
        case Direction::inverse: std::cout << "Checking the correctness of the IFFT against the comparison...\n"; break;
    }

    const int n = FFT_LENGTH;
    int ell = getNumNodes(n);
    biFuncNode<double, 2> root[ell];
    init_fft_tree(root, n);

    auto in = (Complex<double, 2>*) malloc(n*sizeof(Complex<double, 2>));
    auto out_new = (Complex<double, 2>*) malloc(n*sizeof(Complex<double, 2>));
    auto out_comp = (Complex<double, 2>*) malloc(n*sizeof(Complex<double, 2>));
    auto out_ref = (Complex<double, 2>*) malloc(n*sizeof(Complex<double, 2>));
    for(int i = 0; i < n; i++) in[i] = Complex<double, 2>(i, 2.*i);

    reference_composite_FFT(n, in, out_ref, dir);
    root->fptr(in, out_new, 1, 1, root, dir);

    
#if CHECKSUM_COMP
    std::cout << "Looking at the norm squared of all the errors...\n";
    double sum = (out_ref[0] - out_new[0]).modulus()[0];
    for(int i = 1; i < n; i++) sum += (out_ref[i]-out_new[i]).modulus()[0];
    std::cout << "Norm squared of Error: " << sum << "\n";
#endif
#if ENTRYWISE_COMP
    std::cout << "Comparing the values at each entry...\n";
    for(int i = 0; i < n; i++) {
            if(dir == Direction::inverse) {
                out_ref[i] /= (double) n;
                out_new[i] /= (double) n;
            }
            std::cout << "out_ref[" << i << "] = " << out_ref[i]  <<
                       ", out_new[" << i << "] = " << out_new[i]  <<
                       ", Err = " << (out_ref[i] - out_new[i]).modulus()[0] << "\n";
    }
#endif
    free(in); free(out_comp); free(out_ref); free(out_new);
    switch(dir) {
        case Direction::forward: std::cout << "Done checking the FFT! \n\n";  break;
        case Direction::inverse: std::cout << "Done checking the IFFT!\n\n"; break;
    }
}

// Checks how correct the forward Fourier transforms are
void check_fft_multidim(Direction dir) {
    switch(dir) {
        case Direction::forward: std::cout << "Checking the correctness of the multidimensional FFT against the comparison... \n"; break;
        case Direction::inverse: std::cout << "Checking the correctness of the multidimensional IFFT against the comparison...\n"; break;
    }

    const int n = FFT_LENGTH;
    int ell = getNumNodes(n);
    biFuncNode<double, 4> root[ell];
    init_fft_tree(root, n);

    auto al = alignof(Complex<double, 4>);

    auto in = (Complex<double, 4>*) aligned_alloc(al, n*sizeof(Complex<double, 4>));
    auto out_new = (Complex<double, 4>*) aligned_alloc(al, n*sizeof(Complex<double, 4>));
    auto out_comp = (Complex<double, 4>*) aligned_alloc(al, n*sizeof(Complex<double, 4>));
    auto out_ref = (Complex<double, 4>*) aligned_alloc(al, n*sizeof(Complex<double, 4>));
    double tmp_arr[4] = {0., 0., 0., 0.};
    for(int i = 0; i < n; i++) {
        in[i] = Complex<double, 4> (tmp_arr);
        tmp_arr[0] += 1; tmp_arr[1] += 2;
        tmp_arr[2] += 3; tmp_arr[3] += 4;
    }

    reference_composite_FFT(n, in, out_ref, dir);
    root->fptr(in, out_new, 1, 1, root, dir);

    
#if CHECKSUM_COMP
    std::cout << "Looking at the norm squared of all the errors...\n";
    auto tmp = (out_ref[0] - out_new[0]).modulus();
    double sum1 = tmp[0];
    double sum2 = tmp[2];
    for(int i = 1; i < n; i++) {
        tmp = (out_ref[i] - out_new[i]).modulus();
        sum1 += tmp[0];
        sum2 += tmp[2];
    }
    std::cout << "Norm squared of Error: (" << sum1 << ", " << sum2 << ")\n";
#endif
#if ENTRYWISE_COMP
    std::cout << "Comparing the values at each entry...\n";
    for(int i = 0; i < n; i++) {
            if(dir == Direction::inverse) {
                out_ref[i] /= (double) n;
                out_new[i] /= (double) n;
            }
            auto tmp = (out_ref[i] - out_new[i]).modulus();
            std::cout << "out_ref[" << i << "] = " << out_ref[i]  <<
                       ", out_new[" << i << "] = " << out_new[i]  <<
                       ", Err = (" << tmp[0] << ", " << tmp[1] << ")\n";
    }
#endif
    free(in); free(out_comp); free(out_ref); free(out_new);
    switch(dir) {
        case Direction::forward: std::cout << "Done checking the FFT! \n\n";  break;
        case Direction::inverse: std::cout << "Done checking the IFFT!\n\n"; break;
    }
}

// Checks an example of the function call tree developed
void check_fft_tree() {
    std::cout << "Checking the output of call graph node sizes against expected output...\n";
    size_t N = 15120;
    size_t ell = getNumNodes(N);
    biFuncNode<double, 2> root[ell];
    init_fft_tree(root, N);
    std::cout << "\t\t";
    printTree(root);
    std::cout << "\nExpected\t15120: (16, 945: (27, 35: (5, 7)))\n";
    std::cout << "Done checking the call graph!\n\n";
}

/* Functions to check all the arithmetic operations for complex types. 
 * Performs similar checks for float-4, float-8, double-2, and double-4.
 */

// Double-2 checks
void check_complex_ops_double_2() {
    std::cout << "Checking Complex Operations of Complex<double,2>\n\n";
    
    // Initial Values
    std::complex<double> x1 {1., 2.};
    std::complex<double> y1 {3., 4.};
    double alpha = 2.1;

    // Complex Objects
    Complex<double, 2> x {x1};
    Complex<double, 2> y {y1};

    // Test each out-of-place operation
    std::cout << x << "*" << y << " = " << x*y << ", ref = " << x1*y1 << "\n";
    std::cout << x << "+" << y << " = " << x+y << ", ref = " << x1+y1 << "\n";
    std::cout << x << "-" << y << " = " << x-y << ", ref = " << x1-y1 << "\n";
    std::cout << x << "/" << y << " = " << x/y << ", ref = " << x1/y1 << "\n";
    
    // Test each in place operation
    std::cout << x << "*=" << y << " ==> ";
    x *= y; x1 *= y1;
    std::cout << x << ", ref = " << x1 << "\n";

    std::cout << x << "+=" << y << " ==> ";
    x += y; x1 += y1;
    std::cout << x << ", ref = " << x1 << "\n";

    std::cout << x << "/=" << y << " ==> ";
    x /= y; x1 /= y1;
    std::cout << x << ", ref = " << x1 << "\n";

    std::cout << x << "-=" << y << " ==> ";
    x -= y; x1 -= y1;
    std::cout << x << ", ref = " << x1 << "\n";

    // Test each out-of-place operation with real number
    std::cout << x << "*" << alpha << " = " << x*alpha << ", ref = " << x1*alpha << "\n";
    std::cout << x << "+" << alpha << " = " << x+alpha << ", ref = " << x1+alpha << "\n";
    std::cout << x << "-" << alpha << " = " << x-alpha << ", ref = " << x1-alpha << "\n";
    std::cout << x << "/" << alpha << " = " << x/alpha << ", ref = " << x1/alpha << "\n";
    
    // Test each in-place operation with real number
    std::cout << x << "*=" << alpha << " ==> ";
    x *= alpha; x1 *= alpha;
    std::cout << x << ", ref = " << x1 << "\n";

    std::cout << x << "+=" << alpha << " ==> ";
    x += alpha; x1 += alpha;
    std::cout << x << ", ref = " << x1 << "\n";

    std::cout << x << "/=" << alpha << " ==> ";
    x /= alpha; x1 /= alpha;
    std::cout << x << ", ref = " << x1 << "\n";

    std::cout << x << "-=" << alpha << " ==> ";
    x -= alpha; x1 -= alpha;
    std::cout << x << ", ref = " << x1 << "\n";

    std::cout << "\nDone testing Complex<double,2>!\n\n";
}

// Double-4 checks
void check_complex_ops_double_4() {
    std::cout << "Checking Complex Operations of Complex<double,4>\n\n";
    
    // Initial Values
    std::complex<double> x1 {1., 2.};
    std::complex<double> y1 {3., 4.};
    std::complex<double> x2 {5., 6.};
    std::complex<double> y2 {7., 8.};
    double alpha = 2.1;

    // Complex Objects
    Complex<double, 4> x {x1, x2};
    Complex<double, 4> y {y1, y2};

    // Test each out-of-place operation
    std::cout << x << "*" << y << " = " << x*y << ", ref = " << x1*y1 << ", " << x2*y2 << "\n";
    std::cout << x << "+" << y << " = " << x+y << ", ref = " << x1+y1 << ", " << x2+y2 << "\n";
    std::cout << x << "-" << y << " = " << x-y << ", ref = " << x1-y1 << ", " << x2-y2 << "\n";
    std::cout << x << "/" << y << " = " << x/y << ", ref = " << x1/y1 << ", " << x2/y2 << "\n";
    
    // Test each in place operation
    std::cout << x << "*=" << y << " ==> ";
    x *= y; x1 *= y1; x2 *= y2;
    std::cout << x << ", ref = " << x1 << ", " << x2 << "\n";

    std::cout << x << "+=" << y << " ==> ";
    x += y; x1 += y1; x2 += y2;
    std::cout << x << ", ref = " << x1 << ", " << x2 << "\n";

    std::cout << x << "/=" << y << " ==> ";
    x /= y; x1 /= y1; x2 /= y2;
    std::cout << x << ", ref = " << x1 << ", " << x2 << "\n";

    std::cout << x << "-=" << y << " ==> ";
    x -= y; x1 -= y1; x2 -= y2;
    std::cout << x << ", ref = " << x1 << ", " << x2 << "\n";

    // Test each out-of-place operation with real number
    std::cout << x << "*" << alpha << " = " << x*alpha << ", ref = " << x1*alpha << ", " << x2*alpha << "\n";
    std::cout << x << "+" << alpha << " = " << x+alpha << ", ref = " << x1+alpha << ", " << x2+alpha << "\n";
    std::cout << x << "-" << alpha << " = " << x-alpha << ", ref = " << x1-alpha << ", " << x2-alpha << "\n";
    std::cout << x << "/" << alpha << " = " << x/alpha << ", ref = " << x1/alpha << ", " << x2/alpha << "\n";
    
    // Test each in-place operation with real number
    std::cout << x << "*=" << alpha << " ==> ";
    x *= alpha; x1 *= alpha; x2 *= alpha;
    std::cout << x << ", ref = " << x1 << ", " << x2 << "\n";

    std::cout << x << "+=" << alpha << " ==> ";
    x += alpha; x1 += alpha; x2 += alpha;
    std::cout << x << ", ref = " << x1 << ", " << x2 << "\n";

    std::cout << x << "/=" << alpha << " ==> ";
    x /= alpha; x1 /= alpha; x2 /= alpha;
    std::cout << x << ", ref = " << x1 << ", " << x2 << "\n";

    std::cout << x << "-=" << alpha << " ==> ";
    x -= alpha; x1 -= alpha; x2 -= alpha;
    std::cout << x << ", ref = " << x1 << ", " << x2 << "\n";

    std::cout << "\nDone testing Complex<double,4>!\n\n";
}

// Float-4 checks
void check_complex_ops_float_4() {
    std::cout << "Checking Complex Operations of Complex<float,4>\n\n";
    
    // Initial Values
    std::complex<float> x1 {1., 2.};
    std::complex<float> y1 {3., 4.};
    std::complex<float> x2 {5., 6.};
    std::complex<float> y2 {7., 8.};
    float alpha = 2.1f;

    // Complex Objects
    Complex<float, 4> x {x1, x2};
    Complex<float, 4> y {y1, y2};

    // Test each out-of-place operation
    std::cout << x << "*" << y << " = " << x*y << ", ref = " << x1*y1 << ", " << x2*y2 << "\n";
    std::cout << x << "+" << y << " = " << x+y << ", ref = " << x1+y1 << ", " << x2+y2 << "\n";
    std::cout << x << "-" << y << " = " << x-y << ", ref = " << x1-y1 << ", " << x2-y2 << "\n";
    std::cout << x << "/" << y << " = " << x/y << ", ref = " << x1/y1 << ", " << x2/y2 << "\n";
    
    // Test each in place operation
    std::cout << x << "*=" << y << " ==> ";
    x *= y; x1 *= y1; x2 *= y2;
    std::cout << x << ", ref = " << x1 << ", " << x2 << "\n";

    std::cout << x << "+=" << y << " ==> ";
    x += y; x1 += y1; x2 += y2;
    std::cout << x << ", ref = " << x1 << ", " << x2 << "\n";

    std::cout << x << "/=" << y << " ==> ";
    x /= y; x1 /= y1; x2 /= y2;
    std::cout << x << ", ref = " << x1 << ", " << x2 << "\n";

    std::cout << x << "-=" << y << " ==> ";
    x -= y; x1 -= y1; x2 -= y2;
    std::cout << x << ", ref = " << x1 << ", " << x2 << "\n";

    // Test each out-of-place operation with real number
    std::cout << x << "*" << alpha << " = " << x*alpha << ", ref = " << x1*alpha << ", " << x2*alpha << "\n";
    std::cout << x << "+" << alpha << " = " << x+alpha << ", ref = " << x1+alpha << ", " << x2+alpha << "\n";
    std::cout << x << "-" << alpha << " = " << x-alpha << ", ref = " << x1-alpha << ", " << x2-alpha << "\n";
    std::cout << x << "/" << alpha << " = " << x/alpha << ", ref = " << x1/alpha << ", " << x2/alpha << "\n";
    
    // Test each in-place operation with real number
    std::cout << x << "*=" << alpha << " ==> ";
    x *= alpha; x1 *= alpha; x2 *= alpha;
    std::cout << x << ", ref = " << x1 << ", " << x2 << "\n";

    std::cout << x << "+=" << alpha << " ==> ";
    x += alpha; x1 += alpha; x2 += alpha;
    std::cout << x << ", ref = " << x1 << ", " << x2 << "\n";

    std::cout << x << "/=" << alpha << " ==> ";
    x /= alpha; x1 /= alpha; x2 /= alpha;
    std::cout << x << ", ref = " << x1 << ", " << x2 << "\n";

    std::cout << x << "-=" << alpha << " ==> ";
    x -= alpha; x1 -= alpha; x2 -= alpha;
    std::cout << x << ", ref = " << x1 << ", " << x2 << "\n";

    std::cout << "\nDone testing Complex<float,4>!\n\n";
}

// Float-8 checks
void check_complex_ops_float_8() {
    std::cout << "Checking Complex Operations of Complex<float,8>\n\n";
    
    // Initial Values
    std::complex<float> y1 {3., 4.};
    std::complex<float> x1 {1., 2.};
    std::complex<float> x2 {5., 6.};
    std::complex<float> y2 {7., 8.};
    std::complex<float> x3 {9., 10.};
    std::complex<float> y3 {11., 12.};
    std::complex<float> x4 {13., 14.};
    std::complex<float> y4 {15., 16.};
    float alpha = 2.1f;

    // Complex Objects
    Complex<float, 8> x {x1, x2, x3, x4};
    Complex<float, 8> y {y1, y2, y3, y4};

    // Test each out-of-place operation
    std::cout << x << "*" << y << " = " << x*y << ", ref = " << x1*y1 << ", " << x2*y2 << ", " << x3*y3 << ", " << x4*y4 << "\n";
    std::cout << x << "+" << y << " = " << x+y << ", ref = " << x1+y1 << ", " << x2+y2 << ", " << x3+y3 << ", " << x4+y4 << "\n";
    std::cout << x << "-" << y << " = " << x-y << ", ref = " << x1-y1 << ", " << x2-y2 << ", " << x3-y3 << ", " << x4-y4 << "\n";
    std::cout << x << "/" << y << " = " << x/y << ", ref = " << x1/y1 << ", " << x2/y2 << ", " << x3/y3 << ", " << x4/y4 << "\n";
    
    // Test each in place operation
    std::cout << x << "*=" << y << " ==> ";
    x *= y; x1 *= y1; x2 *= y2; x3 *= y3; x4 *= y4;
    std::cout << x << ", ref = " << x1 << ", " << x2 << ", " << x3 << ", " << x4 << "\n";

    std::cout << x << "+=" << y << " ==> ";
    x += y; x1 += y1; x2 += y2; x3 += y3; x4 += y4;
    std::cout << x << ", ref = " << x1 << ", " << x2 << ", " << x3 << ", " << x4 << "\n";

    std::cout << x << "/=" << y << " ==> ";
    x /= y; x1 /= y1; x2 /= y2; x3 /= y3; x4 /= y4;
    std::cout << x << ", ref = " << x1 << ", " << x2 << ", " << x3 << ", " << x4 << "\n";

    std::cout << x << "-=" << y << " ==> ";
    x -= y; x1 -= y1; x2 -= y2; x3 -= y3; x4 -= y4;
    std::cout << x << ", ref = " << x1 << ", " << x2 << ", " << x3 << ", " << x4 << "\n";

    // Test each out-of-place operation with real number
    std::cout << x << "*" << alpha << " = " << x*alpha << ", ref = " << x1*alpha << ", " << x2*alpha << ", " << x3*alpha << ", " << x4*alpha << "\n";
    std::cout << x << "+" << alpha << " = " << x+alpha << ", ref = " << x1+alpha << ", " << x2+alpha << ", " << x3+alpha << ", " << x4+alpha << "\n";
    std::cout << x << "-" << alpha << " = " << x-alpha << ", ref = " << x1-alpha << ", " << x2-alpha << ", " << x3-alpha << ", " << x4-alpha << "\n";
    std::cout << x << "/" << alpha << " = " << x/alpha << ", ref = " << x1/alpha << ", " << x2/alpha << ", " << x3/alpha << ", " << x4/alpha << "\n";
    
    // Test each in-place operation with real number
    std::cout << x << "*=" << alpha << " ==> ";
    x *= alpha; x1 *= alpha; x2 *= alpha; x3 *= alpha; x4 *= alpha;
    std::cout << x << ", ref = " << x1 << ", " << x2 << ", " << x3 << ", " << x4 << "\n";

    std::cout << x << "+=" << alpha << " ==> ";
    x += alpha; x1 += alpha; x2 += alpha; x3 += alpha; x4 += alpha;
    std::cout << x << ", ref = " << x1 << ", " << x2 << ", " << x3 << ", " << x4 << "\n";

    std::cout << x << "/=" << alpha << " ==> ";
    x /= alpha; x1 /= alpha; x2 /= alpha; x3 /= alpha; x4 /= alpha;
    std::cout << x << ", ref = " << x1 << ", " << x2 << ", " << x3 << ", " << x4 << "\n";

    std::cout << x << "-=" << alpha << " ==> ";
    x -= alpha; x1 -= alpha; x2 -= alpha; x3 -= alpha; x4 -= alpha;
    std::cout << x << ", ref = " << x1 << ", " << x2 << ", " << x3 << ", " << x4 << "\n";

    std::cout << "\nDone testing Complex<float,8>!\n\n";
}

// Checking the arithmetic for every complex type.
void check_complex_ops() {
    check_complex_ops_double_2();
    check_complex_ops_double_4();
    check_complex_ops_float_4();
    check_complex_ops_float_8();
}

// Check the FFT 
void check_batch_fft() {
    typedef float T;
    const int P = 7;
    const int N = FFT_LENGTH;
    std::cout << "Checking the batch FFT on " << P << " vectors of length " << N << "...\n\n";
    STOCK_FFT::std_arrvec<T,P> in_d {};
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < P; j++) {
            auto t = std::complex<T>(2*j*i + i, 2*j*i + 2*i);
            in_d[j].push_back(t);
        }
    }
    auto out_d = batch_fft<P>(in_d, Direction::forward);
    #if BATCH_PRINT
    for(int i = 0; i < N; i++) {
        std::cout << out_d[0][i];
        for(int j = 1; j < P; j++) {
            std::cout << ", " << out_d[j][i];
        }
        std::cout << "\n";
    }
    #endif
    std::cout << "Done checking the batch FFT!\n\n";
}

/* Functions to clock my running times. These check in different cases
 * to see if what I'm doing is actually sufficiently fast
 */

// Compares a new FFT to the reference composite FFT
void time_fft() {
    std::cout << "Timing the current FFT tree configuration against reference composite FFT...\n";
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;
    typedef std::chrono::high_resolution_clock clock;

    const int n = FFT_LENGTH;
    int ell = getNumNodes(n);
    biFuncNode<double, 2> root[ell];
    init_fft_tree(root, n);

    auto in = (Complex<double, 2>*) malloc(n*sizeof(Complex<double, 2>));
    auto out_ref = (Complex<double, 2>*) malloc(n*sizeof(Complex<double, 2>));
    auto out_new = (Complex<double, 2>*) malloc(n*sizeof(Complex<double, 2>));

    for(size_t i = 0; i < n; i++) in[i] = Complex<double, 2>(i, 2*i);
    Direction dir = Direction::forward;

    int Ntrials = 10;

    auto start = clock::now();
    root->fptr(in, out_new, 1, 1, root, dir);
    auto end = clock::now();
    auto new_fft = duration_cast<nanoseconds>(end-start).count();
    for(int i = 1; i < Ntrials; i++) {
        start = clock::now();
        root->fptr(in, out_new, 1, 1, root, dir);
        end = clock::now();
        new_fft += duration_cast<nanoseconds>(end-start).count();
        if(in == (Complex<double, 2>*) 0x12345) std::cout << "this shouldn't print\n";
    }
#if 1
    start = clock::now();
    reference_composite_FFT(n, in, out_new, Direction::forward);
    end = clock::now();
    auto ref_fft = duration_cast<nanoseconds>(end-start).count();
    for(int i = 1; i < Ntrials; i++) {
        start = clock::now();
        reference_composite_FFT(n, in, out_new, Direction::forward);
        end = clock::now();
        ref_fft += duration_cast<nanoseconds>(end-start).count();
        if(in == (Complex<double, 2>*) 0x12345) std::cout << "this shouldn't print\n";
    }
    std::cout << "The newest tranform takes " <<
                 (new_fft*1.e-9)/((double) Ntrials) <<
                 "s, where the reference takes " <<
                 (ref_fft*1e-9)/((double) Ntrials) << "s\n";
#endif    
    free(in); free(out_ref); free(out_new);

    std::cout << "Done timing the FFTs!\n\n";
}

// Time the FFT 
void time_batch_fft() {
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;
    typedef std::chrono::high_resolution_clock clock;
    typedef float T;

    const int P = 17;
    const int N = FFT_LENGTH;
    std::cout << "Timing the batch FFT on " << P << " vectors of length " << N << "...\n\n";
    STOCK_FFT::std_arrvec<T,P> in_d {};
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < P; j++) {
            auto t = std::complex<T>(2*j*i + i, 2*j*i + 2*i);
            in_d[j].push_back(t);
        }
    }
    int N_warmup = 0;
    for(int i = 0; i < N_warmup; i++) batch_fft<P>(in_d, Direction::forward);

    auto start = clock::now();
    auto out_d = batch_fft<P>(in_d, Direction::forward);
    auto end = clock::now();
    std::cout << "The newest tranform takes " << (duration_cast<nanoseconds>(end-start).count()*1.e-9) << "\n";
    std::cout << "Done timing the batch FFT!\n\n";
}

void test_FFT_into_csv(std::string filename, int maxsize) {
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;
    typedef std::chrono::high_resolution_clock clock;
    auto rand = std::bind(std::uniform_real_distribution<>{1.001,1.075}, std::default_random_engine{});
    using F = double;
    const int L = 4;

    std::cout << "Timing the current FFT tree configuration against reference...\n";
    std::cout << "Outputting results into " << filename << "\n";
    std::ofstream myfile;

    Complex<F,L> *x = (Complex<F,L>*) aligned_alloc(alignof(Complex<F,L>), maxsize*sizeof(Complex<F,L>));
    Complex<F,L> *y = (Complex<F,L>*) aligned_alloc(alignof(Complex<F,L>), maxsize*sizeof(Complex<F,L>));
    Complex<F,L> *y_dft = (Complex<F,L>*) aligned_alloc(alignof(Complex<F,L>), maxsize*sizeof(Complex<F,L>));

    for(int i = 0; i < maxsize; i++) {
        auto tmp = std::array<std::complex<F>,L/2> {};
        for(int j = 0; j < L/2; j++) tmp[j] = std::complex<F> {(F) (2*j), (F) (2*j+1)};
        x[i] = Complex<F, L> {tmp.data()};
    }

    myfile.open(filename);
    myfile << "N,STOCK_FFT,DFT,err\n";
    
    const int iter_max = 200;
    auto dft_max = 1e9;
    int N = 10;
    bool do_dft = true;
    for(int j = 1; (j < iter_max) && (N < maxsize); j++) {
        int ell = getNumNodes(N);
        biFuncNode<F,L> root[ell];
        init_fft_tree(root, N);
        Direction dir = Direction::forward;

        auto start = clock::now();
        root->fptr(x, y, 1, 1, root, dir);
        auto end = clock::now();
        auto myTime = duration_cast<nanoseconds>(end-start).count();

        if(do_dft){
            start = clock::now();
            reference_DFT(N, x, y_dft, dir);
            end = clock::now();
        }
        auto dftTime = duration_cast<nanoseconds>(end-start).count();

        double err = 0.;
        if(do_dft){
            for(int j = 0; j < N; j++) {
                err += (y[j]-y_dft[j]).modulus()[0];
            }
        }

        myfile << N << "," << myTime;
        if(do_dft){
            myfile << "," << dftTime << "," << err << "\n";
        }
        else {
            myfile << ",0,0\n";
        }
        auto tmp = (int) (((double) N)*rand());
        N = (tmp > N) ? tmp : 2*N;
        if(do_dft && (dftTime > dft_max)) do_dft = false;
    }
    myfile.close();
    std::cout << "Finished!\n";
}

// Compares my complex multiplication to std::complex multiplication
std::pair<long long int,long long int> time_complex_mult(size_t len, bool at_runtime) {
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;
    typedef std::chrono::high_resolution_clock clock;
    using F = float;
    using s_cpx = std::complex<F>;
    const int n_cpx = 4;
    using f_cpx = Complex<F,2*n_cpx>;
    std::vector<std::array<s_cpx, n_cpx>> stdcomp1 {};
    std::vector<std::array<s_cpx, n_cpx>> stdcomp2 {};
    complex_vector<F,2*n_cpx> mycomp1 {};
    complex_vector<F,2*n_cpx> mycomp2 {};
    auto rand = std::bind(std::uniform_real_distribution<>{0.0,10.0}, std::default_random_engine{});
    for(size_t i = 0; i < len; i++) {
        std::array<s_cpx, n_cpx> tmp {};
        for(int j = 0; j < n_cpx; j++) tmp[j] = s_cpx(rand(), rand());
        stdcomp1.push_back(tmp);
        mycomp1.push_back(f_cpx(tmp.data()));
        for(int j = 0; j < n_cpx; j++) tmp[j] = s_cpx(rand(), rand());
        stdcomp2.push_back(tmp);
        mycomp2.push_back(f_cpx(tmp.data()));
    }

    auto start_std = clock::now();
    for(size_t i = 0; i < len; i++) {
        auto r = stdcomp1[i][0] * stdcomp2[i][0];
        for(int j = 1; j < n_cpx; j++) {
            r = stdcomp1[i][j] * stdcomp2[i][j];
        }
        if(at_runtime) std::cout << r << "\n";
    }
    auto end_std = clock::now();
    auto stdmult = duration_cast<nanoseconds>(end_std - start_std).count();
    
    auto start_my = clock::now();
    for(size_t i = 0; i < len; i++) {
        auto r = mycomp1[i] * mycomp2[i];
        if(at_runtime) std::cout << r << "\n";
    }
    auto end_my = clock::now();
    auto mymult = duration_cast<nanoseconds>(end_my-start_my).count();
    return std::pair<long long int, long long int> { stdmult, mymult };
}