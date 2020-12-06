#include "test.hpp"
/**
 * Code Author: Danny Sharp
 * This file is part of FFTE (Fast Fourier Transform Engine)
 */

#define CHECKSUM_COMP 1
#define ENTRYWISE_COMP 0
#define FFT_LENGTH 5*27*64+1

/*
 * Functions to test veracity of outputs. These check against references
 * or expected values to get their answers, so sometimes (due to floating point
 * inaccuracies), there might be slight discrepencies
 */
using namespace FFTE;

// Checks how correct the forward Fourier transforms are
void check_fft(Direction dir) {
    switch(dir) {
        case Direction::forward: std::cout << "Checking the correctness of the FFT against the comparison... \n"; break;
        case Direction::inverse: std::cout << "Checking the correctness of the IFFT against the comparison...\n"; break;
    }

    const int n = FFT_LENGTH;
    int ell = getNumNodes(n);
    biFuncNode root[ell];
    init_fft_tree(root, n);

    auto in = (Complex*) malloc(n*sizeof(Complex));
    auto out_new = (Complex*) malloc(n*sizeof(Complex));
    auto out_comp = (Complex*) malloc(n*sizeof(Complex));
    auto out_ref = (Complex*) malloc(n*sizeof(Complex));
    for(int i = 0; i < n; i++) in[i] = Complex(i, 2.*i);

    reference_composite_FFT(n, in, out_ref, dir);
    Omega w(dir);
    root->fptr(in, out_new, 1, 1, root, w);

    
#if CHECKSUM_COMP
    std::cout << "Looking at the norm squared of all the errors...\n";
    double sum = (out_ref[0] - out_new[0]).modulus();
    for(int i = 1; i < n; i++) sum += (out_ref[i]-out_new[i]).modulus();
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
                       ", Err = " << (out_ref[i] - out_new[i]).modulus() << "\n";
    }
#endif
    free(in); free(out_comp); free(out_ref); free(out_new);
    switch(dir) {
        case Direction::forward: std::cout << "Done checking the FFT! \n\n";  break;
        case Direction::inverse: std::cout << "Done checking the IFFT!\n\n"; break;
    }
}

void check_omega() {
    std::cout << "Checking the correctness of the Omega class...\n";
    size_t N1 = 9;
    size_t N2 = 5;
    size_t N3 = 7;
    size_t N = N1*N2*N3;
    Omega w (N, Direction::forward);
    double sum = 0.;
    for(size_t i = 0; i < N; i++) {
        for(size_t j = 0; j < N; j++) sum += (w(i, j, N)-omega(i*j, N, Direction::forward)).modulus();
    }
    std::cout << "L2 Error in initialization with N = " << N << " is " << sum << "\n";
    
    sum = 0.;
    for(size_t i = 0; i < N1; i++) {
        for(size_t j = 0; j < N1; j++) sum += (w(i, j, N1)-omega(i*j, N1, Direction::forward)).modulus();
    }
    std::cout << "L2 Error in Checking against N1 = " << N1 << " is " << sum << "\n";
    
    sum = 0.;
    for(size_t i = 0; i < N2; i++) {
        for(size_t j = 0; j < N2; j++) sum += (w(i, j, N2)-omega(i*j, N2, Direction::forward)).modulus();
    }
    std::cout << "L2 Error in Checking against N2 = " << N2 << " is " << sum << "\n";
    
    sum = 0.;
    for(size_t i = 0; i < N3; i++) {
        for(size_t j = 0; j < N3; j++) {
            Complex ww = w(i,j, N3);
            Complex om = omega(i*j, N3, Direction::forward);

#if ENTRYWISE_COMP
            std::cout << "w(" << i << "," << j << "," << N3 << ") = " << ww << ", omega(" << i << "*" << j << ", " << N3 << ") = " << om << "\n";
#endif
            sum += (ww-om).modulus();
        }
    }
    std::cout << "L2 Error in Checking against N3 = " << N3 << " is " << sum << "\n";
    std::cout << "Done checking the Omega class!\n\n";
}

void check_fft_tree() {
    std::cout << "Checking the output of call graph node sizes against expected output...\n";
    size_t N = 15120;
    size_t ell = getNumNodes(N);
    biFuncNode root[ell];
    init_fft_tree(root, N);
    std::cout << "\t\t";
    printTree(root);
    std::cout << "\nExpected\t15120: (16, 945: (27, 35: (5, 7)))\n";
    std::cout << "Done checking the call graph!\n\n";
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
    biFuncNode root[ell];
    init_fft_tree(root, n);

    auto in = (Complex*) malloc(n*sizeof(Complex));
    auto out_ref = (Complex*) malloc(n*sizeof(Complex));
    auto out_new = (Complex*) malloc(n*sizeof(Complex));

    for(size_t i = 0; i < n; i++) in[i] = Complex(i, 2*i);
    Omega w(Direction::forward);

    int Ntrials = 10;

    auto start = clock::now();
    root->fptr(in, out_new, 1, 1, root, w);
    auto end = clock::now();
    auto new_fft = duration_cast<nanoseconds>(end-start).count();
    for(int i = 1; i < Ntrials; i++) {
        start = clock::now();
        root->fptr(in, out_new, 1, 1, root, w);
        end = clock::now();
        new_fft += duration_cast<nanoseconds>(end-start).count();
        if(in == (Complex*) 0x12345) std::cout << "this shouldn't print\n";
    }
    
    start = clock::now();
    reference_composite_FFT(n, in, out_new, Direction::forward);
    end = clock::now();
    auto ref_fft = duration_cast<nanoseconds>(end-start).count();
    for(int i = 1; i < Ntrials; i++) {
        start = clock::now();
        reference_composite_FFT(n, in, out_new, Direction::forward);
        end = clock::now();
        ref_fft += duration_cast<nanoseconds>(end-start).count();
        if(in == (Complex*) 0x12345) std::cout << "this shouldn't print\n";
    }
    std::cout << "The newest tranform takes " <<
                 (new_fft*1.e-9)/((double) Ntrials) <<
                 "s, where the reference takes " <<
                 (ref_fft*1e-9)/((double) Ntrials) << "s\n";
    
    free(in); free(out_ref); free(out_new);

    std::cout << "Done timing the FFTs!\n\n";
}


// Compares my complex multiplication to std::complex multiplication
void time_complex_mult() {
    std::cout << "Comparing the performance of my Complex multiplication vs. std::complex...\n";
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;
    typedef std::chrono::high_resolution_clock clock;

    std::vector<std::complex<double>> stdcomp1 {};
    std::vector<std::complex<double>> stdcomp2 {};
    std::vector<Complex> mycomp1 {};
    std::vector<Complex> mycomp2 {};
    auto rand = std::bind(std::uniform_real_distribution<>{0.0,10.0}, std::default_random_engine{});
    size_t len = 5e6;
    double a,b;
    for(size_t i = 0; i < len; i++) {
        a = rand(); b = rand();
        stdcomp1.push_back(std::complex<double>(a, b));
        mycomp1.push_back(Complex(a, b));
        a = rand(); b = rand();
        stdcomp2.push_back(std::complex<double>(a, b));
        mycomp2.push_back(Complex(a, b));
    }

    auto start = clock::now();
    for(size_t i = 0; i < len; i++) {
        auto m = stdcomp1[i] * stdcomp2[i];
        if((void*) (&m) == (void*) 0x123456) std::cout << "aaaa\n";
    }
    auto end = clock::now();
    auto stdmult = duration_cast<nanoseconds>(end-start).count();
    std::cout << "Standard mult took " << stdmult << "ns\n";
    
    start = clock::now();
    for(size_t i = 0; i < len; i++) {
        auto m = mycomp1[i] * mycomp2[i];
        if((void*) (&m) == (void*) 0x123456) std::cout << "aaaa\n";
    }
    end = clock::now();
    auto mymult = duration_cast<nanoseconds>(end-start).count();
    std::cout << "My mult took " << mymult << "ns\n";
    std::cout << "My mult was " << (((float) mymult/ (float) stdmult)) << " times the speed of std\n";
    std::cout << "Done timing complex multiplication!\n\n";
}

// Times how fast making a tree of uints at compile time is vs. runtime
void time_const_tree() {
    std::cout << "Time how long it takes to construct a const tree vs. a tree only known at runtime...\n";
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;
    typedef std::chrono::high_resolution_clock clock;

    size_t N = 1e5;
    std::cout << "Testing compile time performance...\n";
    auto start = clock::now();
    for(size_t i = 0; i < N; i++) {
        size_t ell = i*i - 3*i + 4;
        size_t k = getNumNodes(ell);
        constBiNode<size_t> root[k];
        initUintConstBiTree(root, ell);
        std::cout << "ell = " << ell << ", ";
        printRoot(root, k);
    }
    auto end = clock::now();
    std::cout << "Known at compile time took " << duration_cast<nanoseconds>(end-start).count() << "ns\n";
    std::cout << "Initializing random ints...\n";
    auto rand = std::bind(std::uniform_int_distribution<>{2, 4}, std::default_random_engine{});
    std::vector<uint> v {};
    for(size_t i = 0; i < N; i++) {
        v.push_back(rand());
    }
    std::cout << "Random ints initialized. Testing runtime performance...\n";
    start = clock::now();
    for(size_t i = 0; i < N; i++) {
        size_t ell = i*i - v[i]*i + 5;
        std::cout << "ell = " << ell << ", ";
        size_t k = getNumNodes(ell);
        constBiNode<size_t> root[k];
        initUintConstBiTree(root, ell);
        printRoot(root, k);
    }
    end = clock::now();
    std::cout << "Known at run time took " << duration_cast<nanoseconds>(end-start).count() << "ns\n";
    std::cout << "Done timing tree construction performance!\n\n";
}

// Times how fast the call to the function omega is
size_t omega_fcn_time(std::vector<size_t>& Nvec) {
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;
    typedef std::chrono::high_resolution_clock clock;

    auto start = clock::now();
    for(auto& Nj : Nvec) {
        for(size_t i = 0; i < Nj; i++) {
            for(size_t j = 0; j < Nj; j++) {
                Complex w0 = omega(i*j, Nj, Direction::forward);
                if(Nj == 16) std::cout << w0 << "\n";
            }
        }
    }
    auto end = clock::now();
    return duration_cast<nanoseconds>(end-start).count();
}

// Times how fast using the Omega class is
size_t omega_class_time(std::vector<size_t>& Nvec, Omega& w) {
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;
    typedef std::chrono::high_resolution_clock clock;

    auto start = clock::now();
    for(auto& Nj : Nvec) {
        for(size_t i = 0; i < Nj; i++) {
            for(size_t j = 0; j < Nj; j++) {
                Complex w0 = w(i, j, Nj);
                if(Nj == 16) std::cout << w0 << "\n";
            }
        }
    }
    auto end = clock::now();
    return duration_cast<nanoseconds>(end-start).count();
}

// Compares the time of the omega function to the omega class
void time_omega() {
    std::cout << "Comparing the time to use the Omega class vs. explicitly constructing omega at every iteration...\n";
    size_t N_factors = 10;
    int Njmax = 10;
    
    auto rand = std::bind(std::uniform_int_distribution<>{1,Njmax}, std::default_random_engine{});

    size_t Ntrials = 100;
    
    std::vector<size_t> Nvec {};
    size_t N = 1;
    for(size_t i = 0; i < N_factors; i++) {
        size_t tmp = rand();
        N *= tmp;
        Nvec.push_back(tmp);
    }
    Omega w (N, Direction::forward);
    auto fcn_time = omega_fcn_time(Nvec);
    auto class_time = omega_class_time(Nvec, w);

    for(size_t i = 0; i < Ntrials; i++) {
        Nvec.clear();
        N = 1;
        for(size_t i = 0; i < N_factors; i++) {
            size_t tmp = rand();
            N *= tmp;
            Nvec.push_back(tmp);
        }
        w(N);

        fcn_time += omega_fcn_time(Nvec);
        class_time += omega_class_time(Nvec, w);
    }

    std::cout << "Naive implementation elapsed time is " << (fcn_time/((double) Ntrials))*1.e-9 << "s\n";

    std::cout << "Class-based implementation of seconds is " << (class_time/((double) Ntrials))*1.e-9 << "s\n";

    std::cout << "Done timing the omega construction!\n\n";
}