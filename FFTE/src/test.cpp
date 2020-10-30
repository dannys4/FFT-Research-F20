#include "test.hpp"

#define CHECKSUM_COMP 1
#define ENTRYWISE_COMP 0
#define FFT_LENGTH 27*128*35

/*
 * Functions to test veracity of outputs. These check against references
 * or expected values to get their answers, so sometimes (due to floating point
 * inaccuracies), there might be slight discrepencies
 */
void check_fft() {
    const int n = 25;
    int ell = getNumNodes(n);
    constBiFuncNode root[ell];
    init_fft_tree(root, n);

    auto in = (Complex*) malloc(n*sizeof(Complex));
    auto out_new = (Complex*) malloc(n*sizeof(Complex));
    auto out_comp = (Complex*) malloc(n*sizeof(Complex));
    auto out_ref = (Complex*) malloc(n*sizeof(Complex));
    for(int i = 0; i < n; i++) in[i] = Complex(1., 2.);

    reference_DFT(n, in, out_ref);
    Omega w(Direction::forward);
    root->fptr(in, out_new, 1, 1, root, w);

    
#if CHECKSUM_COMP
    double sum = (out_ref[0] - out_new[0]).modulus();
    for(int i = 1; i < n; i++) sum += (out_ref[i]-out_new[i]).modulus();
    std::cout << "Norm of Error: " << sum << "\n";
#endif
#if ENTRYWISE_COMP
    for(int i = 0; i < n; i++) std::cout << "out_ref[" << i << "] = " << out_ref[i] << ", out_new[" << i << "] = " << out_new[i] << ", Err = " << (out_ref[i] - out_new[i]).modulus() << "\n";
#endif
    free(in); free(out_comp); free(out_ref); free(out_new);
}

void check_omega() {
    uint64_t N1 = 27;
    uint64_t N2 = 5;
    uint64_t N3 = 7;
    uint64_t N = N1*N2*N3;
    Omega w (N, Direction::forward);
    double sum = 0.;
    for(uint64_t i = 0; i < N; i++) {
        for(uint64_t j = 0; j < N; j++) sum += (w(i, j, N)-omega(i*j, N)).modulus();
    }
    std::cout << "L2 Error in initialization with N = " << N << " is " << sum << "\n";
    
    sum = 0.;
    for(uint64_t i = 0; i < N1; i++) {
        for(uint64_t j = 0; j < N1; j++) sum += (w(i, j, N1)-omega(i*j, N1)).modulus();
    }
    std::cout << "L2 Error in Checking against N1 = " << N1 << " is " << sum << "\n";
    
    sum = 0.;
    for(uint64_t i = 0; i < N2; i++) {
        for(uint64_t j = 0; j < N2; j++) sum += (w(i, j, N2)-omega(i*j, N2)).modulus();
    }
    std::cout << "L2 Error in Checking against N2 = " << N2 << " is " << sum << "\n";
    
    sum = 0.;
    for(uint64_t i = 0; i < N3; i++) {
        for(uint64_t j = 0; j < N3; j++) {
            Complex ww = w(i,j, N3);
            Complex om = omega(i*j, N3);

#if ENTRYWISE_COMP
            std::cout << "w(" << i << "," << j << "," << N3 << ") = " << ww << ", omega(" << i << "*" << j << ", " << N3 << ") = " << om << "\n";
#endif
            sum += (ww-om).modulus();
        }
    }
    std::cout << "L2 Error in Checking against N3 = " << N3 << " is " << sum << "\n";
}

void check_fft_tree() {
    uint64_t N = 15120;
    uint ell = getNumNodes(N);
    constBiFuncNode root[ell];
    init_fft_tree(root, N);
    std::cout << "\t\t";
    printTree(root);
    std::cout << "\nExpected\t15120: (16, 945: (27, 35: (5, 7)))\n";
}


/*
 * Functions to clock my running times. These check in different cases
 * to see if what I'm doing is actually sufficiently fast
 */

// Deprecated function evaluations. Will fix in future.
void time_fft() {
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;
    typedef std::chrono::high_resolution_clock clock;

    uint64_t n = 3*3*1594323;
    auto in = (Complex*) malloc(n*sizeof(Complex));
    auto out = (Complex*) malloc(n*sizeof(Complex));
    for(uint64_t i = 0; i < n; i++) in[i] = Complex(i, 2*i);
    auto start = clock::now();
    // pow3_FFT(n, in, out, 1);
    auto end = clock::now();
    auto myfft = duration_cast<nanoseconds>(end-start).count();
    std::cout << "Elapsed time is " << myfft*1.e-9 << "s\n";
    free(in); free(out);
}

void time_complex_mult() {
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;
    typedef std::chrono::high_resolution_clock clock;

    std::vector<std::complex<double>> stdcomp1 {};
    std::vector<std::complex<double>> stdcomp2 {};
    std::vector<Complex> mycomp1 {};
    std::vector<Complex> mycomp2 {};
    auto rand = std::bind(std::uniform_real_distribution<>{0.0,10.0}, std::default_random_engine{});
    uint len = 5e6;
    double a,b;
    for(uint i = 0; i < len; i++) {
        a = rand(); b = rand();
        stdcomp1.push_back(std::complex<double>(a, b));
        mycomp1.push_back(Complex(a, b));
        a = rand(); b = rand();
        stdcomp2.push_back(std::complex<double>(a, b));
        mycomp2.push_back(Complex(a, b));
    }

    auto start = clock::now();
    for(uint i = 0; i < len; i++) {
        auto m = stdcomp1[i] * stdcomp2[i];
        if((void*) (&m) == (void*) 0x123456) std::cout << "aaaa\n";
    }
    auto end = clock::now();
    auto stdmult = duration_cast<nanoseconds>(end-start).count();
    std::cout << "Standard mult took " << stdmult << "ns\n";
    
    start = clock::now();
    for(uint i = 0; i < len; i++) {
        auto m = mycomp1[i] * mycomp2[i];
        if((void*) (&m) == (void*) 0x123456) std::cout << "aaaa\n";
    }
    end = clock::now();
    auto mymult = duration_cast<nanoseconds>(end-start).count();
    std::cout << "My mult took " << mymult << "ns\n";
    std::cout << "My mult was " << (((float) mymult/ (float) stdmult)) << " times the speed of std\n";
}

void time_const_tree() {
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;
    typedef std::chrono::high_resolution_clock clock;

    uint N = 1e5;
    std::cerr << "Testing compile time performance...\n";
    auto start = clock::now();
    for(uint64_t i = 0; i < N; i++) {
        uint64_t ell = i*i - 3*i + 4;
        uint64_t k = getNumNodes(ell);
        constBiNode<uint64_t> root[k];
        initUintConstBiTree(root, ell);
        std::cout << "ell = " << ell << ", ";
        printRoot(root, k);
    }
    auto end = clock::now();
    std::cerr << "Known at compile time took " << duration_cast<nanoseconds>(end-start).count() << "ns\n";
    std::cerr << "Initializing random ints...\n";
    auto rand = std::bind(std::uniform_int_distribution<>{2, 4}, std::default_random_engine{});
    std::vector<uint> v {};
    for(uint64_t i = 0; i < N; i++) {
        v.push_back(rand());
    }
    std::cerr << "Random ints initialized. Testing runtime performance...\n";
    start = clock::now();
    for(uint64_t i = 0; i < N; i++) {
        uint64_t ell = i*i - v[i]*i + 5;
        std::cout << "ell = " << ell << ", ";
        uint64_t k = getNumNodes(ell);
        constBiNode<uint64_t> root[k];
        initUintConstBiTree(root, ell);
        printRoot(root, k);
    }
    end = clock::now();
    std::cerr << "Known at run time took " << duration_cast<nanoseconds>(end-start).count() << "ns\n";
}

auto omega_fcn_time(std::vector<uint64_t>& Nvec) {
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;
    typedef std::chrono::high_resolution_clock clock;

    auto start = clock::now();
    for(auto& Nj : Nvec) {
        for(uint64_t i = 0; i < Nj; i++) {
            for(uint64_t j = 0; j < Nj; j++) {
                Complex w0 = omega(i*j, Nj);
                if(Nj == 16) std::cout << w0 << "\n";
            }
        }
    }
    auto end = clock::now();
    return duration_cast<nanoseconds>(end-start).count();
}

auto omega_class_time(std::vector<uint64_t>& Nvec, Omega w) {
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;
    typedef std::chrono::high_resolution_clock clock;

    auto start = clock::now();
    for(auto& Nj : Nvec) {
        for(uint64_t i = 0; i < Nj; i++) {
            for(uint64_t j = 0; j < Nj; j++) {
                Complex w0 = w(i, j, Nj);
                if(Nj == 16) std::cout << w0 << "\n";
            }
        }
    }
    auto end = clock::now();
    return duration_cast<nanoseconds>(end-start).count();
}

void time_omega() {
    uint64_t N_factors = 10;
    int Njmax = 10;
    
    auto rand = std::bind(std::uniform_int_distribution<>{1,Njmax}, std::default_random_engine{});

    uint Ntrials = 100;
    
    std::vector<uint64_t> Nvec {};
    uint64_t N = 1;
    for(uint i = 0; i < N_factors; i++) {
        uint tmp = rand();
        N *= tmp;
        Nvec.push_back(tmp);
    }
    Omega w (N, Direction::forward);
    auto fcn_time = omega_fcn_time(Nvec);
    auto class_time = omega_class_time(Nvec, w);

    for(uint i = 0; i < Ntrials; i++) {
        Nvec.clear();
        N = 1;
        for(uint i = 0; i < N_factors; i++) {
            uint tmp = rand();
            N *= tmp;
            Nvec.push_back(tmp);
        }
        w(N);

        fcn_time += omega_fcn_time(Nvec);
        class_time += omega_class_time(Nvec, w);
    }

    std::cout << "Naive implementation elapsed time is " << (fcn_time/((double) Ntrials))*1.e-9 << "s\n";

    
    
    std::cout << "Class-based implementation of seconds is " << (class_time/((double) Ntrials))*1.e-9 << "s\n";
}