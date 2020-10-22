#include "test.hpp"
#include "constTree.hpp"

#define CHECKSUM 1
#define ENTRYWISE 0

void check_fft() {
    int n = 1<<10;
    auto in = (Complex*) malloc(n*sizeof(Complex));
    auto out_rec = (Complex*) malloc(n*sizeof(Complex));
    auto out_comp = (Complex*) malloc(n*sizeof(Complex));
    auto out_dft = (Complex*) malloc(n*sizeof(Complex));
    for(int i = 0; i < n; i++) in[i] = Complex(i, 2*i);

    // pow3_FFT(n, in, out_rec, 1);
    pow2_FFT(in, 1, n, out_rec);
    DFT(n, in, out_dft, 1, 1);
#if CHECKSUM
    double sum = (out_dft[0] - out_rec[0]).modulus();
    for(int i = 1; i < n; i++) sum += (out_rec[i]-out_dft[i]).modulus();
    std::cout << "Norm of Error: " << sum << "\n";
#endif
#if ENTRYWISE
    for(int i = 0; i < n; i++) std::cout << "i = " << i << ", Err = " << (out_dft[i] - out_rec[i]).modulus() << "\n";
#endif
    free(in); free(out_comp); free(out_rec); free(out_dft);
}

void time_fft() {
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;
    typedef std::chrono::high_resolution_clock clock;

    uint64_t n = 3*3*1594323;
    auto in = (Complex*) malloc(n*sizeof(Complex));
    auto out = (Complex*) malloc(n*sizeof(Complex));
    for(uint64_t i = 0; i < n; i++) in[i] = Complex(i, 2*i);
    auto start = clock::now();
    pow3_FFT(n, in, out, 1);
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