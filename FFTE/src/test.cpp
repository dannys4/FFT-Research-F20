/**
 * Code Author: Danny Sharp
 * This file is part of FFTE (Fast Fourier Transform Engine)
 */

#include "test.hpp"
#define CHECKSUM_COMP 1
#define ENTRYWISE_COMP 0
#define FFT_LENGTH 32

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

void test_FFT_into_csv(std::string filename, int maxsize) {
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;
    typedef std::chrono::high_resolution_clock clock;
    auto rand = std::bind(std::uniform_real_distribution<>{0.0,10.0}, std::default_random_engine{});

    std::cout << "Timing the current FFT tree configuration against reference...\n";
    std::cout << "Outputting results into " << filename << "\n";
    std::ofstream myfile;
    Complex<double, 2> *x = (Complex<double, 2>*) malloc(maxsize*sizeof(Complex<double, 2>));
    Complex<double, 2> *y_new = (Complex<double, 2>*) malloc(maxsize*sizeof(Complex<double, 2>));
    Complex<double, 2> *y_ctc = (Complex<double, 2>*) malloc(maxsize*sizeof(Complex<double, 2>));
    Complex<double, 2> *y_dft = (Complex<double, 2>*) malloc(maxsize*sizeof(Complex<double, 2>));
    for(int i = 0; i < maxsize; i++) x[i] = Complex<double, 2>(rand(),rand());

    myfile.open(filename);
    myfile << "N,FFTE Setup Time(s),FFTE Time (s),Composite FFT Time (s),DFT Time (s),L2 Error\n";
    for(int N = 10; N < maxsize; N++) {
        auto start = clock::now();
        int ell = getNumNodes(N);
        biFuncNode<double, 2> root[ell];
        init_fft_tree(root, N);
        Direction dir = Direction::forward;
        auto end = clock::now();
        auto setupTime = duration_cast<nanoseconds>(end-start).count();

        start = clock::now();
        root->fptr(x, y_new, 1, 1, root, dir);
        end = clock::now();
        auto myTime = duration_cast<nanoseconds>(end-start).count();

        start = clock::now();
        reference_composite_FFT(N, x, y_ctc, dir);
        end = clock::now();
        auto compositeTime = duration_cast<nanoseconds>(end-start).count();

        start = clock::now();
        reference_DFT(N, x, y_dft, dir);
        end = clock::now();
        auto dftTime = duration_cast<nanoseconds>(end-start).count();

        double err = 0.;
        for(int j = 0; j < N; j++) {
            err += (y_new[j]-y_dft[j]).modulus()[0];
        }

        myfile << N << "," << setupTime << "," << myTime << "," << compositeTime << "," << dftTime << "," << err << "\n";
    }
    myfile.close();
    std::cout << "Finished!\n";
}

// Compares my complex multiplication to std::complex multiplication
void time_complex_mult() {
    std::cout << "Comparing the performance of my Complex<double, 2> multiplication vs. std::complex...\n";
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;
    typedef std::chrono::high_resolution_clock clock;

    std::vector<std::complex<double>> stdcomp1 {};
    std::vector<std::complex<double>> stdcomp2 {};
    std::vector<Complex<double, 2>> mycomp1 {};
    std::vector<Complex<double, 2>> mycomp2 {};
    auto rand = std::bind(std::uniform_real_distribution<>{0.0,10.0}, std::default_random_engine{});
    size_t len = 5e6;
    double a,b;
    for(size_t i = 0; i < len; i++) {
        a = rand(); b = rand();
        stdcomp1.push_back(std::complex<double>(a, b));
        mycomp1.push_back(Complex<double, 2>(a, b));
        a = rand(); b = rand();
        stdcomp2.push_back(std::complex<double>(a, b));
        mycomp2.push_back(Complex<double, 2>(a, b));
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
    std::vector<unsigned long long> v {};
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