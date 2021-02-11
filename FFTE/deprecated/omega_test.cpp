/**
 * Code Author: Danny Sharp
 * This file is part of FFTE (Fast Fourier Transform Engine)
 */

#include "omega_test.hpp"
using namespace FFTE;

void check_omega() {
    std::cout << "Checking the correctness of the Omega class...\n";
    size_t N1 = 9;
    size_t N2 = 5;
    size_t N3 = 7;
    size_t N = N1*N2*N3;
    Omega w (N, FFTE::Direction::forward);
    double sum = 0.;
    for(size_t i = 0; i < N; i++) {
        for(size_t j = 0; j < N; j++) sum += (w(i*j, N)-omega(i*j, N, Direction::forward)).modulus();
    }
    std::cout << "L2 Error in initialization with N = " << N << " is " << sum << "\n";
    
    sum = 0.;
    for(size_t i = 0; i < N1; i++) {
        for(size_t j = 0; j < N1; j++) sum += (w(i*j, N1)-omega(i*j, N1, Direction::forward)).modulus();
    }
    std::cout << "L2 Error in Checking against N1 = " << N1 << " is " << sum << "\n";
    
    sum = 0.;
    for(size_t i = 0; i < N2; i++) {
        for(size_t j = 0; j < N2; j++) sum += (w(i*j, N2)-omega(i*j, N2, Direction::forward)).modulus();
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
