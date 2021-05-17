#include "/home/dsharp/Documents/FFT-Research-F20/FFTE/include/Engine.hpp"
#include <array>
#include <cstring>
#include <chrono>
#include <math.h>
#include <fstream>

// Domain size, extended domain gridpoints
const double L  = 1.;
const size_t Ne = 1<<4;
const size_t m = Ne/2 - 1;
const double eps = 1e-14;

// Step size
double h = L / ((double) Ne/2);
auto h2 = h*h;

double* poisson();
void poissonPeriodic(std::string);

int main(int argc, char** argv) {
    poissonPeriodic("/mnt/c/Users/danny/Downloads/poissonFFTE.csv");
}


double* poisson() {
    // Allocate and initialize x, y, and f
    auto x = (double*) malloc(m*sizeof(double));
    auto y = (double*) malloc(m*sizeof(double));
    auto f = (double*) malloc(m*m*sizeof(double*));

    #pragma omp for simd
    for(size_t j = 0; j < m; j++) {
        x[j] = y[j] = h*((double) (j+1));
        for(size_t i = 0; i < m; i++) f[j*m + i] = 4.;
    }

    // Correct f on the boundaries
    #pragma omp for simd
    for(size_t j = 0; j < m; j++) {
        auto tmpy = y[j]*y[j];
        auto tmpx = x[j]*x[j];
        f[j*m] -= tmpx/h2;
        f[j*m + m-1] -= (1. + tmpx)/h2;
        f[j] -= tmpy/h2;
        f[(m-1)*m + j] -= (1. + tmpy)/h2;
    }

    // Allocate g
    auto g    = FFTE::complex_alloc<double,2>::alloc(Ne*Ne);
    
    // Zero out g
    std::memset(g, 0, Ne*Ne*sizeof(FFTE::Complex<double,2>));

    // Create g as odd extension of f in m^2 time (instead of Ne^2 time)
    for(size_t row = 0; row < m; row++) {
        auto f_row = &f[row*m];
        for(size_t col = 0; col < m; col++) {
            auto tmp = f_row[col];
            auto g_row1 = row+1; auto g_row2 = 2*m+1 - row;
            auto g_col1 = col+1; auto g_col2 = 2*m+1 - col;
            g[g_row1*Ne + g_col1] = g[g_row2*Ne + g_col2] = FFTE::Complex<double,2>(tmp, 0.);
            tmp *= -1;
            g[g_row1*Ne + g_col2] = g[g_row2*Ne + g_col1] = FFTE::Complex<double,2>(tmp, 0.);
        }
    }
    free(f);

    // Put g in the spectral domain
    auto ghat = FFTE::dim2engine<Ne,Ne>::fft2(g, FFTE::Direction::forward, FFTE::Major::row);

    // Done with g
    free(g);

    // Scale ghat according to the properties of the PDE
    for(size_t row = 0; row < Ne; row++) {
        auto s_row = sin(row*M_PI/((double) Ne));
        auto mu = 8*(s_row*s_row)/h2;
        if(row != 0) ghat[row*Ne + row] /= (-mu);

        for(size_t col = row+1; col < Ne; col++) {
            auto s_col = sin(col*M_PI/((double) Ne));
            auto mu = 4*(s_row*s_row + s_col*s_col)/h2;
            ghat[row*Ne + col] /= (-mu);
            ghat[col*Ne + row] /= (-mu);
        }
    }

    // Get the rudimentary solution
    auto v = FFTE::dim2engine<Ne,Ne>::fft2(ghat, FFTE::Direction::inverse, FFTE::Major::row);
    // Done with ghat
    free(ghat);

    // Extract real part and scale accordingly
    auto u = (double*) malloc((m+2)*(m+2)*sizeof(double));
    for(size_t row = 0; row < (m+2); row++) {
        auto v_row = &v[row*Ne];
        auto u_row = &u[row*(m+2)];
        for(size_t col = 0; col < (m+2); col++) {
            u_row[col] = v_row[col].get()[0]/(Ne*Ne);
        }
    }
    free(v);

    auto err = 0.;
    for(size_t row = 0; row < m; row++) {
        for(size_t col = 0; col < m; col++) {
            auto tmp = u[(row+1)*(m+2) + (col+1)] - (x[row]*x[row] + y[col]*y[col]);
            err += std::max<double>(tmp, -tmp);
        }
    }
    std::cout << "L2 Error: " << err << "\n";

    // Free memory allocated
    free(x); free(y);
    return u;
}

// Solve a Poisson BVP with periodic conditions
void poissonPeriodic(std::string filename) {
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;
    typedef std::chrono::high_resolution_clock clock;

    auto start = clock::now();
    const int N = 1<< _P_IDX;
    const double L = 1.;
    const double sig = 0.1;
    
    double h = L/((double) N); // Step size
    
    // Boundary conditions
    auto f_fun = [L,sig](double x, double y){
        double rsq = (x-0.5*L)*(x-0.5*L) + (y-0.5*L)*(y-0.5*L);
        double sigsq = sig*sig;
        return exp(-rsq/(2*sigsq))*(rsq - 2*sigsq)/(sigsq*sigsq);
    };
    
    // Instantiate and initialize A
    auto A = FFTE::complex_alloc<double,2>::alloc(N*N);
    for(int row = 0; row < N; row++) {
        auto A_row = &A[N*row];
        auto x = ((double) row)*h;
        for(int col = 0; col < N; col++) {
            auto y = ((double) col)*h;
            A_row[col] = FFTE::Complex<double,2> {f_fun(x,y),0.};
        }
    }
    
    // Get the transform of A
    auto A_hat = FFTE::dim2engine<N,N>::fft2(A,
                FFTE::Direction::forward, FFTE::Major::row);
    
    // Scale A_hat
    for(int row = 0; row < N; row++) {
        auto A_hat_row = &A_hat[N*row];
        auto xi_1 = (row > N/2 - 1) ? row - N : row;
        auto xi_1_sq = ((double) xi_1)*((double) xi_1);
        for(int col = 0; col < N; col++) {
            auto xi_2 = (col > N/2 - 1) ? col - N : col;
            auto xi_2_sq = ((double) xi_2)*((double) xi_2);
            // Delta Matrix entry (j,k)
            auto D_jk = -L*L/(4.*M_PI*M_PI*(xi_1_sq + xi_2_sq));
            if(row == 0 && col == 0) D_jk = 1.;
            // Scale A_hat
            A_hat_row[col] *= D_jk;
        }
    }
    
    // Get inverse transform
    auto U = FFTE::dim2engine<N,N>::fft2(A_hat,
            FFTE::Direction::inverse, FFTE::Major::row);
    
    auto u = (double*) malloc(N*N*sizeof(double));
    double disp = U[0].get()[0] / ((double) N*N);
    for(int row = 0; row < N; row++) {
        auto u_row = &u[row*N];
        auto U_row = &U[row*N];
        for(int col = 0; col < N; col++) {
            // Get real part
            u_row[col] = U_row[col].get()[0]/((double) N*N) - disp;
        }
    }
    free(A); free(A_hat); free(U);
    auto end = clock::now();
    std::ofstream myfile;
    myfile.open(filename);
    myfile << "u\n";
    for(int i = 0; i < N*N; i++) myfile << u[i] << "\n";
    myfile.close();
    free(u);
    std::cout << "p = " << _P_IDX << ", time = " << duration_cast<nanoseconds>(end-start).count()*(1.e-9) << "\n";
}