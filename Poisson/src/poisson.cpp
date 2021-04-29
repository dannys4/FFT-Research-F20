#include "/home/dsharp/Documents/FFT-Research-F20/FFTE/include/Engine.hpp"
#include<array>
#include<cstring>

// Domain size, extended domain gridpoints
const double L  = 1.;
const size_t Ne = 1<<4;
const size_t m = Ne/2 - 1;
const double eps = 1e-14;

// Step size
double h = L / ((double) Ne/2);
auto h2 = h*h;

double* poisson();

int main(int argc, char** argv) {
    auto u = poisson();
    for(size_t i = 0; i < (m+2); i++) {
        auto row = &u[i*(m+2)];
        for(size_t j = 0; j < (m+2); j++) {
            auto tmp = row[j];
            std::cout << ((tmp*tmp > eps) ? tmp : 0.) << ((j == (m+1)) ? ";\n" : ", ");
        }
    }
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

    auto err = 0.;
    for(size_t row = 0; row < m; row++) {
        for(size_t col = 0; col < m; col++) {
            auto tmp = u[(row+1)*(m+2) + (col+1)] - (x[row]*x[row] + y[col]*y[col]);
            err += std::max<double>(tmp, -tmp);
        }
    }
    std::cout << "L2 Error: " << err << "\n";

    // Free memory allocated
    free(x); free(y); free(v);
    return u;
}