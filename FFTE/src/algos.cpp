#include "algos.hpp"

static constexpr uint8_t BitReverseTable256[256] = 
{
#   define R2(n)     n,     n + 2*64,     n + 1*64,     n + 3*64
#   define R4(n) R2(n), R2(n + 2*16), R2(n + 1*16), R2(n + 3*16)
#   define R6(n) R4(n), R4(n + 2*4 ), R4(n + 1*4 ), R4(n + 3*4 )
    R6(0), R6(2), R6(1), R6(3)
};

constexpr uint64_t reverse(uint64_t v, uint8_t nbits) {
    uint64_t c = 0;
    uint8_t * p = (uint8_t *) &v;
    uint8_t * q = (uint8_t *) &c;
    q[7] = BitReverseTable256[p[0]]; 
    q[6] = BitReverseTable256[p[1]]; 
    q[5] = BitReverseTable256[p[2]]; 
    q[4] = BitReverseTable256[p[3]];
    q[3] = BitReverseTable256[p[4]]; 
    q[2] = BitReverseTable256[p[5]]; 
    q[1] = BitReverseTable256[p[6]]; 
    q[0] = BitReverseTable256[p[7]];
    return c >> (64-nbits);
}

#define FACTORS_LEN 15
static const uint factors[FACTORS_LEN] {4, 2, 3, 5, 7, 11, 13, 16, 17, 19, 23, 29, 31, 37, 41};

Complex omega(uint power, uint N) {
    double a = 2.*M_PI*((double) power)/((double) N);
    Complex ret {cos(a), -sin(a)};
    return ret;
}

uint factor(uint f) {
    uint k;
    for(k = 0; k < FACTORS_LEN; k++) {
        if( f % factors[k] == 0) return factors[k];
    }
    for(k = factors[k - 1]; k*k < f; k+=2) {
        if( f % k == 0 ) return k;
    }
    return f;
}

void pow2_FFT(Complex* sig_in, uint64_t stride, uint64_t size, Complex* sig_out) {
    if(size == 1) {
        *sig_out = *sig_in;
        return;
    }
    int m = size/2;
    pow2_FFT(sig_in, stride*2, m, sig_out);
    pow2_FFT(sig_in+stride, stride*2, m, sig_out+m);
    Complex w {cos(2*M_PI/(double) size), -sin(2*M_PI/(double) size)};
    Complex wj {1, 0};
    for(int j = 0; j < m; j++) {
        Complex sig_out_j = sig_out[j];
        sig_out[j] = sig_out_j + wj*sig_out[j+m];
        sig_out[j+m] = sig_out_j - wj*sig_out[j+m];
        wj = w*wj;
    }
}

void DFT(uint64_t size, Complex* sig_in, Complex* sig_out, uint64_t s_in, uint64_t s_out) {
    Complex w0 = omega(1, size);
    Complex wk = w0;
    Complex wkn = w0;
    Complex tmp = sig_in[0];
    for(uint n = 0; n < size; n++) {
        tmp = tmp + sig_in[n];
    }
    sig_out[0] = tmp;

    for(uint k = 1; k < size; k++) {
        tmp = sig_in[0];
        for(uint n = 1; n < size; n++) {
            tmp = tmp + wkn*sig_in[n*s_in];
            wkn = wkn*wk;
        }
        sig_out[k*s_out] = tmp;
        wk = wk*w0;
        wkn = wk;
    }
}

void composite_FFT(uint64_t N, Complex* x, Complex* y, uint64_t s_in, uint64_t s_out) {
    uint N1 = factor(N);
    uint N2 = N/N1;

    if(N1 == N) {
        DFT(N, x, y, s_in, s_out);
        return;
    }
    Complex* z = (Complex*) malloc(N*sizeof(Complex));

    for(uint n1 = 0; n1 < N1; n1++) {
        composite_FFT(N2, x+n1*s_in, z+N2*n1, N1*s_in, 1);
        for(uint k2 = 1; (k2 < N2) && (n1 > 0); k2++) {
            z[n1*N2 + k2] = z[n1*N2 + k2]*omega(n1*k2, N);
        }
    }

    // Don't need n1 for the second transform as it's just the indexer. Take strides of N2
    for(uint k2 = 0; k2 < N2; k2++) {
        composite_FFT(N1, z+k2, y+k2, N2, N2*s_out);
    }
    free(z);
}

void pow3_FFT(uint64_t N, Complex* x, Complex* y, uint64_t s_in) {
    Complex plus60 {-0.5, -sqrt(3)/2};
    Complex minus60 {-0.5, sqrt(3)/2};
    if(N == 3) {
        y[0] = x[0] + x[s_in] + x[2*s_in];
        y[1] = x[0] + plus60*x[s_in] + minus60*x[2*s_in];
        y[2] = x[0] + minus60*x[s_in] + plus60*x[2*s_in];
        return;
    }
    int Nprime = N/3;
    pow3_FFT(Nprime, x, y, s_in*3);
    pow3_FFT(Nprime, x+s_in, y+Nprime, s_in*3);
    pow3_FFT(Nprime, x+2*s_in, y+2*Nprime, s_in*3);
    Complex w1 {cos(2*M_PI/(double) N), -sin(2*M_PI/(double) N)};
    Complex w2 = w1 * w1;
    Complex wk1 {1, 0};
    Complex wk2 {1, 0};
    for(int k = 0; k < Nprime; k++) {
        Complex tmpk = y[k];
        Complex tmpk_p_1 = y[Nprime + k];
        Complex tmpk_p_2 = y[2*Nprime + k];
        y[k] = tmpk + wk1*tmpk_p_1 + wk2*tmpk_p_2;
        y[Nprime + k] = tmpk + plus60*wk1*tmpk_p_1 + minus60*wk2*tmpk_p_2;
        y[2*Nprime + k] = tmpk + minus60*wk1*tmpk_p_1 + plus60*wk2*tmpk_p_2;
        wk1 = wk1*w1;
        wk2 = wk2*w2;
    }
}

void pow2_FFT_it(Complex* sig_in, uint64_t size, Complex* sig_out) {
    int log_2_sz = 0;
    while(size >> (log_2_sz+1)) log_2_sz++;

    Complex w {cos(2*M_PI/(double) size), -sin(2*M_PI/(double) size)};
    Complex* Omega = (Complex*) malloc(size*sizeof(Complex)>>1);
    *Omega = Complex(1., 0.);
    sig_out[0] = sig_in[0];
    sig_out[size-1] = sig_in[size-1];
    for(uint64_t i = 1; i < (size>>1); i++){
        auto rev_i = reverse(i, log_2_sz);
        sig_out[i] = sig_in[rev_i];
        sig_out[rev_i] = sig_in[i];
        Omega[i] = Omega[i-1]*w;
    }

    for(int i = 1; i <= log_2_sz; i++) {
        uint64_t L_o = 2  << (i-1);
        uint64_t r   = size >> i;
        for(uint64_t j = 0; j < r; j++) {
            for(uint64_t k = 0; k < L_o; k++) {
                Complex sig_out_tmp = sig_out[j*L_o+k];
                sig_out[j*L_o+k]     = sig_out_tmp + Omega[k] * sig_out[(j+1)*L_o+k];
                sig_out[(j+1)*L_o+k] = sig_out_tmp - Omega[k] * sig_out[(j+1)*L_o+k];
            }
        }
    }
}