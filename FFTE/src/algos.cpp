#include "algos.hpp"

static constexpr uint8_t BitReverseTable256[256] = 
{
#   define R2(n)    n,     n + 2*64,     n + 1*64,     n + 3*64
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

Complex omega(uint power, uint N) {
    double a = 2.*M_PI*((double) power)/((double) N);
    Complex ret {cos(a), -sin(a)};
    return ret;
}

void pow2_FFT_helper(uint64_t N, Complex* x, Complex* y, uint64_t s_in, uint64_t s_out, Omega& w) {
    if(N == 1) {
        *y = *x;
        return;
    }
    int m = N/2;
    pow2_FFT_helper(m, x, y, s_in*2, s_out, w);
    pow2_FFT_helper(m, x+s_in, y+s_out*m, s_in*2, s_out, w);
    for(int j = 0; j < m; j++) {
        Complex wj = w(j, N);
        int j_stride = j*s_out;
        int jm_stride = (j+m)*s_out;
        Complex y_j = y[j_stride];
        y[j_stride] = y_j + wj*y[jm_stride];
        y[jm_stride] = y_j - wj*y[jm_stride];
    }
}

void pow2_FFT(Complex* x, Complex* y, uint64_t s_in, uint64_t s_out, constBiFuncNode* sRoot, Omega& w) {
    const uint64_t N = sRoot->sz;
    if(!w()) w(N);
    pow2_FFT_helper(N, x, y, s_in, s_out, w);
}

void DFT_helper(uint64_t size, Complex* sig_in, Complex* sig_out, uint64_t s_in, uint64_t s_out, Omega& w) {
    Complex tmp = sig_in[0];
    for(uint n = 1; n < size; n++) {
        tmp = tmp + sig_in[n*s_in];
    }
    sig_out[0] = tmp;

    for(uint k = 1; k < size; k++) {
        tmp = sig_in[0];
        for(uint n = 1; n < size; n++) {
            tmp = tmp + w(k, n, size)*sig_in[n*s_in];
        }
        sig_out[k*s_out] = tmp;
    }
}

void DFT_helper(uint64_t size, Complex* sig_in, Complex* sig_out, uint64_t s_in, uint64_t s_out) {
    Complex w0 = omega(1, size);
    Complex wk = w0;
    Complex wkn = w0;
    Complex tmp = sig_in[0];
    for(uint n = 1; n < size; n++) {
        tmp = tmp + sig_in[n*s_in];
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

void DFT(Complex* x, Complex* y, uint64_t s_in, uint64_t s_out, constBiFuncNode* sLeaf, Omega& w) {
    if(!w()) DFT_helper(sLeaf->sz, x, y, s_in, s_out);
    else DFT_helper(sLeaf->sz, x, y, s_in, s_out, w);
}

void reference_DFT(uint64_t N, Complex* x, Complex* y) {
    DFT_helper(N, x, y, 1, 1);
}

void composite_FFT(Complex* x, Complex* y, uint64_t s_in, uint64_t s_out, constBiFuncNode* sRoot, Omega& w) {
    uint64_t N = sRoot->sz;
    if(!w()) w(N);
    constBiFuncNode* left = sRoot + sRoot->left;
    constBiFuncNode* right = sRoot + sRoot->right;
    uint N1 = left->sz;
    uint N2 = right->sz;

    if(N1 == N) {
        DFT(x, y, s_in, s_out, sRoot, w);
        return;
    }
    Complex* z = (Complex*) malloc(N*sizeof(Complex));

    for(uint n1 = 0; n1 < N1; n1++) {
        right->fptr(x+n1*s_in, z+N2*n1, N1*s_in, 1, right, w);
        for(uint k2 = 1; (k2 < N2) && (n1 > 0); k2++) {
            z[n1*N2 + k2] = z[n1*N2 + k2]*w(n1, k2, N);
        }
    }

    // Don't need n1 for the second transform as it's just the indexer.
    // Take strides of N2 since z is allocated on the fly in this function for N.
    for(uint k2 = 0; k2 < N2; k2++) {
        left->fptr(z+k2, y+k2*s_out, N2, N2*s_out, left, w);
    }
    free(z);
}

void reference_composite_FFT(uint64_t N, Complex* x, Complex* y, uint64_t s_in, uint64_t s_out) {
    uint N1 = factor(N);
    uint N2 = N/N1;
    if(N1 == N) {
        DFT_helper(N, x, y, s_in, s_out);
        return;
    }
    Complex* z = (Complex*) malloc(N*sizeof(Complex));

    reference_composite_FFT(N2, x, z, N1*s_in, 1);
    for(uint n1 = 1; n1 < N1; n1++) {
        reference_composite_FFT(N2, x+n1*s_in, z+N2*n1, N1*s_in, 1);
        for(uint k2 = 1; k2 < N2; k2++) {
            z[n1*N2 + k2] = z[n1*N2 + k2]*omega(n1*k2, N);
        }
    }
    // Don't need n1 for the second transform as it's just the indexer. Take strides of N2
    for(uint k2 = 0; k2 < N2; k2++) {
        reference_composite_FFT(N1, z+k2, y+k2*s_out, N2, N2*s_out);
    }
    free(z);
}

void pow3_FFT_helper(uint64_t N, Complex* x, Complex* y, uint64_t s_in, uint64_t s_out, Omega& w, Complex& plus120, Complex& minus120) {
    if(N == 3) {
        y[0] = x[0] + x[s_in] + x[2*s_in];
        y[s_out] = x[0] + plus120*x[s_in] + minus120*x[2*s_in];
        y[2*s_out] = x[0] + minus120*x[s_in] + plus120*x[2*s_in];
        return;
    }
    uint64_t Nprime = N/3;
    pow3_FFT_helper(Nprime, x, y, s_in*3, s_out, w, plus120, minus120);
    pow3_FFT_helper(Nprime, x+s_in, y+Nprime*s_out, s_in*3, s_out, w, plus120, minus120);
    pow3_FFT_helper(Nprime, x+2*s_in, y+2*Nprime*s_out, s_in*3, s_out, w, plus120, minus120);
    Complex w1 {cos(2*M_PI/(double) N), -sin(2*M_PI/(double) N)};
    Complex w2 = w1 * w1;
    Complex wk1 {1, 0};
    Complex wk2 {1, 0};
    for(uint64_t k = 0; k < Nprime; k++) {
        auto k1 = k*s_out;
        auto k2 = (Nprime+k)*s_out;
        auto k3 = (2*Nprime+k)*s_out;
        Complex tmpk = y[k1];
        Complex tmpk_p_1 = y[k2];
        Complex tmpk_p_2 = y[k3];
        y[k1] = tmpk + wk1*tmpk_p_1 + wk2*tmpk_p_2;
        y[k2] = tmpk + plus120*wk1*tmpk_p_1 + minus120*wk2*tmpk_p_2;
        y[k3] = tmpk + minus120*wk1*tmpk_p_1 + plus120*wk2*tmpk_p_2;
        wk1 = wk1*w1;
        wk2 = wk2*w2;
    }
}

void pow3_FFT(Complex* x, Complex* y, uint64_t s_in, uint64_t s_out, constBiFuncNode* sRoot, Omega& w) {
    const uint64_t N = sRoot->sz;
    if(!w()) w(N);
    Complex plus120 {-0.5, -sqrt(3)/2.};
    Complex minus120 {-0.5, sqrt(3)/2.};
    switch(w.dir) {
        case Direction::forward: pow3_FFT_helper(N, x, y, s_in, s_out, w, plus120, minus120); break;
        case Direction::inverse: pow3_FFT_helper(N, x, y, s_in, s_out, w, minus120, plus120); break;
    }
}

// DOES NOT WORK, DO NOT USE
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