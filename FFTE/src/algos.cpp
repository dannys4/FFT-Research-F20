#include "algos.hpp"
/**
 * Code Author: Danny Sharp
 * This file is part of FFTE (Fast Fourier Transform Engine)
 */

// The following two parts here are sourced from the following link:
// https://graphics.stanford.edu/~seander/bithacks.html#BitReverseObvious
static constexpr uint8_t BitReverseTable256[256] = 
{
#   define R2(n)    n,     n + 2*64,     n + 1*64,     n + 3*64
#   define R4(n) R2(n), R2(n + 2*16), R2(n + 1*16), R2(n + 3*16)
#   define R6(n) R4(n), R4(n + 2*4 ), R4(n + 1*4 ), R4(n + 3*4 )
    R6(0), R6(2), R6(1), R6(3)
};

/*
constexpr size_t reverse(size_t v, uint8_t nbits) {
    size_t c = 0;
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
*/
namespace FFTE {
    // Calculate an omega factor for the FFT. This isn't used by me in practice.
    Complex omega(uint power, uint N, Direction dir) {
        double a = 2.*M_PI*((double) power)/((double) N);
        Complex ret {cos(a), static_cast<int>(dir)*sin(a)};
        return ret;
    }

    // Recursive helper function implementing a classic C-T FFT
    void pow2_FFT_helper(size_t N, Complex* x, Complex* y, size_t s_in, size_t s_out, Omega& w) {
        
        // Trivial case
        if(N == 1) {
            *y = *x;
            return;
        }

        // Size of sub-problem
        int m = N/2;

        // Divide into two sub-problems
        pow2_FFT_helper(m, x, y, s_in*2, s_out, w);
        pow2_FFT_helper(m, x+s_in, y+s_out*m, s_in*2, s_out, w);

        // Conquer larger problem accordingly
        for(int j = 0; j < m; j++) {
            Complex wj = w(j, N);
            int j_stride = j*s_out;
            int jm_stride = (j+m)*s_out;
            Complex y_j = y[j_stride];
            y[j_stride] = y_j + wj*y[jm_stride];
            y[jm_stride] = y_j - wj*y[jm_stride];
        }
    }

    // External function to call the C-T radix-2 FFT
    void pow2_FFT(Complex* x, Complex* y, size_t s_in, size_t s_out, biFuncNode* sRoot, Omega& w) {
        const size_t N = sRoot->sz; // Size of problem
        if(!w()) w(N);                // Initialize the Omega if necessary
        pow2_FFT_helper(N, x, y, s_in, s_out, w); // Call the radix-2 FFT
    }

    // Internal helper function to perform a DFT in O(n^2) time
    void DFT_helper(size_t size, Complex* sig_in, Complex* sig_out, size_t s_in, size_t s_out, Omega& w) {
        
        // Find the first element 
        Complex tmp = sig_in[0];
        for(uint n = 1; n < size; n++) {
            tmp = tmp + sig_in[n*s_in];
        }
        sig_out[0] = tmp;

        // Find the rest of the elements
        for(uint k = 1; k < size; k++) {
            tmp = sig_in[0];
            for(uint n = 1; n < size; n++) {
                tmp = tmp + w(k, n, size)*sig_in[n*s_in];
            }
            sig_out[k*s_out] = tmp;
        }
    }

    // Internal helper function to perform a DFT without an Omega object (there's no reason to initialize one)
    void DFT_helper(size_t size, Complex* sig_in, Complex* sig_out, size_t s_in, size_t s_out, Direction dir) {
        // Twiddle with smallest numerator
        Complex w0 = omega(1, size, dir);

        // Base twiddle for each outer iteration
        Complex wk = w0;

        // Twiddle for inner iterations
        Complex wkn = w0;
        
        // Calculate first element of output
        Complex tmp = sig_in[0];
        for(uint n = 1; n < size; n++) {
            tmp = tmp + sig_in[n*s_in];
        }
        sig_out[0] = tmp;

        // Initialize rest of output
        for(uint k = 1; k < size; k++) {
            // Initialize kth output
            tmp = sig_in[0];

            // Calculate kth output
            for(uint n = 1; n < size; n++) {
                tmp = tmp + wkn*sig_in[n*s_in];
                wkn = wkn*wk;
            }
            sig_out[k*s_out] = tmp;

            // "Increment" wk and "reset" wkn
            wk = wk*w0;
            wkn = wk;
        }
    }

    // External-facing function to properly call the internal DFT function
    void DFT(Complex* x, Complex* y, size_t s_in, size_t s_out, biFuncNode* sLeaf, Omega& w) {
        if(!w()) {
            DFT_helper(sLeaf->sz, x, y, s_in, s_out, w.dir);
        }
        else DFT_helper(sLeaf->sz, x, y, s_in, s_out, w);
    }

    // External-facing reference DFT for testing purposes
    void reference_DFT(size_t N, Complex* x, Complex* y, Direction dir) {
        DFT_helper(N, x, y, 1, 1, dir);
    }

    // External & Internal function for radix-N1 C-T FFTs
    void composite_FFT(Complex* x, Complex* y, size_t s_in, size_t s_out, biFuncNode* sRoot, Omega& w) {
        // Retrieve N
        size_t N = sRoot->sz;
        
        // Initialize the omega if necessary
        if(!w()) w(N);

        // Find the children on the call-graph
        biFuncNode* left = sRoot + sRoot->left;
        biFuncNode* right = sRoot + sRoot->right;

        // Get the size of the sub-problems
        uint N1 = left->sz;
        uint N2 = right->sz;

        /* Theoretically, this should've ever be called-- the call graph should've
        * told us to perform a DFT here instead. However, I'm calling this just in case.
        */
        if(N1 == N) {
            DFT(x, y, s_in, s_out, sRoot, w);
            return;
        }

        // I'm currently using a temporary storage space malloc'd in recursive calls.
        // This isn't optimal and will change as the engine develops
        Complex* z = (Complex*) malloc(N*sizeof(Complex));

        // Find the FFT of the "rows" of the input signal and twiddle them accordingly
        for(uint n1 = 0; n1 < N1; n1++) {
            right->fptr(x+n1*s_in, z+N2*n1, N1*s_in, 1, right, w);
            for(uint k2 = 1; (k2 < N2) && (n1 > 0); k2++) {
                z[n1*N2 + k2] = z[n1*N2 + k2]*w(n1, k2, N);
            }
        }

        /* Take the FFT of the "columns" of the output from the above "row" FFTs. 
        * Don't need n1 for the second transform as it's just the indexer.
        * Take strides of N2 since z is allocated on the fly in this function for N.
        */
        for(uint k2 = 0; k2 < N2; k2++) {
            left->fptr(z+k2, y+k2*s_out, N2, N2*s_out, left, w);
        }

        // Free z when possible
        free(z);
    }

    // A factoring function for the reference composite FFT
    size_t referenceFactor(const size_t f) {
        // Check if it's even
        if((f & 0x1) == 0) return 2;

        // Check all odd numbers after that
        for(int k = 3; k*k < f; k+=2) {
            if( f % k == 0 ) return k;
        }

        // return f if no factor was found
        return f;
    }

    // Internal-facing radix-N1 C-T FFT for reference. Doesn't use any call graph or Omega class.
    // Used for checking correctness of output
    void reference_composite_FFT_helper(size_t N, Complex* x, Complex* y, size_t s_in, size_t s_out, Direction dir) {
        
        // Find the factors of N
        uint N1 = referenceFactor(N);
        uint N2 = N/N1;
        
        // Execute a DFT if necessary
        if(N1 == N) {
            DFT_helper(N, x, y, s_in, s_out, dir);
            return;
        }

        // Malloc a temporary array to hold intermediate results
        Complex* z = (Complex*) malloc(N*sizeof(Complex));

        // Take the FFT of the "rows" and put output into z.
        reference_composite_FFT_helper(N2, x, z, N1*s_in, 1, dir);
        for(uint n1 = 1; n1 < N1; n1++) {
            reference_composite_FFT_helper(N2, x+n1*s_in, z+N2*n1, N1*s_in, 1, dir);
            for(uint k2 = 1; k2 < N2; k2++) {
                z[n1*N2 + k2] = z[n1*N2 + k2]*omega(n1*k2, N, dir);
            }
        } 

        // Take the FFT of the "columns" of z and put it into y
        // Don't need n1 for the second transform as it's just the indexer.
        for(uint k2 = 0; k2 < N2; k2++) {
            reference_composite_FFT_helper(N1, z+k2, y+k2*s_out, N2, N2*s_out, dir);
        }

        // Free memory used by z
        free(z);
    }

    // External-facing function to calculate a composite FFT using the reference code.
    void reference_composite_FFT(size_t N, Complex* x, Complex* y, Direction dir) {
        reference_composite_FFT_helper(N, x, y, 1, 1, dir);
    }

    // Internal recursive helper-function that calculates the FFT of a signal with length 3^k
    void pow3_FFT_helper(size_t N, Complex* x, Complex* y, size_t s_in, size_t s_out, Omega& w, Complex& plus120, Complex& minus120) {
        
        // Calculate the DFT manually if necessary
        if(N == 3) {
            y[0] = x[0] + x[s_in] + x[2*s_in];
            y[s_out] = x[0] + plus120*x[s_in] + minus120*x[2*s_in];
            y[2*s_out] = x[0] + minus120*x[s_in] + plus120*x[2*s_in];
            return;
        }

        // Calculate the size of the sub-problem
        size_t Nprime = N/3;

        // Divide into sub-problems
        pow3_FFT_helper(Nprime, x, y, s_in*3, s_out, w, plus120, minus120);
        pow3_FFT_helper(Nprime, x+s_in, y+Nprime*s_out, s_in*3, s_out, w, plus120, minus120);
        pow3_FFT_helper(Nprime, x+2*s_in, y+2*Nprime*s_out, s_in*3, s_out, w, plus120, minus120);

        // Combine the sub-problem solutions
        for(size_t k = 0; k < Nprime; k++) {
            // Index calculation
            auto k1 = k * s_out;
            auto k2 = (  Nprime + k) * s_out;
            auto k3 = (2*Nprime + k) * s_out;

            // Storing temporary variables
            Complex tmpk     = y[k1];
            Complex tmpk_p_1 = y[k2];
            Complex tmpk_p_2 = y[k3];

            // Twiddle factors
            Complex wk1 = w(k, N); Complex wk2 = w(2*k, N);

            // Reassigning the output
            y[k1] = tmpk +             wk1 * tmpk_p_1 +            wk2 * tmpk_p_2;
            y[k2] = tmpk + plus120  *  wk1 * tmpk_p_1 + minus120 * wk2 * tmpk_p_2;
            y[k3] = tmpk + minus120 *  wk1 * tmpk_p_1 + plus120  * wk2 * tmpk_p_2;
        }
    }

    // External-facing function for performing an FFT on signal with length N = 3^k
    void pow3_FFT(Complex* x, Complex* y, size_t s_in, size_t s_out, biFuncNode* sRoot, Omega& w) {
        const size_t N = sRoot->sz;
        if(!w()) w(N);
        Complex plus120 {-0.5, -sqrt(3)/2.};
        Complex minus120 {-0.5, sqrt(3)/2.};
        switch(w.dir) {
            case Direction::forward: pow3_FFT_helper(N, x, y, s_in, s_out, w, plus120, minus120); break;
            case Direction::inverse: pow3_FFT_helper(N, x, y, s_in, s_out, w, minus120, plus120); break;
        }
    }

#if 0
    // DOES NOT WORK, DO NOT USE
    void pow2_FFT_it(Complex* sig_in, size_t size, Complex* sig_out) {
        int log_2_sz = 0;
        while(size >> (log_2_sz+1)) log_2_sz++;

        Complex w {cos(2*M_PI/(double) size), -sin(2*M_PI/(double) size)};
        Complex* Omega = (Complex*) malloc(size*sizeof(Complex)>>1);
        *Omega = Complex(1., 0.);
        sig_out[0] = sig_in[0];
        sig_out[size-1] = sig_in[size-1];
        for(size_t i = 1; i < (size>>1); i++){
            auto rev_i = reverse(i, log_2_sz);
            sig_out[i] = sig_in[rev_i];
            sig_out[rev_i] = sig_in[i];
            Omega[i] = Omega[i-1]*w;
        }

        for(int i = 1; i <= log_2_sz; i++) {
            size_t L_o = 2  << (i-1);
            size_t r   = size >> i;
            for(size_t j = 0; j < r; j++) {
                for(size_t k = 0; k < L_o; k++) {
                    Complex sig_out_tmp = sig_out[j*L_o+k];
                    sig_out[j*L_o+k]     = sig_out_tmp + Omega[k] * sig_out[(j+1)*L_o+k];
                    sig_out[(j+1)*L_o+k] = sig_out_tmp - Omega[k] * sig_out[(j+1)*L_o+k];
                }
            }
        }
    }
#endif
}