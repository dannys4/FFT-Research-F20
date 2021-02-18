/**
 * Code Author: Danny Sharp
 * This file is part of FFTE (Fast Fourier Transform Engine)
 */
#include <iostream>
namespace FFTE {

    // Recursive helper function implementing a classic C-T FFT
	template<typename F, int L>
	inline void pow2_FFT_helper(size_t N, Complex<F,L>* x, Complex<F,L>* y, size_t s_in, size_t s_out, Direction dir) {
        
        // Trivial case
        if(N == 1) {
            *y = *x;
            return;
        }

        // Size of sub-problem
        int m = N/2;

        // Divide into two sub-problems
        pow2_FFT_helper(m, x, y, s_in*2, s_out, dir);
        pow2_FFT_helper(m, x+s_in, y+s_out*m, s_in*2, s_out, dir);

        // Twiddle Factor
        double inc = 2.*M_PI/N;
        Complex<F,L> w1 (cos(inc), static_cast<int>(dir)*sin(inc));
        Complex<F,L> wj (1., 0.);

        // Conquer larger problem accordingly
        for(int j = 0; j < m; j++) {
            int j_stride = j*s_out;
            int jm_stride = (j+m)*s_out;
            Complex<F,L> y_j = y[j_stride];
            y[j_stride] = y_j + wj*y[jm_stride];
            y[jm_stride] = y_j - wj*y[jm_stride];
            wj *= w1;
        }
    }

    // External function to call the C-T radix-2 FFT
	template<typename F, int L>
	inline void pow2_FFT(Complex<F,L>* x, Complex<F,L>* y, size_t s_in, size_t s_out, biFuncNode<F,L>* sRoot, Direction dir) {
        const size_t N = sRoot->sz; // Size of problem
        pow2_FFT_helper(N, x, y, s_in, s_out, dir); // Call the radix-2 FFT
    }

    // Internal helper function to perform a DFT
	template<typename F, int L>
	inline void DFT_helper(size_t size, Complex<F,L>* sig_in, Complex<F,L>* sig_out, size_t s_in, size_t s_out, Direction dir) {
        // Twiddle with smallest numerator
        Complex<F,L> w0 = omega<F,L>::get(1, size, dir);

        // Base twiddle for each outer iteration
        Complex<F,L> wk = w0;

        // Twiddle for inner iterations
        Complex<F,L> wkn = w0;
        
        // Calculate first element of output
        Complex<F,L> tmp = sig_in[0];
        for(size_t n = 1; n < size; n++) {
            tmp = tmp + sig_in[n*s_in];
        }
        sig_out[0] = tmp;

        // Initialize rest of output
        for(size_t k = 1; k < size; k++) {
            // Initialize kth output
            tmp = sig_in[0];

            // Calculate kth output
            for(size_t n = 1; n < size; n++) {
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
	template<typename F, int L>
	inline void DFT(Complex<F,L>* x, Complex<F,L>* y, size_t s_in, size_t s_out, biFuncNode<F,L>* sLeaf, Direction dir) {
        DFT_helper(sLeaf->sz, x, y, s_in, s_out, dir);
    }

    // External-facing reference DFT for testing purposes
	template<typename F, int L>
	inline void reference_DFT(size_t N, Complex<F,L>* x, Complex<F,L>* y, Direction dir) {
        DFT_helper(N, x, y, 1, 1, dir);
    }

    // External & Internal function for radix-N1 C-T FFTs
	template<typename F, int L>
	inline void composite_FFT(Complex<F,L>* x, Complex<F,L>* y, size_t s_in, size_t s_out, biFuncNode<F,L>* sRoot, Direction dir) {
        // Retrieve N
        size_t N = sRoot->sz;

        // Find the children on the call-graph
        biFuncNode<F,L>* left = sRoot + sRoot->left;
        biFuncNode<F,L>* right = sRoot + sRoot->right;

        // Get the size of the sub-problems
        size_t N1 = left->sz;
        size_t N2 = right->sz;

        /* Theoretically, this shouldn't ever be called-- the call graph should've
        * told us to perform a DFT here instead. However, I'm calling this just in case.
        */
        if(N1 == N) {
            DFT(x, y, s_in, s_out, sRoot, dir);
            return;
        }

        // I'm currently using a temporary storage space malloc'd in recursive calls.
        // This isn't optimal and will change as the engine develops
        auto z = (Complex<F,L>*) aligned_alloc(sizeof(Complex<F,L>), N*sizeof(Complex<F,L>));

        // Find the FFT of the "rows" of the input signal and twiddle them accordingly
        Complex<F,L> w1 = omega<F,L>::get(1, N, dir);
        Complex<F,L> wn1 = Complex<F,L>(1., 0.);
        for(size_t n1 = 0; n1 < N1; n1++) {
            Complex<F,L> wk2 = wn1;
            right->fptr(x+n1*s_in, z+N2*n1, N1*s_in, 1, right, dir);
            for(size_t k2 = 1; (k2 < N2) && (n1 > 0); k2++) {
                z[n1*N2 + k2] = z[n1*N2 + k2]*wk2;
                wk2 *= wn1;
            }
            wn1 *= w1;
        }

        /* Take the FFT of the "columns" of the output from the above "row" FFTs. 
        * Don't need n1 for the second transform as it's just the indexer.
        * Take strides of N2 since z is allocated on the fly in this function for N.
        */
        for(size_t k2 = 0; k2 < N2; k2++) {
            left->fptr(z+k2, y+k2*s_out, N2, N2*s_out, left, dir);
        }

        // Free z when possible
        free(z);
    }

    // A factoring function for the reference composite FFT
    inline size_t referenceFactor(const size_t f) {
        // Check if it's even
        if((f & 0x1) == 0) return 2;

        // Check all odd numbers after that
        for(int k = 3; k*k < f; k+=2) {
            if( f % k == 0 ) return k;
        }

        // return f if no factor was found
        return f;
    }

    // Implementation for Rader's Algorithm
	template<typename F, int L>
	inline void rader_FFT(Complex<F,L>* x, Complex<F,L>* y, size_t s_in, size_t s_out, biFuncNode<F,L>* sRoot, Direction dir, size_t a, size_t ainv) {
        // Size of the problem
        size_t p = sRoot->sz;
        
        // Find the children on the call-graph
        biFuncNode<F,L>* subFFT = sRoot + sRoot->left;

        // Temporary workspace
        auto z = (Complex<F,L>*) aligned_alloc(sizeof(Complex<F,L>), (p-1)*sizeof(Complex<F,L>));

        // Loop variables
        int ak = 1;
        int akinv = 1;
        Complex<F,L> y0 = x[0];

        // First, "invert" the order of x
        for(size_t k = 0; k < (p-1); k++) {
            y[k*s_out] = x[akinv*s_in];
            y0 = y0 + y[k*s_out];
            ak = (ak*a) % p;
            akinv = (akinv*ainv) % p;
        }

        // Convolve the resulting vector with twiddle vector

        // First fft the resulting shuffled vector
        subFFT->fptr(y, z, s_out, 1, subFFT, dir);

        // Perform cyclic convolution
        for(size_t m = 0; m < (p-1); m++) {
            Complex<F,L> Cm = omega<F,L>::get(1, p, dir);
            ak = a;
            for(size_t k = 1; k < (p-1); k++) {
                Cm = Cm + omega<F,L>::get(p*(k*m+ak) - ak, p*(p-1), dir);
                ak = (ak*a) % p;
            }
            y[m*s_out] = z[m]*Cm;
        }

        // Bring back into signal domain
        subFFT->fptr(y, z, s_out, 1, subFFT, (Direction) (-1*((int) dir)));

        // Shuffle as needed
        ak = 1;
        y[0] = y0;
        for(size_t m = 0; m < (p-1); m++) {
            y[ak*s_out] = x[0] + (z[m]/((double) (p-1)));
            ak = (ak*a) % p;
        }
    }

    // Internal-facing radix-N1 C-T FFT for reference. Doesn't use any call graph.
    // Used for checking correctness of output
	template<typename F, int L>
	inline void reference_composite_FFT_helper(size_t N, Complex<F,L>* x, Complex<F,L>* y, size_t s_in, size_t s_out, Direction dir) {
        
        // Find the factors of N
        size_t N1 = referenceFactor(N);
        size_t N2 = N/N1;
        
        // Execute a DFT if necessary
        if(N1 == N) {
            DFT_helper(N, x, y, s_in, s_out, dir);
            return;
        }

        // Malloc a temporary array to hold intermediate results
        auto z = (Complex<F,L>*) aligned_alloc(sizeof(Complex<F,L>), N*sizeof(Complex<F,L>));

        // Take the FFT of the "rows" and put output into z.
        reference_composite_FFT_helper(N2, x, z, N1*s_in, 1, dir);
        for(size_t n1 = 1; n1 < N1; n1++) {
            reference_composite_FFT_helper(N2, x+n1*s_in, z+N2*n1, N1*s_in, 1, dir);
            for(size_t k2 = 1; k2 < N2; k2++) {
                z[n1*N2 + k2] = z[n1*N2 + k2]*omega<F,L>::get(n1*k2, N, dir);
            }
        } 

        // Take the FFT of the "columns" of z and put it into y
        // Don't need n1 for the second transform as it's just the indexer.
        for(size_t k2 = 0; k2 < N2; k2++) {
            reference_composite_FFT_helper(N1, z+k2, y+k2*s_out, N2, N2*s_out, dir);
        }

        // Free memory used by z
        free(z);
    }

    // External-facing function to calculate a composite FFT using the reference code.
	template<typename F, int L>
	inline void reference_composite_FFT(size_t N, Complex<F,L>* x, Complex<F,L>* y, Direction dir) {
        reference_composite_FFT_helper(N, x, y, 1, 1, dir);
    }

    // Internal recursive helper-function that calculates the FFT of a signal with length 3^k
	template<typename F, int L>
	inline void pow3_FFT_helper(size_t N, Complex<F,L>* x, Complex<F,L>* y, size_t s_in, size_t s_out, Direction dir, Complex<F,L>& plus120, Complex<F,L>& minus120) {
        
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
        pow3_FFT_helper(Nprime, x, y, s_in*3, s_out, dir, plus120, minus120);
        pow3_FFT_helper(Nprime, x+s_in, y+Nprime*s_out, s_in*3, s_out, dir, plus120, minus120);
        pow3_FFT_helper(Nprime, x+2*s_in, y+2*Nprime*s_out, s_in*3, s_out, dir, plus120, minus120);

        // Combine the sub-problem solutions
        Complex<F,L> wk1 (1., 0.); Complex<F,L> wk2 (1., 0.);
        double inc = 2.*M_PI/N;
        Complex<F,L> w1 (cos(  inc), static_cast<int>(dir)*sin(  inc));
        Complex<F,L> w2 (cos(2*inc), static_cast<int>(dir)*sin(2*inc));
        for(size_t k = 0; k < Nprime; k++) {
            // Index calculation
            auto k1 = k * s_out;
            auto k2 = (  Nprime + k) * s_out;
            auto k3 = (2*Nprime + k) * s_out;

            // Storing temporary variables
            Complex<F,L> tmpk     = y[k1];
            Complex<F,L> tmpk_p_1 = y[k2];
            Complex<F,L> tmpk_p_2 = y[k3];

            // Reassigning the output
            y[k1] = tmpk +             wk1 * tmpk_p_1 +            wk2 * tmpk_p_2;
            y[k2] = tmpk + plus120  *  wk1 * tmpk_p_1 + minus120 * wk2 * tmpk_p_2;
            y[k3] = tmpk + minus120 *  wk1 * tmpk_p_1 + plus120  * wk2 * tmpk_p_2;

            // Twiddle factors
            wk1 *= w1; wk2 *= w2;
        }
    }

    // External-facing function for performing an FFT on signal with length N = 3^k
	template<typename F, int L>
	inline void pow3_FFT(Complex<F,L>* x, Complex<F,L>* y, size_t s_in, size_t s_out, biFuncNode<F,L>* sRoot, Direction dir) {
        const size_t N = sRoot->sz;
        Complex<F,L> plus120 {-0.5, -sqrt(3)/2.};
        Complex<F,L> minus120 {-0.5, sqrt(3)/2.};
        switch(dir) {
            case Direction::forward: pow3_FFT_helper(N, x, y, s_in, s_out, dir, plus120, minus120); break;
            case Direction::inverse: pow3_FFT_helper(N, x, y, s_in, s_out, dir, minus120, plus120); break;
        }
    }

}