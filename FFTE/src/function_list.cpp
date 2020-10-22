#include <vector>
#include "complex.hpp"

// void composite_FFT(uint64_t N, Complex* x, Complex* y, uint64_t s_in, uint64_t s_out)
class FunctionList {
    /*
    * Current format for an FFT:
    * void general_fft(uint64_t N, Complex* x, Complex* y, uint64_t s_in, uint64_t s_out)
    * N: Length of signal
    * x: input signal in time domain
    * y: output signal in frequency domain
    * s_in: how far apart each input element is
    * s_out: how far apart each output element should be placed
    */
   using int_type = uint64_t;
   using f_ptr = void (*)(int_type, Complex*, Complex*, int_type, int_type);
   
};