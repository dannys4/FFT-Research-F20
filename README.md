# FFTE: Fast Fourier Transform Engine
## Primarily created by Danny Sharp
## In collaboriation with Miroslav Stoyanov, Jeff Borggaard

This library is primarily designed to be a flexible alternative to libraries like, say, FFTW. This does not necessarily guarantee the speeds of competing FFT libraries, however it is designed to ensure compatibility and usability. It's primarily designed to be compatible with C++11, though there are some nice features regarding `constexpr` when C++14 or higher is used. Currently, though, C++14 and higher does not support the implementation of Rader's algorithm, so beware of large primes in that case.

If there are any questions, concerns, or problems, reach out using the Github Issues for this library.