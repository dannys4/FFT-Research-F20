/**
 * Code Author: Danny Sharp
 * This file is part of the implementation for a stock FFT algorithm intended for HeFFTe
 */

#ifndef DIRECTION_HPP
#define DIRECTION_HPP
namespace stock_fft {
    // An enum for knowing whether the transform is forward or inverse
    enum class Direction {forward = -1, inverse = 1};
}
#endif // END DIRECTION_HPP