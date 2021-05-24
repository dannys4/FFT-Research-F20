/**
 * Code Author: Danny Sharp
 * This file is part of STOCK_FFT (Fast Fourier Transform Engine)
 */

#ifndef DIRECTION_HPP
#define DIRECTION_HPP
namespace STOCK_FFT {
    // An enum for knowing whether the transform is forward or inverse
    enum class Direction {forward = -1, inverse = 1};
}
#endif // END DIRECTION_HPP