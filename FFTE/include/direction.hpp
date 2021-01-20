/**
 * Code Author: Danny Sharp
 * This file is part of FFTE (Fast Fourier Transform Engine)
 */

#ifndef DIRECTION_HPP
#define DIRECTION_HPP
namespace FFTE {
    // An enum for knowing whether the transform is forward or inverse
    enum class Direction {forward = -1, inverse = 1};
}
#endif // END DIRECTION_HPP