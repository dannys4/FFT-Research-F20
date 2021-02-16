#include "algos.hpp"
#include "tree.hpp"
#include <iostream>

/**
 * Code Author: Danny Sharp
 * This file is part of FFTE (Fast Fourier Transform Engine)
 */
using namespace FFTE;
// Main function. Used for calling different parts of the testing code
int main() {
    int n = 20;
    auto c = (Complex<double, 2>*) malloc(n*sizeof(Complex<double, 2>));
    for(int i = 0; i < n; i++) {
        c[i] = Complex<double, 2> {2.*i, 3.*i};
    }
    int ell = getNumNodes(n);
    biFuncNode<double, 2> root[ell];
    init_fft_tree(root, n);
}