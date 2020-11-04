#ifndef CONSTTREE_HPP
#define CONSTTREE_HPP
/**
 * Code Author: Danny Sharp
 * This file is part of FFTE (Fast Fourier Transform Engine)
 */


#include <iostream>
#include "algos.hpp"


/* This is a generic class designed to be used as a
 * compile-time ready tree. It is not used anywhere here
 * except for debugging purposes
 */
template<typename T>
class constBiNode {
    private:
    public:
        T elem;
        uint64_t left = 0;
        uint64_t right = 0;
        constexpr constBiNode() = default;
        constexpr constBiNode(const T e) {elem = e;}
};

constexpr bool power_of(const uint64_t n, const uint64_t pow) {
    uint64_t k = n;
    if(pow == 2) {
        uint8_t sum = 0x1 & n;
        for(;k;k>>=1) sum += 0x1&k;
        return sum == 1;
    }
    while(!(k % pow)) k/= pow;
    return k == 1;
}

// We first check if N is a power of something.
// Then, we check if it has powers of numbers
// Then, if that doesn't work, we just factor it.
constexpr uint64_t numNodesFactorHelper(const uint64_t N) {
    if(power_of(N, 4) ||
       power_of(N, 2) ||
       power_of(N, 3)) {
        return N;
    }
    if((N % 4) == 0) {
        uint64_t k = N;
        while ((k % 4) == 0) k /= 4;
        return N / k;
    }
    if((N % 2) == 0) {
        uint64_t k = N;
        while ((k % 2) == 0) k /= 2;
        return N / k;
    }
    if((N % 3) == 0) {
        uint64_t k = N;
        while ((k % 3) == 0) k /= 3;
        return N / k;
    }
    return factor(N);
}

// This is a placeholder for what we need when
// prime algorithms are introduced
constexpr uint64_t getLeftover(const uint64_t N, const uint64_t k) {
    return N / k;
}

// We first check if N is a power of something.
// Then, we check if it has powers of numbers
// Then, if that doesn't work, we just factor it.
constexpr fft_fptr fptrFactorHelper(const uint64_t N, uint64_t* k) {
    if(power_of(N, 4)) {
        *k = N;
        return pow2_FFT;
    }
    if(power_of(N, 2)) {
        *k = N;
        return pow2_FFT;
    }
    if(power_of(N, 3)) {
        *k = N;
        return pow3_FFT;
    }
    if((N % 4) == 0) {
        uint64_t ell = N;
        while ((ell % 4) == 0) ell /= 4;
        *k = N/ell;
        return composite_FFT;
    }
    if((N % 2) == 0) {
        uint64_t ell = N;
        while ((ell % 2) == 0) ell /= 2;
        *k = N/ell;
        return composite_FFT;
    }
    if((N % 3) == 0) {
        uint64_t ell = N;
        while ((ell % 3) == 0) ell /= 3;
        *k = N/ell;
        return composite_FFT;
    }
    *k = factor(N);
    return (*k == N) ? DFT : composite_FFT;
}

constexpr uint64_t init_fft_tree(constBiFuncNode* sRoot, const uint64_t N) {
    uint64_t k = 0;
    fft_fptr fptr = fptrFactorHelper(N, &k);
    sRoot->sz = N;
    sRoot->fptr = fptr;
    if (k == N) {
        return 1;
    }
    uint64_t q = getLeftover(N, k);
    uint64_t l = init_fft_tree(sRoot + 1, k);
    uint64_t r = init_fft_tree(sRoot + 1 + l, q);
    sRoot->left = 1;
    sRoot->right = 1 + l;
    return 1 + l + r;

}

constexpr uint64_t getNumNodes(const uint64_t N) {
    uint64_t k = numNodesFactorHelper(N);
    if(k == N) return 1;
    return 1 + getNumNodes(k) + getNumNodes(getLeftover(N, k));
}

constexpr uint64_t initUintConstBiTree(constBiNode<uint64_t>* sRoot, const uint64_t N) {
    uint64_t k = numNodesFactorHelper(N);
    (*sRoot).elem = N;
    if (k == N) return 1;

    uint64_t q = getLeftover(N, k);
    uint64_t l = initUintConstBiTree(sRoot + 1, k);
    uint64_t r = initUintConstBiTree(sRoot + 1 + l, q);
    sRoot->left = 1;
    sRoot->right = 1 + l;
    return 1 + l + r;
}

inline void printRoot(constBiNode<uint64_t>* root, uint64_t N) {
    std::cout << "root[" << N << "] = ";
    for(uint64_t i = 0; i < N; i++) std::cout << root[i].elem << ", ";
    std::cout << "\n";
}

inline void printTree(constBiNode<uint64_t>* root) {
    std::cout << root->elem;
    if(!(root->left || root->right)) return;
    std::cout << ": (";
    if(root->left)  printTree(root + root->left);
    std::cout << ", ";
    if(root->right) printTree(root + root->right);
    std::cout << ")";
}

inline void printTree(constBiFuncNode* root) {
    std::cout << root->sz;
    if(!(root->left || root->right)) return;
    std::cout << ": (";
    if(root->left)  printTree(root + root->left);
    std::cout << ", ";
    if(root->right) printTree(root + root->right);
    std::cout << ")";
}
#endif