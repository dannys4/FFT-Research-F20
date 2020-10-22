#ifndef CONSTTREE_HPP
#define CONSTTREE_HPP

#include <iostream>
#include "algos.hpp"

#define FACTORS_LEN 15
static constexpr uint64_t factors[FACTORS_LEN] {4, 2, 3, 5, 7, 11, 13, 16, 17, 19, 23, 29, 31, 37, 41};

// This is the standard fft function call:
// N, x, y, s_in, s_out.
// N     = input signal length
// x     = input signal
// y     = output signal
// s_in  = the distance (stride) between each input
// s_out = the distance (stride) between each output
using fft_fptr = void (*)(uint64_t, Complex*, Complex*, uint64_t, uint64_t);

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

class constBiFuncNode {
    public:
        fft_fptr fptr = nullptr;
        uint64_t sz = 0;
        uint64_t left = 0;
        uint64_t right = 0;
        constexpr constBiFuncNode() = default;
};

// We first check if N is a power of something.
// Then, we check if it has powers of numbers
// Then, if that doesn't work, we just factor it.
constexpr fft_fptr fptrFactorHelper(const uint64_t N, uint64_t* k) {
    if(power_of(N, 4) ||
       power_of(N, 2) ||
       power_of(N, 3)) {
        return N;
    }
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

constexpr fft_fptr init_fft_node(const uint64_t N, const uint64_t k) {
    if()
}

constexpr uint64_t init_fft_tree(constBiFuncNode* sRoot, const uint64_t N) {
    uint64_t k = numNodesFactorHelper(N);
    (*sRoot).sz = N;
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

constexpr uint64_t factor(const uint64_t f) {
    uint64_t k = 0;
    for(; k < FACTORS_LEN; k++) {
        if( f % factors[k] == 0) return factors[k];
    }
    for(k = factors[k - 1]; k*k < f; k+=2) {
        if( f % k == 0 ) return k;
    }
    return f;
}

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

void printRoot(constBiNode<uint64_t>* root, uint64_t N) {
    std::cout << "root[" << N << "] = ";
    for(uint64_t i = 0; i < N; i++) std::cout << root[i].elem << ", ";
    std::cout << "\n";
}

void printTree(constBiNode<uint64_t>* root) {
    std::cout << root->elem;
    if(!(root->left || root->right)) return;
    std::cout << ": (";
    if(root->left)  printTree(root + root->left);
    std::cout << ", ";
    if(root->right) printTree(root + root->right);
    std::cout << ")";
}
#endif