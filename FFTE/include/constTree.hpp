/* Header file for C++14 and up */
#ifndef STOCK_FFT_TREE_HPP
#define STOCK_FFT_TREE_HPP
/**
 * Code Author: Danny Sharp
 * This file is part of STOCK_FFT (Fast Fourier Transform Engine)
 */

// We use this because it's a variable used in the factor function, so by defining it beforehand,
// going out of bounds on memory is much harder
#define STOCK_FFT_FACTORS_LEN 15 

#include <iostream>
#include "algos.hpp"
namespace STOCK_FFT {

    /* Statically allocated array of factors that are known at compile-time. These
    * are not necessarily prime, just ordered in the way that we prioritize.
    */
    static constexpr size_t factors[STOCK_FFT_FACTORS_LEN] {4, 2, 3, 5, 7, 11, 13, 16, 17, 19, 23, 29, 31, 37, 41};

    // Function to find the smallest usable factor of f at compile-time.
    constexpr size_t factor(const size_t f) {
        size_t k = 0;
        // Prioritize factors in the factors array
        for(; k < STOCK_FFT_FACTORS_LEN; k++) {
            if( f % factors[k] == 0) return factors[k];
        }

        // If none of those work, turn to odd numbers that are greater than the last index
        for(k = factors[k - 1]; k*k < f; k+=2) {
            if( f % k == 0 ) return k;
        }

        // return f if no factor was found
        return f;
    }

    /* This is a generic class designed to be used as a
    * compile-time ready tree. It is not used anywhere here
    * except for debugging purposes
    */
    template<typename T>
    class constBiNode {
        private:
        public:
            T elem;
            size_t left = 0;
            size_t right = 0;
            constexpr constBiNode() = default;
            constexpr constBiNode(const T e) {elem = e;}
    };

    // Check if n is a power of pow
    constexpr bool power_of(const size_t n, const size_t pow) {
        size_t k = n;
        if(pow == 2) {
            unsigned char sum = 0x1 & n;
            for(;k;k>>=1) sum += 0x1&k;
            return sum == 1;
        }
        while(!(k % pow)) k/= pow;
        return k == 1;
    }

    // We first check if N is a power of something.
    // Then, we check if it has powers of numbers
    // Then, if that doesn't work, we just factor it.
    constexpr size_t numNodesFactorHelper(const size_t N) {
        if(power_of(N, 4) ||
        power_of(N, 2) ||
        power_of(N, 3)) {
            return N;
        }
        if((N % 4) == 0) {
            size_t k = N;
            while ((k % 4) == 0) k /= 4;
            return N / k;
        }
        if((N % 2) == 0) {
            size_t k = N;
            while ((k % 2) == 0) k /= 2;
            return N / k;
        }
        if((N % 3) == 0) {
            size_t k = N;
            while ((k % 3) == 0) k /= 3;
            return N / k;
        }
        return factor(N);
    }

    // This is a placeholder for what we need when
    // prime algorithms are introduced
    constexpr size_t getLeftover(const size_t N, const size_t k) {
        return N / k;
    }

    // We first check if N is a power of something.
    // Then, we check if it has powers of numbers
    // Then, if that doesn't work, we just factor it.
    constexpr fft_type fptrFactorHelper(const size_t N, size_t* k) {
        if(power_of(N, 4)) {
            *k = N;
            return fft_type::pow4;
        }
        if(power_of(N, 2)) {
            *k = N;
            return fft_type::pow2;
        }
        if(power_of(N, 3)) {
            *k = N;
            return fft_type::pow3;
        }
        if((N % 4) == 0) {
            size_t ell = N;
            while ((ell % 4) == 0) ell /= 4;
            *k = N/ell;
            return fft_type::composite;
        }
        if((N % 2) == 0) {
            size_t ell = N;
            while ((ell % 2) == 0) ell /= 2;
            *k = N/ell;
            return fft_type::composite;
        }
        if((N % 3) == 0) {
            size_t ell = N;
            while ((ell % 3) == 0) ell /= 3;
            *k = N/ell;
            return fft_type::composite;
        }
        *k = factor(N);
        return (*k == N) ? fft_type::discrete : fft_type::composite;
    }

    /* Initialize an fft tree given an appropriately sized empty array of 
    * function nodes to hold the information
    */
    constexpr size_t init_fft_tree(biFuncNode* sRoot, const size_t N) {
        size_t k = 0;
        fft_type fptr = fptrFactorHelper(N, &k);
        sRoot->sz = N;
        sRoot->fptr = Fourier_Transform(fptr);
        if (k == N) {
            return 1;
        }
        size_t q = getLeftover(N, k);
        size_t l = init_fft_tree(sRoot + 1, k);
        size_t r = init_fft_tree(sRoot + 1 + l, q);
        sRoot->left = 1;
        sRoot->right = 1 + l;
        return 1 + l + r;

    }

    // Get the number of nodes appropriate for a given length signal, N
    constexpr size_t getNumNodes(const size_t N) {
        size_t k = numNodesFactorHelper(N);
        if(k == N) return 1;
        return 1 + getNumNodes(k) + getNumNodes(getLeftover(N, k));
    }

    // Initialize an unsigned integer tree (useful for debugging the tree construction)
    constexpr size_t initUintConstBiTree(constBiNode<size_t>* sRoot, const size_t N) {
        size_t k = numNodesFactorHelper(N);
        (*sRoot).elem = N;
        if (k == N) return 1;

        size_t q = getLeftover(N, k);
        size_t l = initUintConstBiTree(sRoot + 1, k);
        size_t r = initUintConstBiTree(sRoot + 1 + l, q);
        sRoot->left = 1;
        sRoot->right = 1 + l;
        return 1 + l + r;
    }

    // Print the nodes of a size_t tree as they're stored in the length N array root
    inline void printRoot(constBiNode<size_t>* root, size_t N) {
        std::cout << "root[" << N << "] = ";
        for(size_t i = 0; i < N; i++) std::cout << root[i].elem << ", ";
        std::cout << "\n";
    }

    // Print the nodes using a pre-order traversal
    inline void printTree(constBiNode<size_t>* root) {
        std::cout << root->elem;
        if(!(root->left || root->right)) return;
        std::cout << ": (";
        if(root->left)  printTree(root + root->left);
        std::cout << ", ";
        if(root->right) printTree(root + root->right);
        std::cout << ")";
    }

    // print the nodes of an FFT tree using a pre-order traversal
    inline void printTree(biFuncNode* root) {
        std::cout << root->sz;
        if(!(root->left || root->right)) return;
        std::cout << ": (";
        if(root->left)  printTree(root + root->left);
        std::cout << ", ";
        if(root->right) printTree(root + root->right);
        std::cout << ")";
    }
}

#endif // END STOCK_FFT_TREE_HPP