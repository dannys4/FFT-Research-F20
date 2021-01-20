/* Header file for C++11 and under */
#ifndef FFTE_TREE_HPP
#define FFTE_TREE_HPP
/**
 * Code Author: Danny Sharp
 * This file is part of FFTE (Fast Fourier Transform Engine)
 */

#include "algos.hpp"
#include <iostream>

// We use this because it's a variable used in the factor function, so by defining it beforehand,
// going out of bounds on memory is much harder
#define FFTE_FACTORS_LEN 200

// This is the minimum signal size for us to use Rader's algorithm
#define FFTE_RADER_MIN 1000000

namespace FFTE {

    /* Initialize an fft tree given an appropriately sized empty array of 
    * function nodes to hold the information
    */
    size_t init_fft_tree(biFuncNode* sRoot, const size_t N);

    // Get the number of nodes appropriate for a given length signal, N
    size_t getNumNodes(const size_t N);

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
            constBiNode() = default;
            constBiNode(const T e) {elem = e;}
    };

    // Print the nodes of a size_t tree as they're stored in the length N array root
    void printRoot(constBiNode<size_t>* root, size_t N);

    // Print the nodes using a pre-order traversal
    void printTree(constBiNode<size_t>* root);

    // print the nodes of an FFT tree using a pre-order traversal
    void printTree(biFuncNode* root);

    // Initialize an unsigned integer tree (useful for debugging the tree construction)
    size_t initUintConstBiTree(constBiNode<size_t>* sRoot, const size_t N);
}
#endif