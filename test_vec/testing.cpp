#include <iostream>
#include <chrono>
#include <random>
#include <functional>
// #include <vector>
// #include <array>

#ifndef LENGTH
#define LENGTH 25
#endif

#define FACTORS_LEN 15
static constexpr uint64_t factors[FACTORS_LEN] {4, 2, 3, 5, 7, 11, 13, 16, 17, 19, 23, 29, 31, 37, 41};

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

constexpr uint64_t initConstBiTree(constBiNode<uint64_t>* sRoot, const uint64_t N) {
    uint64_t k = numNodesFactorHelper(N);
    (*sRoot).elem = N;
    if (k == N) return 1;

    uint64_t q = getLeftover(N, k);
    uint64_t l = initConstBiTree(sRoot + 1, k);
    uint64_t r = initConstBiTree(sRoot + 1 + l, q);
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

void time_compilation() {
    using std::chrono::duration_cast;
    using std::chrono::nanoseconds;
    typedef std::chrono::high_resolution_clock clock;

    uint N = 1e5;
    std::cerr << "Testing compile time performance...\n";
    auto start = clock::now();
    for(uint64_t i = 0; i < N; i++) {
        uint64_t ell = i*i - 3*i + 4;
        uint64_t k = getNumNodes(ell);
        constBiNode<uint64_t> root[k];
        initConstBiTree(root, ell);
        std::cout << "ell = " << ell << ", ";
        printRoot(root, k);
    }
    auto end = clock::now();
    std::cerr << "Known at compile time took " << duration_cast<nanoseconds>(end-start).count() << "ns\n";
    std::cerr << "Initializing random ints...\n";
    auto rand = std::bind(std::uniform_int_distribution<>{2, 4}, std::default_random_engine{});
    std::vector<uint> v {};
    for(uint64_t i = 0; i < N; i++) {
        v.push_back(rand());
    }
    std::cerr << "Random ints initialized. Testing runtime performance...\n";
    start = clock::now();
    for(uint64_t i = 0; i < N; i++) {
        uint64_t ell = i*i - v[i]*i + 5;
        std::cout << "ell = " << ell << ", ";
        uint64_t k = getNumNodes(ell);
        constBiNode<uint64_t> root[k];
        initConstBiTree(root, ell);
        printRoot(root, k);
    }
    end = clock::now();
    std::cerr << "Known at run time took " << duration_cast<nanoseconds>(end-start).count() << "ns\n";
}

int main() {
    const uint64_t N = 135135;
    std::cout << "N = " << N << "\n";
    uint64_t k = getNumNodes(N);
    constBiNode<uint64_t> root[k];
    initConstBiTree(root, N);
    printRoot(root, k);
    printTree(root);
    std::cout << "\n";
    // time_compilation();
}