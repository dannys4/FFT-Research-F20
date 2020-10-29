#include "test.hpp"

void helper(Complex& n) {
    if(n == NULL) {
        std::cout << "successful!\n";
    }
    else
    {
        std::cout << "unsuccessful :(\n";
    }
    
}

int main() {
    time_omega();
}