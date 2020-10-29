#include "test.hpp"

int main() {
    Omega w{};
    if(w()) std::cout << "w is not null\n";
    else std::cout << "w is null\n";
}