#include <type_traits>
#include <iostream>



template<typename Floating, int N = 1, class = typename std::enable_if<std::is_floating_point<Floating>::value>::type>
class mytype {
    public:
        Floating number;
        mytype(Floating n): number(n) {}
};