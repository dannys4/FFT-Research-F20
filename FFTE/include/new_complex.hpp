#include <type_traits>
#include <iostream>

template<typename T, class = typename std::enable_if<std::is_floating_point<T>::value>::type>
class mytype {
    public:
        T number;
        mytype(T n): number(n) {}
};