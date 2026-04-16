#include <quadmath.h>
#include <iostream>
int main() {
    __float128 x = 1.0Q;
    char buf[128];
    quadmath_snprintf(buf, sizeof(buf), "%.40Qg", x);
    std::cout << buf << std::endl;
    return 0;
}
