#include "bessel_improved.hh"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <quadmath.h>

using namespace bessel_improved;

__float128 to_quad(float64x2 const &x) {
    return (__float128)x._limbs[0] + (__float128)x._limbs[1];
}

void check(__float128 expected, float64x2 const &got, const char* name) {
    __float128 got_q = to_quad(got);
    __float128 diff = got_q - expected;
    if (diff < 0) diff = -diff;
    __float128 abs_expected = expected;
    if (abs_expected < 0) abs_expected = -abs_expected;
    
    double rel_err = (double)(diff / abs_expected);
    
    char buf_exp[128], buf_got[128];
    quadmath_snprintf(buf_exp, sizeof(buf_exp), "%.35Qg", expected);
    quadmath_snprintf(buf_got, sizeof(buf_got), "%.35Qg", got_q);

    std::cout << name << "(10.0):" << std::endl;
    std::cout << "  Expected: " << buf_exp << std::endl;
    std::cout << "  Got:      " << buf_got << std::endl;
    std::cout << "  Rel error: " << std::scientific << std::setprecision(4) << rel_err << std::endl;
}

int main() {
    std::cout << "Testing Improved C++ Bessel functions vs libquadmath reference" << std::endl;
    
    float64x2 x10, x1;
    x10._limbs[0] = 10.0; x10._limbs[1] = 0.0;
    x1._limbs[0] = 1.0; x1._limbs[1] = 0.0;

    std::cout << "--- x = 10.0 ---" << std::endl;
    check(j0q(10.0Q), j0(x10), "j0");
    check(j1q(10.0Q), j1(x10), "j1");
    check(y0q(10.0Q), y0(x10), "y0");
    check(y1q(10.0Q), y1(x10), "y1");

    std::cout << "\n--- x = 1.0 ---" << std::endl;
    check(j0q(1.0Q), j0(x1), "j0");
    check(j1q(1.0Q), j1(x1), "j1");
    check(y0q(1.0Q), y0(x1), "y0");
    check(y1q(1.0Q), y1(x1), "y1");

    float64x2 x05, x100;
    x05._limbs[0] = 0.5; x05._limbs[1] = 0.0;
    x100._limbs[0] = 100.0; x100._limbs[1] = 0.0;

    std::cout << "\n--- x = 0.5 ---" << std::endl;
    check(j0q(0.5Q), j0(x05), "j0");
    check(j1q(0.5Q), j1(x05), "j1");
    check(y0q(0.5Q), y0(x05), "y0");
    check(y1q(0.5Q), y1(x05), "y1");

    std::cout << "\n--- x = 100.0 ---" << std::endl;
    check(j0q(100.0Q), j0(x100), "j0");
    check(j1q(100.0Q), j1(x100), "j1");
    check(y0q(100.0Q), y0(x100), "y0");
    check(y1q(100.0Q), y1(x100), "y1");

    return 0;
}
