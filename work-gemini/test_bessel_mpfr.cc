#include "bessel_improved.hh"
#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono>
#include <mpfr.h>
#include <quadmath.h>

using namespace bessel_improved;

// Precision test helper using MPFR
void check_mpfr(double x_val, const char* name, float64x2 (*func)(float64x2 const&), 
                int (*mpfr_func)(mpfr_t, const mpfr_t, mpfr_rnd_t)) {
    mpfr_t mp_x, mp_res;
    mpfr_init2(mp_x, 256);
    mpfr_init2(mp_res, 256);
    mpfr_set_d(mp_x, x_val, MPFR_RNDN);
    mpfr_func(mp_res, mp_x, MPFR_RNDN);

    float64x2 x_dd;
    x_arg_set:
    x_dd._limbs[0] = x_val;
    x_dd._limbs[1] = 0.0;
    
    float64x2 got = func(x_dd);
    
    // Convert DD to string for high precision diff
    // For simplicity, we use __float128 as intermediate
    __float128 got_q = (__float128)got._limbs[0] + (__float128)got._limbs[1];
    
    mpfr_t mp_got, mp_diff;
    mpfr_init2(mp_got, 256);
    mpfr_init2(mp_diff, 256);
    
    char buf[128];
    quadmath_snprintf(buf, sizeof(buf), "%.40Qg", got_q);
    mpfr_set_str(mp_got, buf, 10, MPFR_RNDN);
    
    mpfr_sub(mp_diff, mp_got, mp_res, MPFR_RNDN);
    mpfr_div(mp_diff, mp_diff, mp_res, MPFR_RNDN);
    mpfr_abs(mp_diff, mp_diff, MPFR_RNDN);
    
    double rel_err = mpfr_get_d(mp_diff, MPFR_RNDN);
    
    std::cout << name << "(" << x_val << ") Rel Error vs MPFR: " 
              << std::scientific << std::setprecision(4) << rel_err << std::endl;

    mpfr_clear(mp_x);
    mpfr_clear(mp_res);
    mpfr_clear(mp_got);
    mpfr_clear(mp_diff);
}

// Benchmark helper
template<typename F>
double benchmark(F func, double x_start, int iterations) {
    auto start = std::chrono::high_resolution_clock::now();
    double dummy = 0;
    for (int i = 0; i < iterations; ++i) {
        dummy += (double)func(x_start + (i * 1e-7));
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    return diff.count() / iterations;
}

int main() {
    mpfr_set_default_prec(256);
    
    std::cout << "--- Precision Test vs MPFR (256-bit) ---" << std::endl;
    check_mpfr(1.0,   "j0", j0, mpfr_j0);
    check_mpfr(10.0,  "j0", j0, mpfr_j0);
    check_mpfr(1.0,   "j1", j1, mpfr_j1);
    check_mpfr(10.0,  "j1", j1, mpfr_j1);
    check_mpfr(1.0,   "y0", y0, mpfr_y0);
    check_mpfr(10.0,  "y0", y0, mpfr_y0);
    check_mpfr(1.0,   "y1", y1, mpfr_y1);
    check_mpfr(10.0,  "y1", y1, mpfr_y1);

    std::cout << "\n--- Speed Test: Improved DD vs libquadmath ---" << std::endl;
    const int N = 1000000;
    
    auto bench_j0_dd = [](double x) {
        float64x2 val; val._limbs[0] = x; val._limbs[1] = 0;
        float64x2 res = j0(val);
        return res._limbs[0];
    };
    auto bench_j0_q = [](double x) { return (double)j0q((__float128)x); };

    double t_dd = benchmark(bench_j0_dd, 10.0, N);
    double t_q  = benchmark(bench_j0_q, 10.0, N);

    std::cout << "j0 performance (avg per call):" << std::endl;
    std::cout << "  Improved DD:  " << std::fixed << std::setprecision(2) << t_dd * 1e9 << " ns" << std::endl;
    std::cout << "  libquadmath: " << std::fixed << std::setprecision(2) << t_q * 1e9 << " ns" << std::endl;
    std::cout << "  Speedup:      " << std::fixed << std::setprecision(2) << t_q / t_dd << "x" << std::endl;

    return 0;
}
