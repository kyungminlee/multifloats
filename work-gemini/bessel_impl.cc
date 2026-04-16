#include "bessel_improved.hh"
#include "bessel_constants.hh"
#include <cmath>

namespace bessel_improved {

using MFD2 = float64x2;

inline MFD2 dd_pair(double hi, double lo) {
    MFD2 r;
    r._limbs[0] = hi;
    r._limbs[1] = lo;
    return r;
}

inline MFD2 dd_neval(MFD2 const &x, const double *hi, const double *lo, int n) {
    MFD2 p = dd_pair(hi[n], lo[n]);
    for (int i = n - 1; i >= 0; --i) {
        p = p * x + dd_pair(hi[i], lo[i]);
    }
    return p;
}

inline MFD2 dd_deval(MFD2 const &x, const double *hi, const double *lo, int n) {
    MFD2 p = x + dd_pair(hi[n], lo[n]);
    for (int i = n - 1; i >= 0; --i) {
        p = p * x + dd_pair(hi[i], lo[i]);
    }
    return p;
}

template<std::size_t N>
inline MFD2 dd_neval_arr(MFD2 const &x, const double (&hi)[N], const double (&lo)[N]) {
    return dd_neval(x, hi, lo, N - 1);
}

template<std::size_t N>
inline MFD2 dd_deval_arr(MFD2 const &x, const double (&hi)[N], const double (&lo)[N]) {
    return dd_deval(x, hi, lo, N - 1);
}

float64x2 j0(float64x2 const &x) {
    MFD2 abs_x = (x._limbs[0] < 0.0) ? -x : x;
    double xx = abs_x._limbs[0];
    if (xx == 0.0) return MFD2(1.0);

    if (xx <= 2.0) {
        MFD2 z = abs_x * abs_x;
        MFD2 p = z * z * dd_neval_arr(z, j0_J0_2N_hi, j0_J0_2N_lo) /
                 dd_deval_arr(z, j0_J0_2D_hi, j0_J0_2D_lo);
        return p - z * 0.25 + MFD2(1.0);
    }

    MFD2 angle = abs_x - dd_pair(PI_QUARTER_hi, PI_QUARTER_lo);
    MFD2 s = multifloats::sin(angle);
    MFD2 c = multifloats::cos(angle);

    MFD2 xinv = MFD2(1.0) / abs_x;
    MFD2 z = xinv * xinv;
    double xinv_d = 1.0 / xx;
    MFD2 p, q;
    
    if (xinv_d <= 0.25) {
        if (xinv_d <= 0.125) {
            if (xinv_d <= 0.0625) {
                p = dd_neval_arr(z, j0_P16_IN_hi, j0_P16_IN_lo) / dd_deval_arr(z, j0_P16_ID_hi, j0_P16_ID_lo);
                q = dd_neval_arr(z, j0_Q16_IN_hi, j0_Q16_IN_lo) / dd_deval_arr(z, j0_Q16_ID_hi, j0_Q16_ID_lo);
            } else {
                p = dd_neval_arr(z, j0_P8_16N_hi, j0_P8_16N_lo) / dd_deval_arr(z, j0_P8_16D_hi, j0_P8_16D_lo);
                q = dd_neval_arr(z, j0_Q8_16N_hi, j0_Q8_16N_lo) / dd_deval_arr(z, j0_Q8_16D_hi, j0_Q8_16D_lo);
            }
        } else if (xinv_d <= 0.1875) {
            p = dd_neval_arr(z, j0_P5_8N_hi, j0_P5_8N_lo) / dd_deval_arr(z, j0_P5_8D_hi, j0_P5_8D_lo);
            q = dd_neval_arr(z, j0_Q5_8N_hi, j0_Q5_8N_lo) / dd_deval_arr(z, j0_Q5_8D_hi, j0_Q5_8D_lo);
        } else {
            p = dd_neval_arr(z, j0_P4_5N_hi, j0_P4_5N_lo) / dd_deval_arr(z, j0_P4_5D_hi, j0_P4_5D_lo);
            q = dd_neval_arr(z, j0_Q4_5N_hi, j0_Q4_5N_lo) / dd_deval_arr(z, j0_Q4_5D_hi, j0_Q4_5D_lo);
        }
    } else {
        if (xinv_d <= 0.375) {
            if (xinv_d <= 0.3125) {
                p = dd_neval_arr(z, j0_P3r2_4N_hi, j0_P3r2_4N_lo) / dd_deval_arr(z, j0_P3r2_4D_hi, j0_P3r2_4D_lo);
                q = dd_neval_arr(z, j0_Q3r2_4N_hi, j0_Q3r2_4N_lo) / dd_deval_arr(z, j0_Q3r2_4D_hi, j0_Q3r2_4D_lo);
            } else {
                p = dd_neval_arr(z, j0_P2r7_3r2N_hi, j0_P2r7_3r2N_lo) / dd_deval_arr(z, j0_P2r7_3r2D_hi, j0_P2r7_3r2D_lo);
                q = dd_neval_arr(z, j0_Q2r7_3r2N_hi, j0_Q2r7_3r2N_lo) / dd_deval_arr(z, j0_Q2r7_3r2D_hi, j0_Q2r7_3r2D_lo);
            }
        } else if (xinv_d <= 0.4375) {
            p = dd_neval_arr(z, j0_P2r3_2r7N_hi, j0_P2r3_2r7N_lo) / dd_deval_arr(z, j0_P2r3_2r7D_hi, j0_P2r3_2r7D_lo);
            q = dd_neval_arr(z, j0_Q2r3_2r7N_hi, j0_Q2r3_2r7N_lo) / dd_deval_arr(z, j0_Q2r3_2r7D_hi, j0_Q2r3_2r7D_lo);
        } else {
            p = dd_neval_arr(z, j0_P2_2r3N_hi, j0_P2_2r3N_lo) / dd_deval_arr(z, j0_P2_2r3D_hi, j0_P2_2r3D_lo);
            q = dd_neval_arr(z, j0_Q2_2r3N_hi, j0_Q2_2r3N_lo) / dd_deval_arr(z, j0_Q2_2r3D_hi, j0_Q2_2r3D_lo);
        }
    }
    p = MFD2(1.0) + z * p;
    q = (z * q - 0.125) * xinv;
    MFD2 tpi = dd_pair(TWO_OVER_PI_hi, TWO_OVER_PI_lo);
    return multifloats::sqrt(tpi / abs_x) * (p * c - q * s);
}

float64x2 j1(float64x2 const &x) {
    bool neg = x._limbs[0] < 0.0;
    MFD2 abs_x = neg ? -x : x;
    double xx = abs_x._limbs[0];
    if (xx == 0.0) return MFD2(0.0);

    MFD2 res;
    if (xx <= 2.0) {
        MFD2 z = abs_x * abs_x;
        res = abs_x * 0.5 + abs_x * z * dd_neval_arr(z, j1_J1_2N_hi, j1_J1_2N_lo) /
                           dd_deval_arr(z, j1_J1_2D_hi, j1_J1_2D_lo);
    } else {
        MFD2 angle = abs_x - dd_pair(THREE_PI_QUARTER_hi, THREE_PI_QUARTER_lo);
        MFD2 s = multifloats::sin(angle);
        MFD2 c = multifloats::cos(angle);

        MFD2 xinv = MFD2(1.0) / abs_x;
        MFD2 z = xinv * xinv;
        double xinv_d = 1.0 / xx;
        MFD2 p, q;
        if (xinv_d <= 0.25) {
            if (xinv_d <= 0.125) {
                if (xinv_d <= 0.0625) {
                    p = dd_neval_arr(z, j1_P16_IN_hi, j1_P16_IN_lo) / dd_deval_arr(z, j1_P16_ID_hi, j1_P16_ID_lo);
                    q = dd_neval_arr(z, j1_Q16_IN_hi, j1_Q16_IN_lo) / dd_deval_arr(z, j1_Q16_ID_hi, j1_Q16_ID_lo);
                } else {
                    p = dd_neval_arr(z, j1_P8_16N_hi, j1_P8_16N_lo) / dd_deval_arr(z, j1_P8_16D_hi, j1_P8_16D_lo);
                    q = dd_neval_arr(z, j1_Q8_16N_hi, j1_Q8_16N_lo) / dd_deval_arr(z, j1_Q8_16D_hi, j1_Q8_16D_lo);
                }
            } else if (xinv_d <= 0.1875) {
                p = dd_neval_arr(z, j1_P5_8N_hi, j1_P5_8N_lo) / dd_deval_arr(z, j1_P5_8D_hi, j1_P5_8D_lo);
                q = dd_neval_arr(z, j1_Q5_8N_hi, j1_Q5_8N_lo) / dd_deval_arr(z, j1_Q5_8D_hi, j1_Q5_8D_lo);
            } else {
                p = dd_neval_arr(z, j1_P4_5N_hi, j1_P4_5N_lo) / dd_deval_arr(z, j1_P4_5D_hi, j1_P4_5D_lo);
                q = dd_neval_arr(z, j1_Q4_5N_hi, j1_Q4_5N_lo) / dd_deval_arr(z, j1_Q4_5D_hi, j1_Q4_5D_lo);
            }
        } else {
            if (xinv_d <= 0.375) {
                if (xinv_d <= 0.3125) {
                    p = dd_neval_arr(z, j1_P3r2_4N_hi, j1_P3r2_4N_lo) / dd_deval_arr(z, j1_P3r2_4D_hi, j1_P3r2_4D_lo);
                    q = dd_neval_arr(z, j1_Q3r2_4N_hi, j1_Q3r2_4N_lo) / dd_deval_arr(z, j1_Q3r2_4D_hi, j1_Q3r2_4D_lo);
                } else {
                    p = dd_neval_arr(z, j1_P2r7_3r2N_hi, j1_P2r7_3r2N_lo) / dd_deval_arr(z, j1_P2r7_3r2D_hi, j1_P2r7_3r2D_lo);
                    q = dd_neval_arr(z, j1_Q2r7_3r2N_hi, j1_Q2r7_3r2N_lo) / dd_deval_arr(z, j1_Q2r7_3r2D_hi, j1_Q2r7_3r2D_lo);
                }
            } else if (xinv_d <= 0.4375) {
                p = dd_neval_arr(z, j1_P2r3_2r7N_hi, j1_P2r3_2r7N_lo) / dd_deval_arr(z, j1_P2r3_2r7D_hi, j1_P2r3_2r7D_lo);
                q = dd_neval_arr(z, j1_Q2r3_2r7N_hi, j1_Q2r3_2r7N_lo) / dd_deval_arr(z, j1_Q2r3_2r7D_hi, j1_Q2r3_2r7D_lo);
            } else {
                p = dd_neval_arr(z, j1_P2_2r3N_hi, j1_P2_2r3N_lo) / dd_deval_arr(z, j1_P2_2r3D_hi, j1_P2_2r3D_lo);
                q = dd_neval_arr(z, j1_Q2_2r3N_hi, j1_Q2_2r3N_lo) / dd_deval_arr(z, j1_Q2_2r3D_hi, j1_Q2_2r3D_lo);
            }
        }
        p = MFD2(1.0) + z * p;
        q = (z * q + 0.375) * xinv;
        MFD2 tpi = dd_pair(TWO_OVER_PI_hi, TWO_OVER_PI_lo);
        res = multifloats::sqrt(tpi / abs_x) * (p * c - q * s);
    }
    return neg ? -res : res;
}

float64x2 y0(float64x2 const &x) {
    if (x._limbs[0] <= 0.0) return MFD2(-INFINITY);
    double xx = x._limbs[0];

    if (xx <= 1e-17) {
        MFD2 tpi = dd_pair(TWO_OVER_PI_hi, TWO_OVER_PI_lo);
        return dd_pair(U0_hi, U0_lo) + tpi * multifloats::log(x);
    }

    if (xx <= 2.0) {
        MFD2 z = x * x;
        MFD2 p = dd_neval_arr(z, j0_Y0_2N_hi, j0_Y0_2N_lo) /
                 dd_deval_arr(z, j0_Y0_2D_hi, j0_Y0_2D_lo);
        MFD2 tpi = dd_pair(TWO_OVER_PI_hi, TWO_OVER_PI_lo);
        return tpi * multifloats::log(x) * j0(x) + p;
    }

    MFD2 angle = x - dd_pair(PI_QUARTER_hi, PI_QUARTER_lo);
    MFD2 s = multifloats::sin(angle);
    MFD2 c = multifloats::cos(angle);

    MFD2 xinv = MFD2(1.0) / x;
    MFD2 z = xinv * xinv;
    double xinv_d = 1.0 / xx;
    MFD2 p, q;
    // same p, q as j0 for x > 2
    if (xinv_d <= 0.25) {
        if (xinv_d <= 0.125) {
            if (xinv_d <= 0.0625) {
                p = dd_neval_arr(z, j0_P16_IN_hi, j0_P16_IN_lo) / dd_deval_arr(z, j0_P16_ID_hi, j0_P16_ID_lo);
                q = dd_neval_arr(z, j0_Q16_IN_hi, j0_Q16_IN_lo) / dd_deval_arr(z, j0_Q16_ID_hi, j0_Q16_ID_lo);
            } else {
                p = dd_neval_arr(z, j0_P8_16N_hi, j0_P8_16N_lo) / dd_deval_arr(z, j0_P8_16D_hi, j0_P8_16D_lo);
                q = dd_neval_arr(z, j0_Q8_16N_hi, j0_Q8_16N_lo) / dd_deval_arr(z, j0_Q8_16D_hi, j0_Q8_16D_lo);
            }
        } else if (xinv_d <= 0.1875) {
            p = dd_neval_arr(z, j0_P5_8N_hi, j0_P5_8N_lo) / dd_deval_arr(z, j0_P5_8D_hi, j0_P5_8D_lo);
            q = dd_neval_arr(z, j0_Q5_8N_hi, j0_Q5_8N_lo) / dd_deval_arr(z, j0_Q5_8D_hi, j0_Q5_8D_lo);
        } else {
            p = dd_neval_arr(z, j0_P4_5N_hi, j0_P4_5N_lo) / dd_deval_arr(z, j0_P4_5D_hi, j0_P4_5D_lo);
            q = dd_neval_arr(z, j0_Q4_5N_hi, j0_Q4_5N_lo) / dd_deval_arr(z, j0_Q4_5D_hi, j0_Q4_5D_lo);
        }
    } else {
        if (xinv_d <= 0.375) {
            if (xinv_d <= 0.3125) {
                p = dd_neval_arr(z, j0_P3r2_4N_hi, j0_P3r2_4N_lo) / dd_deval_arr(z, j0_P3r2_4D_hi, j0_P3r2_4D_lo);
                q = dd_neval_arr(z, j0_Q3r2_4N_hi, j0_Q3r2_4N_lo) / dd_deval_arr(z, j0_Q3r2_4D_hi, j0_Q3r2_4D_lo);
            } else {
                p = dd_neval_arr(z, j0_P2r7_3r2N_hi, j0_P2r7_3r2N_lo) / dd_deval_arr(z, j0_P2r7_3r2D_hi, j0_P2r7_3r2D_lo);
                q = dd_neval_arr(z, j0_Q2r7_3r2N_hi, j0_Q2r7_3r2N_lo) / dd_deval_arr(z, j0_Q2r7_3r2D_hi, j0_Q2r7_3r2D_lo);
            }
        } else if (xinv_d <= 0.4375) {
            p = dd_neval_arr(z, j0_P2r3_2r7N_hi, j0_P2r3_2r7N_lo) / dd_deval_arr(z, j0_P2r3_2r7D_hi, j0_P2r3_2r7D_lo);
            q = dd_neval_arr(z, j0_Q2r3_2r7N_hi, j0_Q2r3_2r7N_lo) / dd_deval_arr(z, j0_Q2r3_2r7D_hi, j0_Q2r3_2r7D_lo);
        } else {
            p = dd_neval_arr(z, j0_P2_2r3N_hi, j0_P2_2r3N_lo) / dd_deval_arr(z, j0_P2_2r3D_hi, j0_P2_2r3D_lo);
            q = dd_neval_arr(z, j0_Q2_2r3N_hi, j0_Q2_2r3N_lo) / dd_deval_arr(z, j0_Q2_2r3D_hi, j0_Q2_2r3D_lo);
        }
    }
    p = MFD2(1.0) + z * p;
    q = (z * q - 0.125) * xinv;
    MFD2 tpi = dd_pair(TWO_OVER_PI_hi, TWO_OVER_PI_lo);
    return multifloats::sqrt(tpi / x) * (p * s + q * c);
}

float64x2 y1(float64x2 const &x) {
    if (x._limbs[0] <= 0.0) return MFD2(-INFINITY);
    double xx = x._limbs[0];

    if (xx <= 1e-30) {
        MFD2 tpi = dd_pair(TWO_OVER_PI_hi, TWO_OVER_PI_lo);
        return -tpi / x;
    }

    if (xx <= 2.0) {
        MFD2 z = x * x;
        MFD2 p = x * dd_neval_arr(z, j1_Y1_2N_hi, j1_Y1_2N_lo) /
                     dd_deval_arr(z, j1_Y1_2D_hi, j1_Y1_2D_lo);
        MFD2 tpi = dd_pair(TWO_OVER_PI_hi, TWO_OVER_PI_lo);
        return tpi * (multifloats::log(x) * j1(x) - MFD2(1.0) / x) + p;
    }

    MFD2 angle = x - dd_pair(THREE_PI_QUARTER_hi, THREE_PI_QUARTER_lo);
    MFD2 s = multifloats::sin(angle);
    MFD2 c = multifloats::cos(angle);

    MFD2 xinv = MFD2(1.0) / x;
    MFD2 z = xinv * xinv;
    double xinv_d = 1.0 / xx;
    MFD2 p, q;
    // same p, q as j1 for x > 2
    if (xinv_d <= 0.25) {
        if (xinv_d <= 0.125) {
            if (xinv_d <= 0.0625) {
                p = dd_neval_arr(z, j1_P16_IN_hi, j1_P16_IN_lo) / dd_deval_arr(z, j1_P16_ID_hi, j1_P16_ID_lo);
                q = dd_neval_arr(z, j1_Q16_IN_hi, j1_Q16_IN_lo) / dd_deval_arr(z, j1_Q16_ID_hi, j1_Q16_ID_lo);
            } else {
                p = dd_neval_arr(z, j1_P8_16N_hi, j1_P8_16N_lo) / dd_deval_arr(z, j1_P8_16D_hi, j1_P8_16D_lo);
                q = dd_neval_arr(z, j1_Q8_16N_hi, j1_Q8_16N_lo) / dd_deval_arr(z, j1_Q8_16D_hi, j1_Q8_16D_lo);
            }
        } else if (xinv_d <= 0.1875) {
            p = dd_neval_arr(z, j1_P5_8N_hi, j1_P5_8N_lo) / dd_deval_arr(z, j1_P5_8D_hi, j1_P5_8D_lo);
            q = dd_neval_arr(z, j1_Q5_8N_hi, j1_Q5_8N_lo) / dd_deval_arr(z, j1_Q5_8D_hi, j1_Q5_8D_lo);
        } else {
            p = dd_neval_arr(z, j1_P4_5N_hi, j1_P4_5N_lo) / dd_deval_arr(z, j1_P4_5D_hi, j1_P4_5D_lo);
            q = dd_neval_arr(z, j1_Q4_5N_hi, j1_Q4_5N_lo) / dd_deval_arr(z, j1_Q4_5D_hi, j1_Q4_5D_lo);
        }
    } else {
        if (xinv_d <= 0.375) {
            if (xinv_d <= 0.3125) {
                p = dd_neval_arr(z, j1_P3r2_4N_hi, j1_P3r2_4N_lo) / dd_deval_arr(z, j1_P3r2_4D_hi, j1_P3r2_4D_lo);
                q = dd_neval_arr(z, j1_Q3r2_4N_hi, j1_Q3r2_4N_lo) / dd_deval_arr(z, j1_Q3r2_4D_hi, j1_Q3r2_4D_lo);
            } else {
                p = dd_neval_arr(z, j1_P2r7_3r2N_hi, j1_P2r7_3r2N_lo) / dd_deval_arr(z, j1_P2r7_3r2D_hi, j1_P2r7_3r2D_lo);
                q = dd_neval_arr(z, j1_Q2r7_3r2N_hi, j1_Q2r7_3r2N_lo) / dd_deval_arr(z, j1_Q2r7_3r2D_hi, j1_Q2r7_3r2D_lo);
            }
        } else if (xinv_d <= 0.4375) {
            p = dd_neval_arr(z, j1_P2r3_2r7N_hi, j1_P2r3_2r7N_lo) / dd_deval_arr(z, j1_P2r3_2r7D_hi, j1_P2r3_2r7D_lo);
            q = dd_neval_arr(z, j1_Q2r3_2r7N_hi, j1_Q2r3_2r7N_lo) / dd_deval_arr(z, j1_Q2r3_2r7D_hi, j1_Q2r3_2r7D_lo);
        } else {
            p = dd_neval_arr(z, j1_P2_2r3N_hi, j1_P2_2r3N_lo) / dd_deval_arr(z, j1_P2_2r3D_hi, j1_P2_2r3D_lo);
            q = dd_neval_arr(z, j1_Q2_2r3N_hi, j1_Q2_2r3N_lo) / dd_deval_arr(z, j1_Q2_2r3D_hi, j1_Q2_2r3D_lo);
        }
    }
    p = MFD2(1.0) + z * p;
    q = (z * q + 0.375) * xinv;
    MFD2 tpi = dd_pair(TWO_OVER_PI_hi, TWO_OVER_PI_lo);
    return multifloats::sqrt(tpi / x) * (p * s + q * c);
}

} // namespace bessel_improved