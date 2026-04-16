#pragma once
#include "multifloats.hh"

namespace bessel_improved {

using float64x2 = multifloats::MultiFloat<double, 2>;

float64x2 j0(float64x2 const &x);
float64x2 j1(float64x2 const &x);
float64x2 y0(float64x2 const &x);
float64x2 y1(float64x2 const &x);

} // namespace bessel_improved
