! bind(c) wrappers over the NATIVE float64x2 elemental operators. Used
! by the cross-language differential fuzz driver (test/crosscheck.cc) to
! exercise both the C++ and Fortran implementations of the DD kernels
! that are INDEPENDENTLY implemented on each side — add, sub, mul, div,
! neg, abs, and the ordered/equality comparisons. A bit-for-bit
! divergence on both limbs flags real algorithmic drift between
! src/multifloats_math.cc and fsrc/multifloats.fypp.
!
! Distinct from dd_bindc.f90 in intent: dd_bindc reimplements the DD
! algorithm in pure Fortran to isolate ABI overhead in bench_abi, while
! fnat_* delegates straight to `a + b` / sqrt(a) / etc. so the
! differential compares two production code paths rather than a
! hand-mirrored reimplementation.
module crosscheck_bindings
  ! Import the whole module so the overloaded arithmetic / comparison
  ! operators and the `sqrt` / `abs` generics are visible. An `only:`
  ! list would drop the operator interfaces and the compiler would
  ! reject `type(float64x2) + type(float64x2)`.
  use multifloats
  use, intrinsic :: iso_c_binding, only: c_double, c_int, c_long
  implicit none
  private

  ! Layout matches float64x2 in include/multifloats.h and dd_c in
  ! dd_bindc.f90. Declared locally so this module does not depend on
  ! dd_bindc (separate concern).
  type, bind(c), public :: dd_c
    real(c_double) :: limbs(2)
  end type dd_c

  public :: fnat_add, fnat_sub, fnat_mul, fnat_div
  public :: fnat_neg, fnat_abs, fnat_sqrt
  public :: fnat_eq,  fnat_ne,  fnat_lt
  ! fnat_le / fnat_ge intentionally omitted: the C++ `<=` is defined as
  ! `!(r < *this)` which is true on NaN pairs (IEEE `<` is false, negation
  ! flips it), whereas Fortran `<=` is IEEE-strict (false on NaN). Any
  ! cross-check would fire on every NaN-bearing iteration for a
  ! pre-existing semantic difference that is out of scope for PR-1.
  !
  ! Integer-returning rounding / decomposition ops — common routines
  ! implemented natively on BOTH sides (not routed through C_DELEGATE_*):
  !
  !   fnat_floor    ↔  C++ (long) mf::floor(x).hi
  !   fnat_ceiling  ↔  C++ (long) mf::ceil(x).hi
  !   fnat_nint     ↔  C++ mf::lround(x)
  !   fnat_exponent ↔  C++ mf::ilogb(x) + 1   (Fortran convention: fraction
  !                                            in [0.5, 1), so exp is 1
  !                                            greater than the C ilogb.)
  public :: fnat_floor, fnat_ceiling, fnat_nint, fnat_exponent

  ! DD-returning scaling ops (Fortran-native, parallel to C++ mf::ldexp /
  ! mf::frexp+mf::ldexp composition):
  !
  !   fnat_scale(x, i)        ↔  C++ mf::ldexp(x, i)
  !   fnat_set_exponent(x, i) ↔  C++ mf::ldexp(mf::frexp(x, &e), i)
  !                              (both compute  fraction(x) * 2^i  via
  !                              exponent-extract + rescale on the hi limb)
  public :: fnat_scale, fnat_set_exponent

contains

  pure function c_to_f(c) result(r)
    type(dd_c), intent(in) :: c
    type(float64x2) :: r
    r%limbs = c%limbs
  end function

  pure function f_to_c(r) result(c)
    type(float64x2), intent(in) :: r
    type(dd_c) :: c
    c%limbs = r%limbs
  end function

  pure function fnat_add(a, b) result(res) bind(c, name='fnat_add')
    type(dd_c), intent(in), value :: a, b
    type(dd_c) :: res
    res = f_to_c(c_to_f(a) + c_to_f(b))
  end function

  pure function fnat_sub(a, b) result(res) bind(c, name='fnat_sub')
    type(dd_c), intent(in), value :: a, b
    type(dd_c) :: res
    res = f_to_c(c_to_f(a) - c_to_f(b))
  end function

  pure function fnat_mul(a, b) result(res) bind(c, name='fnat_mul')
    type(dd_c), intent(in), value :: a, b
    type(dd_c) :: res
    res = f_to_c(c_to_f(a) * c_to_f(b))
  end function

  pure function fnat_div(a, b) result(res) bind(c, name='fnat_div')
    type(dd_c), intent(in), value :: a, b
    type(dd_c) :: res
    res = f_to_c(c_to_f(a) / c_to_f(b))
  end function

  pure function fnat_neg(a) result(res) bind(c, name='fnat_neg')
    type(dd_c), intent(in), value :: a
    type(dd_c) :: res
    res = f_to_c(-c_to_f(a))
  end function

  pure function fnat_abs(a) result(res) bind(c, name='fnat_abs')
    type(dd_c), intent(in), value :: a
    type(dd_c) :: res
    res = f_to_c(abs(c_to_f(a)))
  end function

  ! sqrt is included as an ABI smoke test: on the Fortran side sqrt is
  ! delegated to the C++ sqrtdd kernel (see C_DELEGATE_UNARY_MAP in
  ! fsrc/multifloats.fypp), so any mismatch vs. C++ mf::sqrt indicates
  ! ABI corruption rather than algorithmic drift.
  pure function fnat_sqrt(a) result(res) bind(c, name='fnat_sqrt')
    type(dd_c), intent(in), value :: a
    type(dd_c) :: res
    res = f_to_c(sqrt(c_to_f(a)))
  end function

  pure function fnat_eq(a, b) result(res) bind(c, name='fnat_eq')
    type(dd_c), intent(in), value :: a, b
    integer(c_int) :: res
    res = merge(1, 0, c_to_f(a) == c_to_f(b))
  end function

  pure function fnat_ne(a, b) result(res) bind(c, name='fnat_ne')
    type(dd_c), intent(in), value :: a, b
    integer(c_int) :: res
    res = merge(1, 0, c_to_f(a) /= c_to_f(b))
  end function

  pure function fnat_lt(a, b) result(res) bind(c, name='fnat_lt')
    type(dd_c), intent(in), value :: a, b
    integer(c_int) :: res
    res = merge(1, 0, c_to_f(a) < c_to_f(b))
  end function

  ! floor / ceiling / nint return default integer in the DD-overloaded
  ! interfaces; widen to c_long via explicit cast so the C side can use a
  ! uniform `long` return type. Caller is expected to gate on |x| < 2e9
  ! so the default-integer computation doesn't overflow int32.
  pure function fnat_floor(a) result(res) bind(c, name='fnat_floor')
    type(dd_c), intent(in), value :: a
    integer(c_long) :: res
    res = int(floor(c_to_f(a)), c_long)
  end function

  pure function fnat_ceiling(a) result(res) bind(c, name='fnat_ceiling')
    type(dd_c), intent(in), value :: a
    integer(c_long) :: res
    res = int(ceiling(c_to_f(a)), c_long)
  end function

  pure function fnat_nint(a) result(res) bind(c, name='fnat_nint')
    type(dd_c), intent(in), value :: a
    integer(c_long) :: res
    res = int(nint(c_to_f(a)), c_long)
  end function

  pure function fnat_exponent(a) result(res) bind(c, name='fnat_exponent')
    type(dd_c), intent(in), value :: a
    integer(c_int) :: res
    res = exponent(c_to_f(a))
  end function

  ! scale / set_exponent take a second integer argument. Use c_int on the
  ! boundary so the C++ side can pass a plain `int`. Inside Fortran both
  ! intrinsics accept any integer kind.
  pure function fnat_scale(a, n) result(res) bind(c, name='fnat_scale')
    type(dd_c), intent(in), value :: a
    integer(c_int), intent(in), value :: n
    type(dd_c) :: res
    res = f_to_c(scale(c_to_f(a), n))
  end function

  pure function fnat_set_exponent(a, n) result(res) bind(c, name='fnat_set_exponent')
    type(dd_c), intent(in), value :: a
    integer(c_int), intent(in), value :: n
    type(dd_c) :: res
    res = f_to_c(set_exponent(c_to_f(a), n))
  end function

end module crosscheck_bindings
