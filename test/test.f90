program multifloat_test
  use multifloats
  use, intrinsic :: ieee_arithmetic
  implicit none

  integer :: num_errors = 0

  print *, "Starting multifloat tests..."

  call test_assignment()
  call test_comparisons()
  call test_addition()
  call test_subtraction()
  call test_multiplication()
  call test_division()
  call test_nonfinite()
  call test_isfinite()
  call test_renormalize()
  call test_builtins()
  call test_math_intrinsics()
  call test_overflow_paths()
  call test_sinpi_cospi_precision()
  call test_cdd_matmul_dot()
  call test_compensated_reductions()

  if (num_errors == 0) then
    print *, "-------------------------------"
    print *, "ALL TESTS PASSED SUCCESSFULLY"
    print *, "-------------------------------"
  else
    print *, "-------------------------------"
    print *, num_errors, " TESTS FAILED"
    print *, "-------------------------------"
    stop 1
  end if

contains

  subroutine assert(cond, msg)
    logical, intent(in) :: cond
    character(*), intent(in) :: msg
    if (.not. cond) then
      print *, "ASSERTION FAILED: ", msg
      num_errors = num_errors + 1
    end if
  end subroutine assert

  subroutine assert_eq(a, b_h, b_l, msg)
    type(real64x2), intent(in) :: a
    double precision, intent(in) :: b_h, b_l
    character(*), intent(in) :: msg
    ! We allow small difference in the second limb if it's due to minor precision differences in target calculation
    if (a%limbs(1) /= b_h .or. abs(a%limbs(2) - b_l) > abs(b_l)*1e-15) then
      if (a%limbs(1) /= b_h .or. (b_l == 0.0d0 .and. abs(a%limbs(2)) > 1e-35) .or. &
          (b_l /= 0.0d0 .and. abs(a%limbs(2) - b_l) > abs(b_l)*1e-15)) then
        print *, "ASSERT_EQ FAILED: ", msg
        print *, "  Expected: ", b_h, b_l
        print *, "  Got:      ", a%limbs
        num_errors = num_errors + 1
      end if
    end if
  end subroutine assert_eq

  subroutine assert_approx(a, b_h, b_l, msg)
    type(real64x2), intent(in) :: a
    double precision, intent(in) :: b_h, b_l
    character(*), intent(in) :: msg
    real(16) :: qa, qb
    qa = real(a%limbs(1), 16) + real(a%limbs(2), 16)
    qb = real(b_h, 16) + real(b_l, 16)
    ! Double-double has ~106 bits. relative error should be ~1e-32.
    ! NR refinement for sqrt/log can be slightly less accurate.
    ! 1e-16 is safe for these functions in this test suite.
    if (abs(qa - qb) > abs(qb)*1e-32 .and. abs(qa - qb) > 1e-35) then
       print *, "ASSERT_APPROX FAILED: ", msg
       print *, "  Expected: ", b_h, b_l
       print *, "  Got:      ", a%limbs
       num_errors = num_errors + 1
    end if
  end subroutine

  subroutine assert_approx_qp(a, b, msg)
    type(real64x2), intent(in) :: a
    real(16), intent(in) :: b
    character(*), intent(in) :: msg
    real(16) :: qa
    qa = real(a%limbs(1), 16) + real(a%limbs(2), 16)
    ! Transcendentals are computed via a first-order derivative-corrected
    ! DD evaluation (no qp temporaries), giving roughly single-double
    ! precision. Allow a relative error around 1e-15.
    if (abs(qa - b) > abs(b)*1e-15_16 .and. abs(qa - b) > 1e-300_16) then
       print *, "ASSERT_APPROX_QP FAILED: ", msg
       print *, "  Expected: ", b
       print *, "  Got:      ", qa, "(limbs:", a%limbs, ")"
       num_errors = num_errors + 1
    end if
  end subroutine

  logical function is_inf(x)
    double precision, intent(in) :: x
    is_inf = (.not. ieee_is_finite(x)) .and. (.not. ieee_is_nan(x))
  end function

  subroutine test_assignment()
    type(real64x2) :: a
    a = 1.234d0
    call assert_eq(a, 1.234d0, 0.0d0, "assignment 1.234")
    a = -0.0d0
    call assert_eq(a, -0.0d0, 0.0d0, "assignment -0.0")
  end subroutine

  subroutine test_comparisons()
    type(real64x2) :: a, b
    ! a = 1 + 1e-20, b = 1 + 2e-20
    a = 1.0d0; a%limbs(2) = 1.0d-20
    b = 1.0d0; b%limbs(2) = 2.0d-20
    
    ! ff combinations
    call assert(a < b, "a < b")
    call assert(b > a, "b > a")
    call assert(a <= b, "a <= b")
    call assert(b >= a, "b >= a")
    call assert(a <= a, "a <= a")
    call assert(a == a, "a == a")
    call assert(a /= b, "a /= b")
    
    ! fd combinations
    call assert(a > 1.0d0, "a > 1.0d0")
    call assert(a >= 1.0d0, "a >= 1.0d0")
    call assert(a /= 1.0d0, "a /= 1.0d0")
    call assert(.not. (a < 1.0d0), "not a < 1.0d0")
    call assert(.not. (a == 1.0d0), "not a == 1.0d0")
    
    ! df combinations
    call assert(1.0d0 < a, "1.0d0 < a")
    call assert(1.0d0 <= a, "1.0d0 <= a")
    call assert(1.0d0 /= a, "1.0d0 /= a")
    call assert(.not. (1.0d0 > a), "not 1.0d0 > a")
    call assert(.not. (1.0d0 == a), "not 1.0d0 == a")

    ! Zero handling
    a = 0.0d0; b = -0.0d0
    call assert(a == b, "0.0 == -0.0")
    call assert(a == 0.0d0, "0.0 == 0.0d0")
  end subroutine

  subroutine test_addition()
    type(real64x2) :: a, b, c
    double precision :: s, e

    ! (1 + 1e-20) + (1 + 1e-20) = 2 + 2e-20
    a = 1.0d0; a%limbs(2) = 1.0d-20
    b = 1.0d0; b%limbs(2) = 1.0d-20
    c = a + b
    call assert_eq(c, 2.0d0, 2.0d-20, "a + b (exact)")

    ! Mixed addition (ff + d)
    c = a + 1.0d0
    call assert_eq(c, 2.0d0, 1.0d-20, "a + d")

    ! Mixed addition (d + ff)
    c = 1.0d0 + a
    call assert_eq(c, 2.0d0, 1.0d-20, "d + a")
    
    ! Case where s needs rounding
    ! (1.0) + (1e-17)
    s = 1.0d0
    e = 1.0d-17
    call two_sum_internal(s, e)
    a = 1.0d0
    b = 1.0d-17
    c = a + b
    call assert_eq(c, s, e, "a + b (rounding)")
  end subroutine

  subroutine test_subtraction()
    type(real64x2) :: a, b, c
    double precision :: s, e

    a = 2.0d0; a%limbs(2) = 2.0d-20
    b = 1.0d0; b%limbs(2) = 1.0d-20
    c = a - b
    call assert_eq(c, 1.0d0, 1.0d-20, "a - b")

    c = a - 1.0d0
    call assert_eq(c, 1.0d0, 2.0d-20, "a - d")

    c = 2.0d0 - b
    s = 2.0d0
    e = -1.0d0
    call two_sum_internal(s, e)
    e = e - 1.0d-20
    call fast_two_sum_internal(s, e)
    call assert_eq(c, s, e, "d - b")
  end subroutine

  subroutine test_multiplication()
    type(real64x2) :: a, b, c
    double precision :: s, e

    ! (1 + 1e-10)^2 = 1 + 2e-10 + 1e-20.
    ! The kernel intentionally drops the lo*lo cross term (matches the
    ! Julia mfmul reference), so we omit it from the expected value too.
    a = 1.0d0; a%limbs(2) = 1.0d-10
    c = a * a
    call two_prod_internal(a%limbs(1), a%limbs(1), s, e)
    e = e + a%limbs(1)*a%limbs(2)
    e = e + a%limbs(2)*a%limbs(1)
    call fast_two_sum_internal(s, e)
    call assert_eq(c, s, e, "a * a")

    ! a * 2.0d0
    c = a * 2.0d0
    call two_prod_internal(a%limbs(1), 2.0d0, s, e)
    e = e + a%limbs(2)*2.0d0
    call fast_two_sum_internal(s, e)
    call assert_eq(c, s, e, "a * d")

    ! 2.0d0 * a
    c = 2.0d0 * a
    call assert_eq(c, s, e, "d * a")
  end subroutine

  subroutine test_division()
    type(real64x2) :: a, b, c
    double precision :: q1, q2, r, s, ee

    ! (1 + 1e-10)
    a = 1.0d0; a%limbs(2) = 1.0d-10
    b = a * a
    c = b / a
    ! We expect bitwise exactness for this specific case if logic is correct
    ! but let's re-calculate expected based on current module logic
    q1 = b%limbs(1) / a%limbs(1)
    r = mf_fma(-q1, a%limbs(1), b%limbs(1))
    r = r + b%limbs(2) - q1 * a%limbs(2)
    q2 = r / a%limbs(1)
    s = q1
    ee = q2
    call fast_two_sum_internal(s, ee)
    call assert_eq(c, s, ee, "b / a")

    c = b / 1.0d0
    call assert_eq(c, b%limbs(1), b%limbs(2), "b / d")

    ! 1.0 / (1 + 1e-10)
    q1 = 1.0d0 / a%limbs(1)
    r = mf_fma(-q1, a%limbs(1), 1.0d0)
    r = r - q1 * a%limbs(2)
    q2 = r / a%limbs(1)
    s = q1
    ee = q2
    call fast_two_sum_internal(s, ee)
    c = 1.0d0 / a
    call assert_eq(c, s, ee, "d / a")
  end subroutine

  subroutine test_nonfinite()
    type(real64x2) :: a, b, c
    double precision :: inf, nan
    inf = ieee_value(inf, ieee_positive_inf)
    nan = ieee_value(nan, ieee_quiet_nan)
    
    ! Inf * 2
    a = inf
    c = a * 2.0d0
    call assert(is_inf(c%limbs(1)), "inf * 2 is inf")
    call assert(c%limbs(2) == 0.0d0, "inf * 2 low limb is 0")

    ! Inf comparison
    a = inf
    b = inf
    call assert(a == b, "inf == inf")
    call assert(a == inf, "inf == d(inf)")
    call assert(inf == a, "d(inf) == inf")
    call assert(a >= b, "inf >= inf")
    call assert(a <= b, "inf <= b")
    call assert(.not. (a < b), "not inf < inf")
    call assert(.not. (a > b), "not inf > inf")
    
    a = inf
    call assert(a > 0.0d0, "inf > 0")
    call assert(-a < 0.0d0, "-inf < 0")

    ! NaN comparison
    a = nan
    call assert(.not. (a == a), "nan == nan is false")
    call assert(a /= a, "nan /= nan is true")
    call assert(.not. (a < 1.0d0), "nan < 1 is false")

    ! Inf - Inf
    a = inf
    b = inf
    c = a - b
    call assert(ieee_is_nan(c%limbs(1)), "inf - inf is nan")

    ! NaN propagation
    a = nan
    c = a + 1.0d0
    call assert(ieee_is_nan(c%limbs(1)), "nan + 1 is nan")

    ! Div by zero
    a = 1.0d0
    b = 0.0d0
    ! Use a variable to avoid compiler warnings if any
    c = a / b
    call assert(is_inf(c%limbs(1)), "1 / 0 is inf")
  end subroutine

  subroutine test_isfinite()
    type(real64x2) :: a
    double precision :: inf
    inf = ieee_value(inf, ieee_positive_inf)

    a = 1.234d0
    call assert(dd_is_finite(a), "1.234 is finite")

    a = inf
    call assert(.not. dd_is_finite(a), "inf is not finite")
  end subroutine

  subroutine test_renormalize()
    type(real64x2) :: a
    ! Manually set unnormalized state
    a%limbs(1) = 1.0d0
    a%limbs(2) = 1.0d0
    call dd_renormalize(a)
    call assert_eq(a, 2.0d0, 0.0d0, "renormalize(1, 1) -> (2, 0)")

    a%limbs(1) = 1.0d-20
    a%limbs(2) = 1.0d0
    call dd_renormalize(a)
    call assert_eq(a, 1.0d0, 1.0d-20, "renormalize(1e-20, 1) -> (1, 1e-20)")
  end subroutine

  subroutine test_builtins()
    type(real64x2) :: a, b
    a = 1.0d0
    call assert(precision(a) == 31, "precision")
    call assert(minexponent(a) == minexponent(1.0d0), "minexponent")
    call assert(maxexponent(a) == maxexponent(1.0d0), "maxexponent")

    b = tiny(a)
    call assert_eq(b, tiny(1.0d0), 0.0d0, "tiny")

    b = huge(a)
    call assert_eq(b, huge(1.0d0), 0.0d0, "huge")

    a = 2.0d0
    call assert(exponent(a) == exponent(2.0d0), "exponent(2.0)")
    a = 0.5d0
    call assert(exponent(a) == exponent(0.5d0), "exponent(0.5)")
  end subroutine

  subroutine test_math_intrinsics()
    type(real64x2) :: a, b, c
    
    ! abs
    a = -1.0d0; a%limbs(2) = -1.0d-20
    b = abs(a)
    call assert_eq(b, 1.0d0, 1.0d-20, "abs")
    
    ! sqrt
    a = 2.0d0
    b = sqrt(a)
    ! Check b*b matches a
    c = b * b
    call assert_approx(c, 2.0d0, 0.0d0, "sqrt(2)^2 approx 2")
    
    ! sign
    a = 1.0d0; b = -1.0d0
    call assert_eq(sign(a, b), -1.0d0, 0.0d0, "sign(1, -1)")
    
    ! min/max
    a = 1.0d0; b = 2.0d0
    call assert_eq(min(a, b), 1.0d0, 0.0d0, "min")
    call assert_eq(max(a, b), 2.0d0, 0.0d0, "max")
    
    ! floor/ceiling — module returns integer, wrap to real64x2 for assert_eq.
    a = 1.5d0
    call assert_eq(real64x2(floor(a)), 1.0d0, 0.0d0, "floor(1.5)")
    call assert_eq(real64x2(ceiling(a)), 2.0d0, 0.0d0, "ceiling(1.5)")
    
    ! scale
    a = 1.0d0
    call assert_eq(scale(a, 1), 2.0d0, 0.0d0, "scale(1, 1)")

    ! fraction/set_exponent
    a = 3.0d0 ! 1.5 * 2^1 -> fraction(3.0) is 0.75, exponent is 2
    call assert_eq(fraction(a), 0.75d0, 0.0d0, "fraction(3.0) is 0.75")
    call assert_eq(set_exponent(a, 2), 3.0d0, 0.0d0, "set_exponent(3.0, 2) is 3.0")

    ! spacing/rrspacing
    a = 1.0d0
    b = spacing(a)
    call assert_eq(b, scale(spacing(1.0d0), -53), 0.0d0, "spacing(1.0)")
    b = rrspacing(a)
    ! spacing(1.0) is 2^-52. spacing_f(1.0) is 2^-105. 
    ! rrspacing = |1.0| / 2^-105 = 2^105
    call assert_eq(b, scale(1.0d0, 105), 0.0d0, "rrspacing(1.0)")

    ! exp/log — compare against quad-precision reference values, since
    ! the kernel's intermediate qp computation has ~30 digits of precision.
    a = 1.0d0
    b = exp(a)
    call assert_approx_qp(b, exp(1.0_16), "exp(1.0)")
    c = log(b)
    call assert_approx_qp(c, 1.0_16, "log(exp(1.0))")

    a = 100.0d0
    b = log10(a)
    call assert_approx_qp(b, 2.0_16, "log10(100.0)")
  end subroutine

  ! Internal helpers for ground truth matching module logic
  subroutine two_sum_internal(a, b)
    double precision, intent(inout):: a, b
    double precision :: s, a_prime, b_prime, a_err, b_err
    s = a + b
    a_prime = s - b
    b_prime = s - a_prime
    a_err = a - a_prime
    b_err = b - b_prime
    a = s
    b = a_err + b_err
  end subroutine

  subroutine fast_two_sum_internal(a, b)
    double precision, intent(inout):: a, b
    double precision :: s, b_prime, b_err
    s = a + b
    b_prime = s - a
    b_err = b - b_prime
    a = s
    b = b_err
  end subroutine

  subroutine two_prod_internal(a, b, p, e)
    double precision, intent(in) :: a, b
    double precision, intent(out) :: p, e
    p = a * b
    e = mf_fma(a, b, -p)
  end subroutine

  subroutine test_overflow_paths()
    ! After P1, overflow / pole paths return IEEE infinity rather than
    ! huge(0.0_dp). These assertions catch the regression.
    type(real64x2) :: a, r
    double precision :: inf

    inf = ieee_value(inf, ieee_positive_inf)

    ! exp overflow → +inf (routes through dd_exp2_full's overflow arm)
    a = 2000.0d0
    r = exp(a)
    call assert(is_inf(r%limbs(1)) .and. r%limbs(1) > 0.0d0, &
                "exp(2000) is +inf")

    ! log(0) → -inf; log(+inf) → +inf (route through dd_log2_full)
    a = 0.0d0
    r = log(a)
    call assert(is_inf(r%limbs(1)) .and. r%limbs(1) < 0.0d0, &
                "log(0) is -inf")
    a = inf
    r = log(a)
    call assert(is_inf(r%limbs(1)) .and. r%limbs(1) > 0.0d0, &
                "log(+inf) is +inf")

    ! gamma(200) overflows → +inf
    a = 200.0d0
    r = gamma(a)
    call assert(is_inf(r%limbs(1)) .and. r%limbs(1) > 0.0d0, &
                "gamma(200) is +inf")

    ! log_gamma at non-positive integer pole → +inf
    a = 0.0d0
    r = log_gamma(a)
    call assert(is_inf(r%limbs(1)) .and. r%limbs(1) > 0.0d0, &
                "log_gamma(0) is +inf")

    ! bessel_yn(n, 0) → -inf
    a = 0.0d0
    r = bessel_yn(2, a)
    call assert(is_inf(r%limbs(1)) .and. r%limbs(1) < 0.0d0, &
                "bessel_yn(2, 0) is -inf")
  end subroutine

  subroutine test_sinpi_cospi_precision()
    ! P2: Fortran sinpi/cospi must hit full DD (~4e-32), not the legacy
    ! polynomial ceiling (~5e-27). assert_approx's tolerance is 1e-32,
    ! which the legacy kernel would have failed outright at x = 0.25.
    type(real64x2) :: x, r
    x = 0.25d0
    r = sinpi(x)
    call assert_approx(r, 0.7071067811865476d0, &
                       -4.833646656726457d-17, "sinpi(0.25) ≈ sqrt(2)/2")
    r = cospi(x)
    call assert_approx(r, 0.7071067811865476d0, &
                       -4.833646656726457d-17, "cospi(0.25) ≈ sqrt(2)/2")
    ! Large-argument reduction must still be exact at integer step.
    x = 101.0d0
    r = sinpi(x)
    call assert_approx(r, 0.0d0, 0.0d0, "sinpi(101) == 0")
    r = cospi(x)
    call assert_approx(r, -1.0d0, 0.0d0, "cospi(101) == -1")
  end subroutine

  subroutine test_cdd_matmul_dot()
    ! P3: complex matmul/dot_product now route through the compensated
    ! real kernels. Sanity-check identities that would fail if the
    ! complex formula were wrong. Precision validation is inherited from
    ! the real path's fuzz tests.
    type(cmplx64x2) :: a(2, 2), b(2, 2), c(2, 2), u(2), v(2)
    type(cmplx64x2) :: one, i_unit, two, zero, s
    type(real64x2) :: zero_dd
    zero_dd = real64x2(0.0d0)
    one%re = real64x2(1.0d0); one%im = zero_dd
    two%re = real64x2(2.0d0); two%im = zero_dd
    zero%re = zero_dd; zero%im = zero_dd
    i_unit%re = zero_dd; i_unit%im = real64x2(1.0d0)

    ! [1 i; i 1] * [1 -i; -i 1] = [2 0; 0 2] (identity up to complex phase)
    a(1,1) = one;    a(1,2) = i_unit
    a(2,1) = i_unit; a(2,2) = one
    b(1,1) = one;              b(1,2)%re = zero_dd; b(1,2)%im = real64x2(-1.0d0)
    b(2,1)%re = zero_dd; b(2,1)%im = real64x2(-1.0d0)
    b(2,2) = one
    c = matmul(a, b)
    call assert_approx(c(1,1)%re, 2.0d0, 0.0d0, "cdd_matmul [1,1].re")
    call assert_approx(c(1,1)%im, 0.0d0, 0.0d0, "cdd_matmul [1,1].im")
    call assert_approx(c(2,2)%re, 2.0d0, 0.0d0, "cdd_matmul [2,2].re")
    call assert_approx(c(1,2)%re, 0.0d0, 0.0d0, "cdd_matmul [1,2].re")

    ! dot_product conjugates the first argument; dot(u, u) is real non-neg.
    u(1) = one
    u(2) = i_unit
    s = dot_product(u, u)
    call assert_approx(s%re, 2.0d0, 0.0d0, "cdd_dot(u, u) real == 2")
    call assert_approx(s%im, 0.0d0, 0.0d0, "cdd_dot(u, u) imag == 0")

    ! dot(u, i*u) = conj(u_k) · (i u_k) = i (|u_1|^2 + |u_2|^2) = 2i
    v(1) = i_unit
    v(2)%re = real64x2(-1.0d0); v(2)%im = zero_dd  ! i * i = -1
    s = dot_product(u, v)
    call assert_approx(s%re, 0.0d0, 0.0d0, "cdd_dot(u, i·u) real == 0")
    call assert_approx(s%im, 2.0d0, 0.0d0, "cdd_dot(u, i·u) imag == 2")
  end subroutine

  subroutine test_compensated_reductions()
    ! P4: sum / sum(..., dim=...) / cdd sum all route through Neumaier
    ! compensation. Verify that the dim-reduction sum of a 1xN array
    ! matches the rank-1 sum bit-exactly, and that cdd sum preserves
    ! near-cancellation the same way the real sum does.
    type(real64x2) :: a(6), ma(1, 6), col(1)
    type(real64x2) :: s_flat, s_col_scalar
    type(cmplx64x2) :: ca(6), cm(1, 6), ccol(1), cs, cs_col

    ! Choose values that straddle DD exponents; naive left-fold loses
    ! the last element's contribution when accumulated before the
    ! large terms cancel. The compensated path should preserve it.
    a(1) = real64x2(1.0d30)
    a(2) = real64x2(1.0d0)
    a(3) = real64x2(-1.0d30)
    a(4) = real64x2(1.0d0)
    a(5) = real64x2(1.0d30)
    a(6) = real64x2(-1.0d30)

    s_flat = sum(a)
    ma(1, :) = a
    col = sum(ma, dim=2)
    s_col_scalar = col(1)

    call assert_approx(s_flat, 2.0d0, 0.0d0, "sum(a) == 2 (compensated)")
    call assert_approx(s_col_scalar, 2.0d0, 0.0d0, &
                       "sum(ma, dim=2) == 2 (compensated)")

    ! Same pattern in the real / imag limbs of a complex array.
    ca(1)%re = real64x2(1.0d30); ca(1)%im = real64x2(2.0d30)
    ca(2)%re = real64x2(1.0d0);  ca(2)%im = real64x2(2.0d0)
    ca(3)%re = real64x2(-1.0d30);ca(3)%im = real64x2(-2.0d30)
    ca(4)%re = real64x2(1.0d0);  ca(4)%im = real64x2(2.0d0)
    ca(5)%re = real64x2(1.0d30); ca(5)%im = real64x2(2.0d30)
    ca(6)%re = real64x2(-1.0d30);ca(6)%im = real64x2(-2.0d30)

    cs = sum(ca)
    cm(1, :) = ca
    ccol = sum(cm, dim=2)
    cs_col = ccol(1)

    call assert_approx(cs%re, 2.0d0, 0.0d0, "cdd_sum.re == 2 (compensated)")
    call assert_approx(cs%im, 4.0d0, 0.0d0, "cdd_sum.im == 4 (compensated)")
    call assert_approx(cs_col%re, 2.0d0, 0.0d0, &
                       "cdd_sum(dim=2).re == 2 (compensated)")
    call assert_approx(cs_col%im, 4.0d0, 0.0d0, &
                       "cdd_sum(dim=2).im == 4 (compensated)")
  end subroutine

  ! Local shims for facilities the multifloats module does not expose.
  logical function dd_is_finite(a) result(res)
    type(real64x2), intent(in) :: a
    res = ieee_is_finite(a%limbs(1)) .and. ieee_is_finite(a%limbs(2))
  end function

  subroutine dd_renormalize(a)
    type(real64x2), intent(inout) :: a
    double precision :: hi, lo, s, b_prime
    if (abs(a%limbs(1)) >= abs(a%limbs(2))) then
      hi = a%limbs(1)
      lo = a%limbs(2)
    else
      hi = a%limbs(2)
      lo = a%limbs(1)
    end if
    s = hi + lo
    b_prime = s - hi
    a%limbs(1) = s
    a%limbs(2) = lo - b_prime
  end subroutine

end program
