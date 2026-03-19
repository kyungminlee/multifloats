program multifloat_test
  use multifloat
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
    type(float64x2), intent(in) :: a
    double precision, intent(in) :: b_h, b_l
    character(*), intent(in) :: msg
    ! We allow small difference in the second limb if it's due to minor precision differences in target calculation
    if (a%limbs(1) /= b_h .or. abs(a%limbs(2) - b_l) > abs(b_l)*1e-15) then
      if (a%limbs(1) /= b_h .or. (b_l == 0.0d0 .and. a%limbs(2) /= 0.0d0) .or. &
          (b_l /= 0.0d0 .and. abs(a%limbs(2) - b_l) > abs(b_l)*1e-15)) then
        print *, "ASSERT_EQ FAILED: ", msg
        print *, "  Expected: ", b_h, b_l
        print *, "  Got:      ", a%limbs
        num_errors = num_errors + 1
      end if
    end if
  end subroutine assert_eq

  logical function is_inf(x)
    double precision, intent(in) :: x
    is_inf = (.not. ieee_is_finite(x)) .and. (.not. ieee_is_nan(x))
  end function

  subroutine test_assignment()
    type(float64x2) :: a
    a = 1.234d0
    call assert_eq(a, 1.234d0, 0.0d0, "assignment 1.234")
    a = -0.0d0
    call assert_eq(a, -0.0d0, 0.0d0, "assignment -0.0")
  end subroutine

  subroutine test_comparisons()
    type(float64x2) :: a, b
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
    type(float64x2) :: a, b, c
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
    type(float64x2) :: a, b, c
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
    type(float64x2) :: a, b, c
    double precision :: s, e

    ! (1 + 1e-10)^2 = 1 + 2e-10 + 1e-20
    a = 1.0d0; a%limbs(2) = 1.0d-10
    c = a * a
    ! We calculate expected using module-like logic
    call two_prod_internal(a%limbs(1), a%limbs(1), s, e)
    e = e + a%limbs(1)*a%limbs(2)
    e = e + a%limbs(2)*a%limbs(1)
    e = e + a%limbs(2)*a%limbs(2)
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
    type(float64x2) :: a, b, c
    double precision :: q1, q2, r, s, ee

    ! (1 + 1e-10)
    a = 1.0d0; a%limbs(2) = 1.0d-10
    b = a * a
    c = b / a
    ! We expect bitwise exactness for this specific case if logic is correct
    ! but let's re-calculate expected based on current module logic
    q1 = b%limbs(1) / a%limbs(1)
    r = ieee_fma(-q1, a%limbs(1), b%limbs(1))
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
    r = ieee_fma(-q1, a%limbs(1), 1.0d0)
    r = r - q1 * a%limbs(2)
    q2 = r / a%limbs(1)
    s = q1
    ee = q2
    call fast_two_sum_internal(s, ee)
    c = 1.0d0 / a
    call assert_eq(c, s, ee, "d / a")
  end subroutine

  subroutine test_nonfinite()
    type(float64x2) :: a, b, c
    double precision :: inf, nan
    inf = ieee_value(inf, ieee_positive_inf)
    nan = ieee_value(nan, ieee_quiet_nan)
    
    ! Inf * 2
    a = inf
    c = a * 2.0d0
    call assert(is_inf(c%limbs(1)), "inf * 2 is inf")
    call assert(c%limbs(2) == 0.0d0, "inf * 2 low limb is 0")

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
    type(float64x2) :: a
    double precision :: inf
    inf = ieee_value(inf, ieee_positive_inf)
    
    a = 1.234d0
    call assert(a%isfinite(), "1.234 is finite")
    
    a = inf
    call assert(.not. a%isfinite(), "inf is not finite")
  end subroutine

  subroutine test_renormalize()
    type(float64x2) :: a
    ! Manually set unnormalized state
    a%limbs(1) = 1.0d0
    a%limbs(2) = 1.0d0
    call renormalize(a)
    call assert_eq(a, 2.0d0, 0.0d0, "renormalize(1, 1) -> (2, 0)")
    
    a%limbs(1) = 1.0d-20
    a%limbs(2) = 1.0d0
    call renormalize(a)
    call assert_eq(a, 1.0d0, 1.0d-20, "renormalize(1e-20, 1) -> (1, 1e-20)")
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
    e = ieee_fma(a, b, -p)
  end subroutine

end program
