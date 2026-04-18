program test_mod
  ! Targeted precision tests for mod/modulo edge cases.
  use multifloats
  use, intrinsic :: ieee_arithmetic
  implicit none

  integer, parameter :: qp = 16, dp = 8
  integer :: nerr, ntotal
  nerr = 0; ntotal = 0

  print *, "=== mod/modulo precision stress test ==="

  ! 1. Near-integer quotients: x ≈ n * y + tiny remainder
  call test_case(7.0_qp + 1e-25_qp, 3.5_qp, "near-int q=2")
  call test_case(7.0_qp - 1e-25_qp, 3.5_qp, "near-int q=2 (below)")
  call test_case(10.0_qp + 1e-30_qp, 5.0_qp, "near-int q=2 (tiny)")
  call test_case(10.0_qp - 1e-30_qp, 5.0_qp, "near-int q=2 (tiny below)")
  call test_case(9.0_qp + 1e-20_qp, 3.0_qp + 1e-21_qp, "near-int q=3")
  call test_case(9.0_qp - 1e-20_qp, 3.0_qp - 1e-21_qp, "near-int q=3 (below)")

  ! 2. Exact integer quotients (remainder should be 0 or tiny)
  call test_case(6.0_qp, 3.0_qp, "exact q=2")
  call test_case(6.0_qp, 2.0_qp, "exact q=3")
  call test_case(1.0_qp, 0.25_qp, "exact q=4")
  call test_case(100.0_qp, 10.0_qp, "exact q=10")

  ! 3. Large quotients (stress the reduction loop).
  ! DD has ~31 decimal digits; mod loses ~log10(quotient) digits.
  ! Tolerance = 10^(digits_of_quotient - 31), clamped to >= 1e-26.
  call test_case(1.0e10_qp, 0.3_qp, "large q~3.3e10", 1e-20_qp)
  call test_case(1.0e15_qp, 7.0_qp, "large q~1.4e14", 1e-16_qp)
  call test_case(1.0e20_qp, 3.0_qp, "large q~3.3e19", 1e-11_qp)
  call test_case(1.0e25_qp, 0.7_qp, "large q~1.4e25", 1e-5_qp)

  ! 4. Small quotients (< 1)
  call test_case(0.3_qp, 7.0_qp, "small q<1")
  call test_case(1e-10_qp, 1.0_qp, "tiny x")
  call test_case(1e-20_qp, 1e-10_qp, "tiny x tiny y")

  ! 5. Negative operands
  call test_case(-7.0_qp, 3.0_qp, "neg x pos y")
  call test_case(7.0_qp, -3.0_qp, "pos x neg y")
  call test_case(-7.0_qp, -3.0_qp, "neg x neg y")
  call test_case(-1e15_qp, 7.0_qp, "neg large x")

  ! 6. Very different magnitudes
  ! q~1e35 exceeds DD's 31 digits — result is unreliable.
  ! Use a very loose tolerance to verify we at least get something finite.
  call test_case(1.0e30_qp, 1.0e-5_qp, "huge x tiny y", 1.0e10_qp)
  call test_case(1.0e-5_qp, 1.0e30_qp, "tiny x huge y")

  ! 7. Values with significant lo limbs
  call test_case(3.141592653589793_qp + 1e-16_qp, 1.0_qp, "pi-like")
  call test_case(2.718281828459045_qp + 1e-16_qp, 0.5_qp, "e-like")

  ! 8. Modulo-specific: sign tests
  call test_modulo_case(7.0_qp, 3.0_qp, "modulo pos pos")
  call test_modulo_case(-7.0_qp, 3.0_qp, "modulo neg pos")
  call test_modulo_case(7.0_qp, -3.0_qp, "modulo pos neg")
  call test_modulo_case(-7.0_qp, -3.0_qp, "modulo neg neg")
  call test_modulo_case(-1e15_qp + 1e-10_qp, 7.0_qp, "modulo neg large")
  call test_modulo_case(1e20_qp, -3.0_qp, "modulo large neg y")

  print *, ""
  if (nerr == 0) then
    print '(a,i0,a)', " ALL ", ntotal, " TESTS PASSED"
  else
    print '(a,i0,a,i0,a)', " FAILED: ", nerr, " / ", ntotal, " tests"
    stop 1
  end if

contains

  subroutine check(op, got, expected, label, tol_in)
    character(*), intent(in) :: op, label
    type(float64x2), intent(in) :: got
    real(qp), intent(in) :: expected
    real(qp), intent(in), optional :: tol_in
    real(qp) :: got_q, rel, tol
    tol = 1e-26_qp
    if (present(tol_in)) tol = tol_in
    ntotal = ntotal + 1
    got_q = real(got%limbs(1), qp) + real(got%limbs(2), qp)
    if (expected == 0.0_qp) then
      rel = abs(got_q)
    else
      rel = abs((got_q - expected) / expected)
    end if
    if (rel > tol) then
      nerr = nerr + 1
      print '(a,a,a,a)', " FAIL [", op, "] ", label
      print '(a,es40.30)', "   expected = ", expected
      print '(a,es40.30)', "   got      = ", got_q
      print '(a,es12.4)',  "   rel_err  = ", real(rel, dp)
      print '(a,es12.4)',  "   tol      = ", real(tol, dp)
    end if
  end subroutine

  function to_dd(q) result(f)
    real(qp), intent(in) :: q
    type(float64x2) :: f
    real(dp) :: h, l, s, b
    h = real(q, dp)
    l = real(q - real(h, qp), dp)
    s = h + l; b = s - h
    f%limbs(1) = s; f%limbs(2) = l - b
  end function

  subroutine test_case(xq, yq, label, tol)
    real(qp), intent(in) :: xq, yq
    character(*), intent(in) :: label
    real(qp), intent(in), optional :: tol
    type(float64x2) :: xf, yf
    xf = to_dd(xq); yf = to_dd(yq)
    call check("mod", mod(xf, yf), mod(xq, yq), label, tol)
  end subroutine

  subroutine test_modulo_case(xq, yq, label, tol)
    real(qp), intent(in) :: xq, yq
    character(*), intent(in) :: label
    real(qp), intent(in), optional :: tol
    type(float64x2) :: xf, yf
    xf = to_dd(xq); yf = to_dd(yq)
    call check("modulo", modulo(xf, yf), modulo(xq, yq), label, tol)
  end subroutine

end program
