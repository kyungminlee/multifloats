program multifloat_fuzz
  use multifloat
  use, intrinsic :: ieee_arithmetic
  implicit none

  integer, parameter :: qp = 16
  integer, parameter :: num_iterations = 1000000
  integer :: i, num_errors = 0
  real(qp) :: q1, q2, qres
  type(float64x2) :: f1, f2, fres
  double precision :: d1, d2

  real(qp) :: inf_q, nan_q
  double precision :: inf_d, nan_d

  inf_q = ieee_value(inf_q, ieee_positive_inf)
  nan_q = ieee_value(nan_q, ieee_quiet_nan)
  inf_d = ieee_value(inf_d, ieee_positive_inf)
  nan_d = ieee_value(nan_d, ieee_quiet_nan)

  print *, "Starting fuzz tests with", num_iterations, "iterations..."

  do i = 1, num_iterations
    ! Generate q1, q2
    call generate_pair(q1, q2)
    f1 = to_f64x2(q1)
    f2 = to_f64x2(q2)
    d1 = real(q1, 8)
    d2 = real(q2, 8)

    ! Test Addition
    qres = q1 + q2
    fres = f1 + f2
    call check(fres, qres, "add", q1, q2, num_errors)

    ! Test Subtraction
    qres = q1 - q2
    fres = f1 - f2
    call check(fres, qres, "sub", q1, q2, num_errors)

    ! Test Multiplication
    qres = q1 * q2
    fres = f1 * f2
    call check(fres, qres, "mul", q1, q2, num_errors)

    ! Test Division
    if (q2 /= 0.0_qp) then
      qres = q1 / q2
      fres = f1 / f2
      call check(fres, qres, "div", q1, q2, num_errors)
    end if

    ! Test Sqrt
    if (q1 >= 0.0_qp) then
      qres = sqrt(q1)
      fres = sqrt(f1)
      call check(fres, qres, "sqrt", q1, 0.0_qp, num_errors)
    end if

    ! Test Exp/Log periodically (slower)
    if (mod(i, 100) == 0) then
       qres = exp(q1)
       fres = exp(f1)
       if (ieee_is_finite(real(qres, 8))) call check(fres, qres, "exp", q1, 0.0_qp, num_errors)
       
       if (q1 > 0.0_qp) then
         qres = log(q1)
         fres = log(f1)
         call check(fres, qres, "log", q1, 0.0_qp, num_errors)
       end if
    end if

    ! Test Comparisons
    call check_comp(f1, f2, q1, q2, num_errors)

    ! Mixed mode fuzzing periodically
    if (mod(i, 10) == 0) then
       ! ff + d
       qres = q1 + real(d2, qp)
       fres = f1 + d2
       call check(fres, qres, "add_fd", q1, real(d2, qp), num_errors)

       ! d * ff
       qres = real(d1, qp) * q2
       fres = d1 * f2
       call check(fres, qres, "mul_df", real(d1, qp), q2, num_errors)
    end if

    if (mod(i, 100000) == 0) print *, "Completed", i, "iterations..."
    if (num_errors > 100) then
      print *, "Too many errors, stopping fuzz test."
      exit
    end if
  end do

  if (num_errors == 0) then
    print *, "FUZZ TEST PASSED"
  else
    print *, "FUZZ TEST FAILED WITH", num_errors, "ERRORS"
    stop 1
  end if

contains

  subroutine generate_pair(q1, q2)
    real(qp), intent(out) :: q1, q2
    real(8) :: r(4)
    call random_number(r)

    ! Strategy for interesting cases
    select case (int(r(1) * 10))
    case (0) ! Non-finite
      q1 = pick_nonfinite(r(2))
      q2 = pick_nonfinite(r(3))
    case (1) ! Close numbers
      q1 = (r(2) - 0.5d0) * 10.0d0**int(r(3)*20 - 10)
      q2 = q1 * (1.0_qp + real(r(4)*1d-15, qp))
    case (2) ! Sum near zero
      q1 = (r(2) - 0.5d0) * 10.0d0**int(r(3)*20 - 10)
      q2 = -q1 + real((r(4)-0.5d0)*1d-25, qp) * q1
    case (3) ! Near max double
      q1 = huge(1.0d0) * (0.9d0 + 0.1d0 * r(2))
      q2 = huge(1.0d0) * (0.9d0 + 0.1d0 * r(3))
    case (4) ! Near tiny double
      q1 = tiny(1.0d0) * (1.0d0 + 10.0d0 * r(2))
      q2 = tiny(1.0d0) * (1.0d0 + 10.0d0 * r(3))
    case default ! Random
      q1 = (r(2) - 0.5d0) * 10.0d0**int(r(3)*60 - 30)
      q2 = (r(4) - 0.5d0) * 10.0d0**int(r(1)*60 - 30)
    end select
  end subroutine

  real(qp) function pick_nonfinite(r)
    real(8), intent(in) :: r
    if (r < 0.33d0) then
      pick_nonfinite = inf_q
    else if (r < 0.66d0) then
      pick_nonfinite = -inf_q
    else
      pick_nonfinite = nan_q
    end if
  end function

  type(float64x2) function to_f64x2(q)
    real(qp), intent(in) :: q
    real(8) :: h, l
    if (.not. ieee_is_finite(real(q, 8))) then
      to_f64x2%limbs(1) = real(q, 8)
      to_f64x2%limbs(2) = 0.0d0
      return
    end if
    h = real(q, 8)
    l = real(q - real(h, qp), 8)
    to_f64x2%limbs(1) = h
    to_f64x2%limbs(2) = l
    call renormalize(to_f64x2)
  end function

  subroutine check(f, q, op, i1, i2, errs)
    type(float64x2), intent(in) :: f
    real(qp), intent(in) :: q
    character(*), intent(in) :: op
    real(qp), intent(in) :: i1, i2
    integer, intent(inout) :: errs
    real(qp) :: f_q, diff, rel_err, tol, input_mag
    logical :: failed
    
    failed = .false.
    if (ieee_is_nan(q)) then
      if (.not. ieee_is_nan(f%limbs(1))) failed = .true.
    else if (is_inf_q(q)) then
      if (.not. is_inf_d(f%limbs(1)) .or. &
          sign(1.0_qp, q) /= sign(1.0_qp, real(f%limbs(1), qp))) failed = .true.
    else
      f_q = real(f%limbs(1), qp) + real(f%limbs(2), qp)
      diff = abs(f_q - q)
      input_mag = max(abs(i1), abs(i2), 1e-300_qp)
      
      if (abs(q) > input_mag * 1e-10_qp) then
        rel_err = diff / abs(q)
        if (op == "div" .or. op == "sqrt" .or. op == "exp" .or. op == "log") then
          tol = 1e-15_qp
        else
          tol = 1e-26_qp
        end if
      else
        ! For results near zero, use error relative to input magnitude
        rel_err = diff / input_mag
        tol = 1e-28_qp
      end if

      if (rel_err > tol) then
        if (abs(q) > huge(1.0d0)*0.99_qp .or. &
            (abs(q) < tiny(1.0d0) .and. abs(q) > 0.0_qp)) return
        if (diff < 1e-35_qp) return
        failed = .true.
      end if
    end if

    if (failed) then
      errs = errs + 1
      print *, "FAIL: ", op
      print *, "  Inputs: ", i1, i2
      print *, "  Expected: ", q
      print *, "  Got:      ", f%limbs, " (as qp sum: ", real(f%limbs(1), qp) + real(f%limbs(2), qp), ")"
      print *, "  Rel Err:  ", rel_err
    end if
  end subroutine

  subroutine check_comp(f1, f2, q1, q2, errs)
    type(float64x2), intent(in) :: f1, f2
    real(qp), intent(in) :: q1, q2
    integer, intent(inout) :: errs
    real(qp) :: diff, mag
    
    if (ieee_is_nan(q1) .or. ieee_is_nan(q2)) return

    diff = abs(q1 - q2)
    mag = max(abs(q1), abs(q2), 1e-300_qp)

    if ((f1 < f2) .neqv. (q1 < q2)) then
       if (diff > mag*1e-20_qp .and. diff > 1e-32_qp) then
         call report_comp_fail("lt", q1, q2, f1 < f2, q1 < q2, errs)
       end if
    end if
    if ((f1 == f2) .neqv. (q1 == q2)) then
       if (diff > mag*1e-20_qp .and. diff > 1e-32_qp) then
         call report_comp_fail("eq", q1, q2, f1 == f2, q1 == q2, errs)
       end if
    end if
  end subroutine

  subroutine report_comp_fail(op, q1, q2, fres, qres, errs)
    character(*), intent(in) :: op
    real(qp), intent(in) :: q1, q2
    logical, intent(in) :: fres, qres
    integer, intent(inout) :: errs
    errs = errs + 1
    print *, "FAIL COMP: ", op
    print *, "  Inputs: ", q1, q2
    print *, "  Expected: ", qres
    print *, "  Got:      ", fres
  end subroutine

  logical function is_inf_q(x)
    real(qp), intent(in) :: x
    is_inf_q = (.not. ieee_is_finite(x)) .and. (.not. ieee_is_nan(x))
  end function

  logical function is_inf_d(x)
    double precision, intent(in) :: x
    is_inf_d = (.not. ieee_is_finite(x)) .and. (.not. ieee_is_nan(x))
  end function

end program
