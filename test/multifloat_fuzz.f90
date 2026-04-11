program multifloat_fuzz
  use multifloats
  use, intrinsic :: ieee_arithmetic
  use, intrinsic :: iso_fortran_env, only: int8, int16, int32, int64
  implicit none

  integer, parameter :: qp = 16
  integer, parameter :: dp = 8
  integer, parameter :: sp = 4
  integer, parameter :: num_iterations = 1000000
  integer :: i, num_errors = 0
  real(qp) :: q1, q2, qres
  type(float64x2) :: f1, f2, fres
  double precision :: d1, d2

  real(qp) :: inf_q, nan_q
  double precision :: inf_d, nan_d

  ! Per-op precision statistics (tracked across all check() calls).
  integer, parameter :: max_stats = 256
  type :: stat_entry
    character(len=20) :: name = ''
    real(qp) :: max_rel = 0.0_qp
    real(qp) :: sum_rel = 0.0_qp
    integer(8) :: count = 0_8
  end type stat_entry
  type(stat_entry) :: stats(max_stats)
  integer :: nstats = 0

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

    ! ----------------------------------------------------------------
    ! Hot loop: arithmetic + sqrt + comparisons + unary basics
    ! ----------------------------------------------------------------

    qres = q1 + q2
    fres = f1 + f2
    call check(fres, qres, "add", q1, q2, num_errors)

    qres = q1 - q2
    fres = f1 - f2
    call check(fres, qres, "sub", q1, q2, num_errors)

    qres = q1 * q2
    fres = f1 * f2
    call check(fres, qres, "mul", q1, q2, num_errors)

    if (q2 /= 0.0_qp) then
      qres = q1 / q2
      fres = f1 / f2
      call check(fres, qres, "div", q1, q2, num_errors)
    end if

    if (q1 >= 0.0_qp) then
      qres = sqrt(q1)
      fres = sqrt(f1)
      call check(fres, qres, "sqrt", q1, 0.0_qp, num_errors)
    end if

    qres = abs(q1)
    fres = abs(f1)
    call check(fres, qres, "abs", q1, 0.0_qp, num_errors)

    qres = -q1
    fres = -f1
    call check(fres, qres, "neg", q1, 0.0_qp, num_errors)

    ! Bit-exact rounding/manipulation (full DD)
    if (ieee_is_finite(real(q1, 8)) .and. abs(q1) < 1.0e15_qp) then
      call check(aint(f1), aint(q1), "aint", q1, 0.0_qp, num_errors)
      call check(anint(f1), anint(q1), "anint", q1, 0.0_qp, num_errors)
    end if
    if (ieee_is_finite(real(q1, 8)) .and. q1 /= 0.0_qp) then
      call check(fraction(f1), fraction(q1), "fraction", q1, 0.0_qp, num_errors)
    end if

    call check_comp(f1, f2, q1, q2, num_errors)

    ! ----------------------------------------------------------------
    ! Periodic (every 10): mixed mode + binary functions
    ! ----------------------------------------------------------------
    if (mod(i, 10) == 0) then
       qres = q1 + real(d2, qp)
       fres = f1 + d2
       call check(fres, qres, "add_fd", q1, real(d2, qp), num_errors)

       qres = real(d1, qp) * q2
       fres = d1 * f2
       call check(fres, qres, "mul_df", real(d1, qp), q2, num_errors)

       ! min, max
       call check(min(f1, f2), min(q1, q2), "min", q1, q2, num_errors)
       call check(max(f1, f2), max(q1, q2), "max", q1, q2, num_errors)

       ! sign, dim
       call check(sign(f1, f2), sign(q1, q2), "sign", q1, q2, num_errors)
       call check(dim(f1, f2), dim(q1, q2), "dim", q1, q2, num_errors)

       ! mod, modulo (only for finite divisors away from 0)
       if (ieee_is_finite(real(q1, 8)) .and. ieee_is_finite(real(q2, 8)) .and. &
           q2 /= 0.0_qp .and. abs(q1) < 1e20_qp .and. abs(q2) > 1e-20_qp) then
         call check(mod(f1, f2), mod(q1, q2), "mod", q1, q2, num_errors)
         call check(modulo(f1, f2), modulo(q1, q2), "modulo", q1, q2, num_errors)
       end if

       ! hypot (skip when both inputs are non-finite to avoid spurious NaN traps)
       call check(hypot(f1, f2), hypot(q1, q2), "hypot", q1, q2, num_errors)
    end if

    ! ----------------------------------------------------------------
    ! Periodic (every 100): unary and binary transcendentals
    ! ----------------------------------------------------------------
    if (mod(i, 100) == 0) then
       ! exp / log already in this block
       if (ieee_is_finite(real(q1, 8))) then
         qres = exp(q1)
         fres = exp(f1)
         if (ieee_is_finite(real(qres, 8))) call check(fres, qres, "exp", q1, 0.0_qp, num_errors)
       end if

       if (q1 > 0.0_qp .and. ieee_is_finite(real(q1, 8))) then
         qres = log(q1)
         fres = log(f1)
         call check(fres, qres, "log", q1, 0.0_qp, num_errors)

         qres = log10(q1)
         fres = log10(f1)
         call check(fres, qres, "log10", q1, 0.0_qp, num_errors)
       end if

       ! Trig: keep magnitudes moderate so ULP-of-input doesn't dominate
       if (ieee_is_finite(real(q1, 8)) .and. abs(q1) < 1e6_qp) then
         call check(sin(f1), sin(q1), "sin", q1, 0.0_qp, num_errors)
         call check(cos(f1), cos(q1), "cos", q1, 0.0_qp, num_errors)
         if (abs(cos(q1)) > 1e-12_qp) then
           call check(tan(f1), tan(q1), "tan", q1, 0.0_qp, num_errors)
         end if
       end if

       ! Inverse trig
       if (ieee_is_finite(real(q1, 8)) .and. abs(q1) <= 1.0_qp) then
         call check(asin(f1), asin(q1), "asin", q1, 0.0_qp, num_errors)
         call check(acos(f1), acos(q1), "acos", q1, 0.0_qp, num_errors)
       end if
       if (ieee_is_finite(real(q1, 8))) then
         call check(atan(f1), atan(q1), "atan", q1, 0.0_qp, num_errors)
       end if

       ! Hyperbolic (avoid overflow region)
       if (ieee_is_finite(real(q1, 8)) .and. abs(q1) < 700.0_qp) then
         call check(sinh(f1), sinh(q1), "sinh", q1, 0.0_qp, num_errors)
         call check(cosh(f1), cosh(q1), "cosh", q1, 0.0_qp, num_errors)
         call check(tanh(f1), tanh(q1), "tanh", q1, 0.0_qp, num_errors)
       end if
       if (ieee_is_finite(real(q1, 8))) then
         call check(asinh(f1), asinh(q1), "asinh", q1, 0.0_qp, num_errors)
       end if
       if (ieee_is_finite(real(q1, 8)) .and. q1 >= 1.0_qp) then
         call check(acosh(f1), acosh(q1), "acosh", q1, 0.0_qp, num_errors)
       end if
       if (ieee_is_finite(real(q1, 8)) .and. abs(q1) < 1.0_qp) then
         call check(atanh(f1), atanh(q1), "atanh", q1, 0.0_qp, num_errors)
       end if

       ! Error / gamma family
       if (ieee_is_finite(real(q1, 8)) .and. abs(q1) < 100.0_qp) then
         call check(erf(f1), erf(q1), "erf", q1, 0.0_qp, num_errors)
         call check(erfc(f1), erfc(q1), "erfc", q1, 0.0_qp, num_errors)
         call check(erfc_scaled(f1), erfc_scaled(q1), "erfc_scaled", q1, 0.0_qp, num_errors)
       end if
       if (ieee_is_finite(real(q1, 8)) .and. q1 > 0.0_qp .and. q1 < 100.0_qp) then
         call check(gamma(f1), gamma(q1), "gamma", q1, 0.0_qp, num_errors)
         call check(log_gamma(f1), log_gamma(q1), "lgamma", q1, 0.0_qp, num_errors)
       end if

       ! Bessel (gfortran intrinsics; expect single-double precision)
       if (ieee_is_finite(real(q1, 8)) .and. abs(q1) < 1e3_qp) then
         call check(bessel_j0(f1), bessel_j0(q1), "bj0", q1, 0.0_qp, num_errors)
         call check(bessel_j1(f1), bessel_j1(q1), "bj1", q1, 0.0_qp, num_errors)
         call check(bessel_jn(3, f1), bessel_jn(3, q1), "bjn", q1, 0.0_qp, num_errors)
         if (q1 > 0.0_qp) then
           call check(bessel_y0(f1), bessel_y0(q1), "by0", q1, 0.0_qp, num_errors)
           call check(bessel_y1(f1), bessel_y1(q1), "by1", q1, 0.0_qp, num_errors)
           call check(bessel_yn(3, f1), bessel_yn(3, q1), "byn", q1, 0.0_qp, num_errors)
         end if
       end if

       ! atan2 (works for any pair, including non-finite)
       call check(atan2(f1, f2), atan2(q1, q2), "atan2", q1, q2, num_errors)

       ! Power operator: x**y for positive base, modest exponent.
       ! Range chosen so q1**q2 stays well inside dp range, otherwise
       ! the kernel under/overflows and the recorded rel_err is ~1.
       if (ieee_is_finite(real(q1, 8)) .and. q1 > 1.0e-3_qp .and. q1 < 1.0e3_qp .and. &
           ieee_is_finite(real(q2, 8)) .and. abs(q2) < 30.0_qp) then
         call check(f1 ** f2, q1 ** q2, "pow", q1, q2, num_errors)
         ! Mixed-type pow: mf**dp, dp**mf
         call check(f1 ** real(q2, 8), q1 ** real(q2, 8), "pow_md", q1, q2, num_errors)
         call check(real(q1, 8) ** f2, real(q1, 8) ** q2, "pow_dm", q1, q2, num_errors)
       end if
       if (ieee_is_finite(real(q1, 8)) .and. abs(q1) < 1e10_qp) then
         call check(f1 ** 3, q1 ** 3, "pow_int", q1, 3.0_qp, num_errors)
       end if

       ! scale(x, k) = x * 2^k — exact, full DD. Note: cast the dp
       ! reference up to qp so check() can compare in its native type.
       if (ieee_is_finite(real(q1, 8))) then
         call check(scale(f1, 5), real(scale(real(q1, 8), 5), qp), &
             "scale", q1, 0.0_qp, num_errors)
         if (q1 /= 0.0_qp) then
           call check(set_exponent(f1, 5), &
               real(set_exponent(real(q1, 8), 5), qp), &
               "set_exponent", q1, 0.0_qp, num_errors)
         end if
       end if

       ! 3-argument min/max
       if (ieee_is_finite(real(q1, 8)) .and. ieee_is_finite(real(q2, 8))) then
         block
           type(float64x2) :: f3
           real(qp) :: q3
           q3 = (q1 + q2) * 0.5_qp
           f3 = to_f64x2(q3)
           call check(min(f1, f2, f3), min(q1, q2, q3), "min3", q1, q2, num_errors)
           call check(max(f1, f2, f3), max(q1, q2, q3), "max3", q1, q2, num_errors)
         end block
       end if
    end if

    ! ----------------------------------------------------------------
    ! Periodic (every 200): complex arithmetic and transcendentals
    ! ----------------------------------------------------------------
    if (mod(i, 200) == 0) then
       call fuzz_complex(f1, f2, q1, q2, num_errors)
    end if

    ! ----------------------------------------------------------------
    ! Periodic (every 1000): small-array reductions
    ! ----------------------------------------------------------------
    if (mod(i, 1000) == 0) then
       call fuzz_arrays(num_errors)
    end if

    ! ----------------------------------------------------------------
    ! Periodic (every 100): assignment-operator round-trips
    ! ----------------------------------------------------------------
    if (mod(i, 100) == 0 .and. ieee_is_finite(real(q1, 8))) then
       call fuzz_assignments(f1, q1, num_errors)
    end if

    if (mod(i, 100000) == 0) print *, "Completed", i, "iterations..."
    if (num_errors > 100) then
      print *, "Too many errors, stopping fuzz test."
      exit
    end if
  end do

  call print_all_stats()

  if (num_errors == 0) then
    print *, "FUZZ TEST PASSED"
  else
    print *, "FUZZ TEST FAILED WITH", num_errors, "ERRORS"
    stop 1
  end if

contains

  subroutine update_stats(op, rel_err)
    character(*), intent(in) :: op
    real(qp), intent(in) :: rel_err
    integer :: k
    if (.not. ieee_is_finite(real(rel_err, 8))) return
    do k = 1, nstats
      if (stats(k)%name == op) then
        if (rel_err > stats(k)%max_rel) stats(k)%max_rel = rel_err
        stats(k)%sum_rel = stats(k)%sum_rel + rel_err
        stats(k)%count = stats(k)%count + 1_8
        return
      end if
    end do
    if (nstats >= max_stats) return  ! out of slots; silently drop
    nstats = nstats + 1
    stats(nstats)%name = op
    stats(nstats)%max_rel = rel_err
    stats(nstats)%sum_rel = rel_err
    stats(nstats)%count = 1_8
  end subroutine

  subroutine print_all_stats()
    integer :: k
    real(qp) :: mean_rel
    print *, ""
    print *, "Per-operation precision report (1M iterations):"
    print *, "  ", "op                  ", "      n  ", "    max_rel  ", "    mean_rel"
    do k = 1, nstats
      if (stats(k)%count > 0_8) then
        mean_rel = stats(k)%sum_rel / real(stats(k)%count, qp)
        write(*, '(2x, a20, 1x, i10, 1x, es14.4, 1x, es14.4)') &
            stats(k)%name, stats(k)%count, real(stats(k)%max_rel, 8), real(mean_rel, 8)
      else
        write(*, '(2x, a20, 1x, i10, 1x, a)') stats(k)%name, stats(k)%count, "   (no data)"
      end if
    end do
    print *, ""
  end subroutine

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
    block
      ! Inline fast_two_sum normalization (assumes |h| >= |l|).
      real(8) :: s, b_prime
      s = h + l
      b_prime = s - h
      to_f64x2%limbs(1) = s
      to_f64x2%limbs(2) = l - b_prime
    end block
  end function

  ! Operations that should produce full ~106-bit DD precision results.
  ! cx_div and cx_sqrt are intentionally NOT here: they involve internal
  ! subtractive cancellation that costs a few digits from the upper bound.
  logical function is_full_dd(op)
    character(*), intent(in) :: op
    select case (op)
    case ("add", "sub", "mul", "div", "sqrt", "abs", "neg", &
          "add_fd", "mul_df", "min", "max", "sign", "dim", &
          "mod", "modulo", "hypot", "pow_int", "pow", "pow_md", "pow_dm", &
          "exp", "log", "log10", "sinh", "cosh", &
          "asin", "acos", "acosh", "atan2", &
          "cx_log_re", "cx_sinh_re", "cx_sinh_im", &
          "cx_cosh_re", "cx_cosh_im", &
          "cx_sin_re", "cx_sin_im", "cx_cos_re", "cx_cos_im", &
          "cx_tan_re", "cx_tan_im", "cx_tanh_re", "cx_tanh_im", &
          "aint", "anint", "fraction", "scale", "set_exponent", &
          "min3", "max3", &
          "arr_sum", "arr_max", "arr_min", "arr_dot", "arr_norm2", "arr_matmul", &
          "cx_add_re", "cx_add_im", "cx_sub_re", "cx_sub_im", &
          "cx_mul_re", "cx_mul_im", "cx_conjg_re", "cx_conjg_im", &
          "cx_abs", "cx_aimag", &
          "asn_mf_dp", "asn_mf_sp", "asn_mf_int", "asn_mf_i8", &
          "asn_mf_i16", "asn_mf_i64", "asn_mf_cdp", "asn_mf_csp", &
          "asn_cx_dp_re", "asn_cx_dp_im", "asn_cx_int_re", "asn_cx_int_im", &
          "asn_cx_cdp_re", "asn_cx_cdp_im", &
          "asn_dp_mf", "asn_sp_mf", "asn_cdp_mf_re", "asn_cdp_mf_im", &
          "asn_csp_mf_re", "asn_int_mf", "asn_i64_mf", &
          "ctor_mf_dp", "ctor_mf_sp", "ctor_mf_int", "ctor_mf_i8", &
          "ctor_mf_i16", "ctor_mf_i64", "ctor_mf_cdp", "ctor_mf_csp", &
          "ctor_cx_dp_re", "ctor_cx_dp_im", "ctor_cx_sp_re", &
          "ctor_cx_cdp_re", "ctor_cx_cdp_im", &
          "ctor_cx_dp_dp_re", "ctor_cx_dp_dp_im", &
          "ctor_cx_mf_dp_re", "ctor_cx_mf_dp_im", &
          "ctor_cx_dp_mf_re", "ctor_cx_dp_mf_im")
      is_full_dd = .true.
    case default
      is_full_dd = .false.
    end select
  end function

  ! Compound transcendentals like pow = exp(b*log(a)) and complex
  ! transcendentals chain two derivative-corrected steps and accumulate
  ! ~10× more error than a single transcendental. Bessel intrinsics in
  ! libm are also typically less than full single-double accurate.
  logical function is_compound(op)
    character(*), intent(in) :: op
    select case (op)
    case ("arr_prod", &
          "bj0", "bj1", "by0", "by1", "bjn", "byn", &
          "gamma", "lgamma", "erfc_scaled", &
          "cx_exp_re", "cx_exp_im", "cx_log_im", &
          "cx_sin_re", "cx_sin_im", "cx_cos_re", "cx_cos_im", &
          "cx_sinh_re", "cx_sinh_im", "cx_cosh_re", "cx_cosh_im", &
          "cx_tan_re", "cx_tan_im", "cx_tanh_re", "cx_tanh_im", &
          "cx_atan_re", "cx_atan_im", "cx_asin_re", "cx_asin_im", &
          "cx_acos_re", "cx_acos_im", "cx_asinh_re", "cx_asinh_im", &
          "cx_acosh_re", "cx_acosh_im", "cx_atanh_re", "cx_atanh_im", &
          "cx_div_re", "cx_div_im", "cx_sqrt_re", "cx_sqrt_im")
      is_compound = .true.
    case default
      is_compound = .false.
    end select
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
    rel_err = 0.0_qp
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
        if (is_full_dd(op)) then
          tol = 1e-26_qp
        else if (is_compound(op)) then
          ! Compound transcendentals (pow, bessel, complex non-elementary)
          ! chain two derivative-corrected steps and are subject to
          ! cancellation; allow ~1e-10 relative.
          tol = 1e-10_qp
        else
          tol = 1e-15_qp
        end if
      else
        ! For results near zero, use error relative to input magnitude.
        rel_err = diff / input_mag
        if (is_full_dd(op)) then
          tol = 1e-28_qp
        else if (is_compound(op)) then
          tol = 1e-10_qp
        else
          tol = 1e-15_qp
        end if
      end if

      ! Record the rel_err in the per-op precision report. Skip the
      ! subnormal-input range where the lo limb is below the dp normal
      ! range and DD precision is unattainable in principle.
      if (.not. ((abs(i1) > 0.0_qp .and. abs(i1) < 1.0e-290_qp) .or. &
                 (abs(i2) > 0.0_qp .and. abs(i2) < 1.0e-290_qp) .or. &
                 (abs(q)  > 0.0_qp .and. abs(q)  < 1.0e-290_qp))) then
        call update_stats(op, rel_err)
      end if

      if (rel_err > tol) then
        if (abs(q) > huge(1.0d0)*0.99_qp .or. &
            (abs(q) < tiny(1.0d0) .and. abs(q) > 0.0_qp)) return
        ! When |x| < 2^-970 ≈ 1e-292, the dp ulp at x becomes subnormal,
        ! so the DD low limb cannot hold a full 53-bit error term and the
        ! effective DD precision degrades to single-double. Skip the
        ! precision check (NaN/Inf propagation has already been verified).
        if ((abs(i1) > 0.0_qp .and. abs(i1) < 1.0e-290_qp) .or. &
            (abs(i2) > 0.0_qp .and. abs(i2) < 1.0e-290_qp) .or. &
            (abs(q) > 0.0_qp .and. abs(q) < 1.0e-290_qp)) return
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

  ! Complex variant: check that real and imaginary parts both pass.
  subroutine check_cx(c, q, op, ai_re, ai_im, bi_re, bi_im, errs)
    type(complex128x2), intent(in) :: c
    complex(qp), intent(in) :: q
    character(*), intent(in) :: op
    real(qp), intent(in) :: ai_re, ai_im, bi_re, bi_im
    integer, intent(inout) :: errs
    real(qp) :: input_mag
    input_mag = max(abs(ai_re), abs(ai_im), abs(bi_re), abs(bi_im), 1e-300_qp)
    call check(c%re, real(q, qp), op//"_re", input_mag, 0.0_qp, errs)
    call check(c%im, aimag(q), op//"_im", input_mag, 0.0_qp, errs)
  end subroutine

  subroutine fuzz_complex(f1, f2, q1, q2, errs)
    type(float64x2), intent(in) :: f1, f2
    real(qp), intent(in) :: q1, q2
    integer, intent(inout) :: errs
    type(complex128x2) :: cf1, cf2, cfres
    complex(qp) :: cq1, cq2, cqres

    cf1 = complex128x2(f1, f2)
    cf2 = complex128x2(f2, f1)
    cq1 = cmplx(q1, q2, qp)
    cq2 = cmplx(q2, q1, qp)

    cqres = cq1 + cq2
    cfres = cf1 + cf2
    call check_cx(cfres, cqres, "cx_add", q1, q2, q2, q1, errs)

    cqres = cq1 - cq2
    cfres = cf1 - cf2
    call check_cx(cfres, cqres, "cx_sub", q1, q2, q2, q1, errs)

    cqres = cq1 * cq2
    cfres = cf1 * cf2
    call check_cx(cfres, cqres, "cx_mul", q1, q2, q2, q1, errs)

    if (cq2 /= cmplx(0, 0, qp) .and. ieee_is_finite(real(q1, 8)) .and. &
        ieee_is_finite(real(q2, 8))) then
      cqres = cq1 / cq2
      cfres = cf1 / cf2
      call check_cx(cfres, cqres, "cx_div", q1, q2, q2, q1, errs)
    end if

    if (ieee_is_finite(real(q1, 8)) .and. ieee_is_finite(real(q2, 8))) then
      cfres = sqrt(cf1)
      cqres = sqrt(cq1)
      call check_cx(cfres, cqres, "cx_sqrt", q1, q2, 0.0_qp, 0.0_qp, errs)

      if (abs(q1) < 100.0_qp .and. abs(q2) < 100.0_qp) then
        cfres = exp(cf1)
        cqres = exp(cq1)
        if (ieee_is_finite(real(real(cqres), 8)) .and. &
            ieee_is_finite(real(aimag(cqres), 8))) then
          call check_cx(cfres, cqres, "cx_exp", q1, q2, 0.0_qp, 0.0_qp, errs)
        end if
      end if

      if (cq1 /= cmplx(0, 0, qp)) then
        cfres = log(cf1)
        cqres = log(cq1)
        call check_cx(cfres, cqres, "cx_log", q1, q2, 0.0_qp, 0.0_qp, errs)
      end if

      if (abs(q1) < 100.0_qp .and. abs(q2) < 100.0_qp) then
        cfres = sin(cf1)
        cqres = sin(cq1)
        call check_cx(cfres, cqres, "cx_sin", q1, q2, 0.0_qp, 0.0_qp, errs)

        cfres = cos(cf1)
        cqres = cos(cq1)
        call check_cx(cfres, cqres, "cx_cos", q1, q2, 0.0_qp, 0.0_qp, errs)

        cfres = sinh(cf1)
        cqres = sinh(cq1)
        call check_cx(cfres, cqres, "cx_sinh", q1, q2, 0.0_qp, 0.0_qp, errs)

        cfres = cosh(cf1)
        cqres = cosh(cq1)
        call check_cx(cfres, cqres, "cx_cosh", q1, q2, 0.0_qp, 0.0_qp, errs)

        ! tan/tanh: skip when sinh(b)/cosh(b) overwhelms the leading
        ! limb so the real-part cancellation isn't catastrophic.
        if (abs(q1) < 3.0_qp .and. abs(q2) < 3.0_qp .and. &
            abs(cos(real(q2, 8))) > 0.01d0) then
          cfres = tan(cf1)
          cqres = tan(cq1)
          call check_cx(cfres, cqres, "cx_tan", q1, q2, 0.0_qp, 0.0_qp, errs)

          cfres = tanh(cf1)
          cqres = tanh(cq1)
          call check_cx(cfres, cqres, "cx_tanh", q1, q2, 0.0_qp, 0.0_qp, errs)
        end if
      end if

      ! Inverse trig / hyperbolic — small but non-tiny input range to
      ! avoid the cancellation regimes near branch cuts and near zero
      ! that the first-order derivative correction can't follow.
      if (abs(q1) < 0.5_qp .and. abs(q2) < 0.5_qp .and. &
          max(abs(q1), abs(q2)) > 0.05_qp) then
        cfres = atan(cf1)
        cqres = atan(cq1)
        call check_cx(cfres, cqres, "cx_atan", q1, q2, 0.0_qp, 0.0_qp, errs)

        cfres = asin(cf1)
        cqres = asin(cq1)
        call check_cx(cfres, cqres, "cx_asin", q1, q2, 0.0_qp, 0.0_qp, errs)

        cfres = acos(cf1)
        cqres = acos(cq1)
        call check_cx(cfres, cqres, "cx_acos", q1, q2, 0.0_qp, 0.0_qp, errs)

        cfres = asinh(cf1)
        cqres = asinh(cq1)
        call check_cx(cfres, cqres, "cx_asinh", q1, q2, 0.0_qp, 0.0_qp, errs)

        cfres = acosh(cf1)
        cqres = acosh(cq1)
        call check_cx(cfres, cqres, "cx_acosh", q1, q2, 0.0_qp, 0.0_qp, errs)

        cfres = atanh(cf1)
        cqres = atanh(cq1)
        call check_cx(cfres, cqres, "cx_atanh", q1, q2, 0.0_qp, 0.0_qp, errs)
      end if

      ! conjg, abs, aimag — full DD
      cfres = conjg(cf1)
      cqres = conjg(cq1)
      call check_cx(cfres, cqres, "cx_conjg", q1, q2, 0.0_qp, 0.0_qp, errs)

      call check(abs(cf1), abs(cq1), "cx_abs", q1, q2, errs)
      call check(aimag(cf1), aimag(cq1), "cx_aimag", q1, q2, errs)
    end if
  end subroutine

  ! Assignment-operator round-trips. Each direction must match the
  ! semantics of the equivalent dp ↔ X assignment.
  subroutine fuzz_assignments(f, q, errs)
    type(float64x2), intent(in) :: f
    real(qp), intent(in) :: q
    integer, intent(inout) :: errs
    type(float64x2) :: x_dp, x_sp, x_int, x_i8, x_i16, x_i64, x_cdp, x_csp
    type(complex128x2) :: z_dp, z_int, z_cdp
    real(8) :: out_dp
    real(sp) :: out_sp
    integer :: out_int
    integer(int8) :: out_i8
    integer(int16) :: out_i16
    integer(int64) :: out_i64
    complex(dp) :: out_cdp
    complex(sp) :: out_csp
    real(8) :: hi
    integer(int64) :: q_trunc

    hi = f%limbs(1)

    ! ---------------- mf <- built-in (forward) ----------------
    x_dp = hi
    call check(x_dp, real(hi, qp), "asn_mf_dp", q, 0.0_qp, errs)

    x_sp = real(hi, sp)
    call check(x_sp, real(real(hi, sp), qp), "asn_mf_sp", q, 0.0_qp, errs)

    if (abs(hi) < 2.0e9_dp) then
      x_int = int(hi)
      call check(x_int, real(int(hi), qp), "asn_mf_int", q, 0.0_qp, errs)
    end if
    if (abs(hi) < 100.0_dp) then
      x_i8 = int(hi, int8)
      call check(x_i8, real(int(hi, int8), qp), "asn_mf_i8", q, 0.0_qp, errs)
    end if
    if (abs(hi) < 30000.0_dp) then
      x_i16 = int(hi, int16)
      call check(x_i16, real(int(hi, int16), qp), "asn_mf_i16", q, 0.0_qp, errs)
    end if
    if (abs(hi) < 1.0e18_dp) then
      x_i64 = int(hi, int64)
      call check(x_i64, real(int(hi, int64), qp), "asn_mf_i64", q, 0.0_qp, errs)
    end if

    x_cdp = cmplx(hi, 1.5_dp, dp)
    call check(x_cdp, real(hi, qp), "asn_mf_cdp", q, 0.0_qp, errs)
    x_csp = cmplx(real(hi, sp), 1.5_sp, sp)
    call check(x_csp, real(real(hi, sp), qp), "asn_mf_csp", q, 0.0_qp, errs)

    ! ---------------- cx <- built-in (forward) ----------------
    z_dp = hi
    call check(z_dp%re, real(hi, qp), "asn_cx_dp_re", q, 0.0_qp, errs)
    call check(z_dp%im, 0.0_qp, "asn_cx_dp_im", q, 0.0_qp, errs)
    if (abs(hi) < 2.0e9_dp) then
      z_int = int(hi)
      call check(z_int%re, real(int(hi), qp), "asn_cx_int_re", q, 0.0_qp, errs)
      call check(z_int%im, 0.0_qp, "asn_cx_int_im", q, 0.0_qp, errs)
    end if
    z_cdp = cmplx(hi, 2.5_dp, dp)
    call check(z_cdp%re, real(hi, qp), "asn_cx_cdp_re", q, 0.0_qp, errs)
    call check(z_cdp%im, 2.5_qp, "asn_cx_cdp_im", q, 0.0_qp, errs)

    ! ---------------- mf -> built-in (reverse) ----------------
    out_dp = f
    call check(float64x2(out_dp), real(hi, qp), "asn_dp_mf", q, 0.0_qp, errs)

    out_sp = f
    call check(float64x2(real(out_sp, dp)), real(real(hi, sp), qp), &
        "asn_sp_mf", q, 0.0_qp, errs)

    out_cdp = f
    call check(float64x2(real(out_cdp, dp)), real(hi, qp), "asn_cdp_mf_re", q, 0.0_qp, errs)
    call check(float64x2(aimag(out_cdp)), 0.0_qp, "asn_cdp_mf_im", q, 0.0_qp, errs)

    out_csp = f
    call check(float64x2(real(real(out_csp, sp), dp)), &
        real(real(hi, sp), qp), "asn_csp_mf_re", q, 0.0_qp, errs)

    ! Integer truncation: must match Fortran's `int_var = real_dp_var`
    ! semantics on the leading limb.
    if (abs(hi) < 1.0e18_dp) then
      out_i64 = f
      q_trunc = int(hi, int64)
      ! Adjust for the lo limb the same way the assignment kernel does.
      if (hi == real(q_trunc, dp)) then
        if (q_trunc > 0_int64 .and. f%limbs(2) < 0.0_dp) q_trunc = q_trunc - 1_int64
        if (q_trunc < 0_int64 .and. f%limbs(2) > 0.0_dp) q_trunc = q_trunc + 1_int64
      end if
      call check(float64x2(real(out_i64, dp)), real(q_trunc, qp), &
          "asn_i64_mf", q, 0.0_qp, errs)
    end if
    if (abs(hi) < 2.0e9_dp) then
      out_int = f
      call check(float64x2(real(out_int, dp)), real(int(hi), qp), &
          "asn_int_mf", q, 0.0_qp, errs)
    end if

    ! ---------------- float64x2(...) constructors ----------------
    call check(float64x2(hi), real(hi, qp), "ctor_mf_dp", q, 0.0_qp, errs)
    call check(float64x2(real(hi, sp)), real(real(hi, sp), qp), &
        "ctor_mf_sp", q, 0.0_qp, errs)
    if (abs(hi) < 2.0e9_dp) then
      call check(float64x2(int(hi)), real(int(hi), qp), &
          "ctor_mf_int", q, 0.0_qp, errs)
    end if
    if (abs(hi) < 100.0_dp) then
      call check(float64x2(int(hi, int8)), real(int(hi, int8), qp), &
          "ctor_mf_i8", q, 0.0_qp, errs)
    end if
    if (abs(hi) < 30000.0_dp) then
      call check(float64x2(int(hi, int16)), real(int(hi, int16), qp), &
          "ctor_mf_i16", q, 0.0_qp, errs)
    end if
    if (abs(hi) < 1.0e18_dp) then
      call check(float64x2(int(hi, int64)), real(int(hi, int64), qp), &
          "ctor_mf_i64", q, 0.0_qp, errs)
    end if
    call check(float64x2(cmplx(hi, 1.5_dp, dp)), real(hi, qp), &
        "ctor_mf_cdp", q, 0.0_qp, errs)
    call check(float64x2(cmplx(real(hi, sp), 1.5_sp, sp)), &
        real(real(hi, sp), qp), "ctor_mf_csp", q, 0.0_qp, errs)

    ! ---------------- complex128x2(...) constructors ----------------
    block
      type(complex128x2) :: ztmp
      ztmp = complex128x2(hi)
      call check(ztmp%re, real(hi, qp), "ctor_cx_dp_re", q, 0.0_qp, errs)
      call check(ztmp%im, 0.0_qp, "ctor_cx_dp_im", q, 0.0_qp, errs)
      ztmp = complex128x2(real(hi, sp))
      call check(ztmp%re, real(real(hi, sp), qp), "ctor_cx_sp_re", q, 0.0_qp, errs)
      ztmp = complex128x2(cmplx(hi, 2.5_dp, dp))
      call check(ztmp%re, real(hi, qp), "ctor_cx_cdp_re", q, 0.0_qp, errs)
      call check(ztmp%im, 2.5_qp, "ctor_cx_cdp_im", q, 0.0_qp, errs)
      ! Two-arg matching kinds
      ztmp = complex128x2(hi, hi)
      call check(ztmp%re, real(hi, qp), "ctor_cx_dp_dp_re", q, 0.0_qp, errs)
      call check(ztmp%im, real(hi, qp), "ctor_cx_dp_dp_im", q, 0.0_qp, errs)
      ! Mixed (mf, dp)
      ztmp = complex128x2(f, hi)
      call check(ztmp%re, real(hi, qp), "ctor_cx_mf_dp_re", q, 0.0_qp, errs)
      call check(ztmp%im, real(hi, qp), "ctor_cx_mf_dp_im", q, 0.0_qp, errs)
      ! Mixed (dp, mf)
      ztmp = complex128x2(hi, f)
      call check(ztmp%re, real(hi, qp), "ctor_cx_dp_mf_re", q, 0.0_qp, errs)
      call check(ztmp%im, real(hi, qp), "ctor_cx_dp_mf_im", q, 0.0_qp, errs)
    end block
  end subroutine

  ! Small-array reductions exercised every 1000 iterations.
  subroutine fuzz_arrays(errs)
    integer, intent(inout) :: errs
    integer, parameter :: nn = 8
    type(float64x2) :: a(nn), b(nn), m(nn,nn)
    real(qp) :: qa(nn), qb(nn), qm(nn,nn)
    type(float64x2) :: mvres(nn)
    real(qp) :: qmvres(nn)
    real(qp) :: r(2)
    integer :: k, l

    do k = 1, nn
      call random_number(r)
      qa(k) = (r(1) - 0.5_qp) * 10.0_qp ** int(r(2)*20.0_qp - 10.0_qp)
      a(k) = to_f64x2(qa(k))
      call random_number(r)
      qb(k) = (r(1) - 0.5_qp) * 10.0_qp ** int(r(2)*20.0_qp - 10.0_qp)
      b(k) = to_f64x2(qb(k))
    end do
    do l = 1, nn
      do k = 1, nn
        call random_number(r)
        qm(k,l) = (r(1) - 0.5_qp) * 10.0_qp ** int(r(2)*20.0_qp - 10.0_qp)
        m(k,l) = to_f64x2(qm(k,l))
      end do
    end do

    call check(sum(a), sum(qa), "arr_sum", maxval(abs(qa)), 0.0_qp, errs)
    call check(product(a), product(qa), "arr_prod", &
               maxval(abs(qa))**nn, 0.0_qp, errs)
    call check(maxval(a), maxval(qa), "arr_max", maxval(abs(qa)), 0.0_qp, errs)
    call check(minval(a), minval(qa), "arr_min", maxval(abs(qa)), 0.0_qp, errs)
    call check(dot_product(a, b), dot_product(qa, qb), "arr_dot", &
               maxval(abs(qa)) * maxval(abs(qb)) * nn, 0.0_qp, errs)
    call check(norm2(a), norm2(qa), "arr_norm2", maxval(abs(qa)), 0.0_qp, errs)

    mvres = matmul(m, a)
    qmvres = matmul(qm, qa)
    do k = 1, nn
      call check(mvres(k), qmvres(k), "arr_matmul", &
                 maxval(abs(qm)) * maxval(abs(qa)) * nn, 0.0_qp, errs)
    end do
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
