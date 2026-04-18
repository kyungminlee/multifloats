program test_multifloats_precision
  use multifloats
  use, intrinsic :: ieee_arithmetic
  use, intrinsic :: iso_fortran_env, only: int8, int16, int32, int64
  implicit none

  integer, parameter :: qp = 16
  integer, parameter :: dp = 8
  integer, parameter :: sp = 4
  real(qp), parameter :: tol = 1.0e5_qp * epsilon(1.0_qp)
  integer :: failures

  failures = 0
  write(*,'(a)') 'Native Float64x2 vs quad-precision reference'

  call test_scalar_constructors(failures)
  call test_native_scalar_arithmetic(failures)
  call test_native_unary_abs(failures)
  call test_native_comparisons(failures)
  call test_complex_arithmetic(failures)
  call test_array_ops(failures)
  call test_rounding_and_random(failures)
  call test_signed_zero(failures)
  call test_infinities(failures)
  call test_nan_propagation(failures)
  call test_subnormal_boundary(failures)
  call test_huge_boundary(failures)
  call test_ulp_boundary(failures)
  call test_reduction_mask(failures)
  call test_reduction_dim(failures)
  call test_loc_back(failures)
  call test_assignments_to_dd(failures)
  call test_assignments_to_cdd(failures)
  call test_assignments_from_dd(failures)
  call test_assignments_from_cdd(failures)
  call test_constructors(failures)

  if (failures /= 0) then
    write(*,'(a,i0)') 'FAIL: ', failures
    error stop 1
  end if

  write(*,'(a)') 'PASS'

contains

  pure elemental function dd_to_qp(x) result(y)
    type(float64x2), intent(in) :: x
    real(qp) :: y
    y = real(x%limbs(1), qp) + real(x%limbs(2), qp)
  end function

  pure elemental function cdd_to_cqp(x) result(y)
    type(complex64x2), intent(in) :: x
    complex(qp) :: y
    y = cmplx(dd_to_qp(x%re), dd_to_qp(x%im), kind=qp)
  end function

  subroutine check_true(label, cond, failures)
    character(*), intent(in) :: label
    logical, intent(in) :: cond
    integer, intent(inout) :: failures
    if (.not. cond) then
      failures = failures + 1
      write(*,'(a)') 'FAIL '//trim(label)
    end if
  end subroutine

  subroutine check_real_close(label, got, expect, failures, scale)
    character(*), intent(in) :: label
    type(float64x2), intent(in) :: got
    real(qp), intent(in) :: expect
    integer, intent(inout) :: failures
    real(qp), intent(in), optional :: scale
    real(qp) :: err, bound, mag, rel
    if (present(scale)) then
      mag = scale
    else
      mag = max(1.0_qp, abs(expect))
    end if
    bound = tol * mag
    err = abs(dd_to_qp(got) - expect)
    rel = err / max(1.0_qp, abs(expect))
    write(*,'(a,1x,es24.16,1x,es24.16,1x,es24.16)') trim(label), err, rel, bound
    if (err > bound) then
      failures = failures + 1
      write(*,'(a,1x,es24.16,1x,es24.16,1x,es24.16)') 'FAIL '//trim(label), err, bound, expect
    end if
  end subroutine

  subroutine check_complex_close(label, got, expect, failures, scale)
    character(*), intent(in) :: label
    type(complex64x2), intent(in) :: got
    complex(qp), intent(in) :: expect
    integer, intent(inout) :: failures
    real(qp), intent(in), optional :: scale
    real(qp) :: err, bound, mag, rel
    if (present(scale)) then
      mag = scale
    else
      mag = max(1.0_qp, abs(expect))
    end if
    bound = tol * mag
    err = abs(cdd_to_cqp(got) - expect)
    rel = err / max(1.0_qp, abs(expect))
    write(*,'(a,1x,es24.16,1x,es24.16,1x,es24.16)') trim(label), err, rel, bound
    if (err > bound) then
      failures = failures + 1
      write(*,'(a,1x,es24.16,1x,es24.16)') 'FAIL '//trim(label), err, bound
    end if
  end subroutine

  ! Looser tolerance (~ulp(double)) for complex transcendentals, which
  ! are derived from the first-order derivative-corrected real DD kernels.
  subroutine check_complex_approx(label, got, expect, failures)
    character(*), intent(in) :: label
    type(complex64x2), intent(in) :: got
    complex(qp), intent(in) :: expect
    integer, intent(inout) :: failures
    real(qp) :: err, bound
    bound = max(abs(expect), 1.0_qp) * 1.0e-14_qp
    err = abs(cdd_to_cqp(got) - expect)
    write(*,'(a,1x,es24.16,1x,es24.16)') trim(label), err, bound
    if (err > bound) then
      failures = failures + 1
      write(*,'(a,1x,es24.16,1x,es24.16)') 'FAIL '//trim(label), err, bound
    end if
  end subroutine

  subroutine test_scalar_constructors(failures)
    integer, intent(inout) :: failures
    type(float64x2) :: x
    type(complex64x2) :: z

    x = float64x2(1.25_dp)
    call check_real_close('ctor dp', x, 1.25_qp, failures)

    x = float64x2(real(1.25_sp, sp))
    call check_real_close('ctor sp', x, real(1.25_sp, qp), failures)

    x = float64x2(17)
    call check_real_close('ctor int', x, 17.0_qp, failures)

    x = float64x2('1.2345678901234567890123456789')
    call check_real_close('ctor char', x, 1.2345678901234567890123456789_qp, failures)

    x = cmplx(2.5_dp, -9.0_dp, kind=dp)
    call check_real_close('assign complex->real', x, 2.5_qp, failures)

    z = complex64x2(float64x2(1.5_dp), float64x2(-0.25_dp))
    call check_complex_close('ctor complex dd', z, cmplx(1.5_qp, -0.25_qp, qp), failures)

    z = complex64x2(3.0_dp, -2.0_dp)
    call check_complex_close('ctor complex dp', z, cmplx(3.0_qp, -2.0_qp, qp), failures)
  end subroutine

  subroutine test_native_scalar_arithmetic(failures)
    integer, intent(inout) :: failures
    type(float64x2) :: x, y
    real(qp) :: qx, qy

    x%limbs = [1.0_dp, 2.0e-30_dp]
    y%limbs = [2.0_dp, -3.0e-30_dp]
    qx = dd_to_qp(x)
    qy = dd_to_qp(y)

    call check_real_close('add dd-dd', x + y, qx + qy, failures)
    call check_real_close('sub dd-dd', x - y, qx - qy, failures)
    call check_real_close('mul dd-dd', x * y, qx * qy, failures)
    call check_real_close('div dd-dd', x / y, qx / qy, failures)

    call check_real_close('add dd-dp', x + 2.0_dp, qx + 2.0_qp, failures)
    call check_real_close('add dp-dd', 2.0_dp + x, 2.0_qp + qx, failures)
    call check_real_close('sub dd-int', x - 2, qx - 2.0_qp, failures)
    call check_real_close('sub int-dd', 2 - x, 2.0_qp - qx, failures)
    call check_real_close('mul dd-sp', x * real(0.5_sp, sp), qx * 0.5_qp, failures)
    call check_real_close('mul sp-dd', real(0.5_sp, sp) * x, 0.5_qp * qx, failures)
    call check_real_close('div dd-dp', x / 2.0_dp, qx / 2.0_qp, failures)
    call check_real_close('div dp-dd', 2.0_dp / x, 2.0_qp / qx, failures)
    call check_real_close('pow int', x ** 3, qx ** 3, failures, scale=abs(qx ** 3))

  end subroutine

  subroutine test_native_unary_abs(failures)
    integer, intent(inout) :: failures
    type(float64x2) :: x
    real(qp) :: qx

    x%limbs = [1.25_dp, 1.0e-30_dp]
    qx = dd_to_qp(x)

    call check_real_close('abs', abs(-x), abs(-qx), failures)
  end subroutine

  subroutine test_native_comparisons(failures)
    integer, intent(inout) :: failures
    type(float64x2) :: x, y

    x%limbs = [1.0_dp, 2.0e-30_dp]
    y%limbs = [2.0_dp, -3.0e-30_dp]

    call check_true('comparison lt', x < y, failures)
    call check_true('comparison gt', y > x, failures)
    call check_true('comparison eq', float64x2(2) == 2, failures)
    call check_true('comparison ne', x /= y, failures)
  end subroutine

  subroutine test_complex_arithmetic(failures)
    integer, intent(inout) :: failures
    type(complex64x2) :: z1, z2
    complex(qp) :: q1, q2

    z1 = complex64x2(float64x2(1.25_dp), float64x2(-0.5_dp))
    z2 = complex64x2(0.75_dp, 0.25_dp)
    q1 = cdd_to_cqp(z1)
    q2 = cdd_to_cqp(z2)

    call check_complex_close('cdd add', z1 + z2, q1 + q2, failures)
    call check_complex_close('cdd sub mixed', z1 - cmplx(0.5_dp, -0.125_dp, dp), q1 - cmplx(0.5_qp, -0.125_qp, qp), failures)
    call check_complex_close('cdd mul', z1 * z2, q1 * q2, failures)
    call check_complex_close('cdd div', z1 / z2, q1 / q2, failures)
    ! sqrt is full DD precision; exp / log / trig / hyperbolic are
    ! built from the first-order derivative-corrected real kernels and
    ! only carry single-double precision.
    call check_complex_close('cdd sqrt', sqrt(z1), sqrt(q1), failures)
    call check_complex_approx('cdd exp', exp(z1), exp(q1), failures)

    call check_real_close('cdd abs', abs(z1), abs(q1), failures)
    call check_real_close('cdd aimag', aimag(z1), aimag(q1), failures)
  end subroutine

  subroutine test_array_ops(failures)
    integer, intent(inout) :: failures
    type(float64x2) :: a(2), m1(2,2), m2(2,2), mv(2), mm(2,2)
    type(complex64x2) :: cz(2)
    real(qp) :: qa(2), qm1(2,2), qm2(2,2), qv(2), qmm(2,2)
    complex(qp) :: qcz(2)

    a(1)%limbs = [1.0_dp, 1.0e-30_dp]
    a(2)%limbs = [2.0_dp, -2.0e-30_dp]
    qa = [dd_to_qp(a(1)), dd_to_qp(a(2))]

    call check_real_close('sum', sum(a), sum(qa), failures)
    call check_real_close('product', product(a), product(qa), failures, scale=abs(product(qa)))
    call check_real_close('dot_product', dot_product(a, a), dot_product(qa, qa), failures, scale=dot_product(qa, qa))
    call check_real_close('norm2', norm2(a), norm2(qa), failures)
    call check_real_close('maxval', maxval(a), maxval(qa), failures)
    call check_true('maxloc', all(maxloc(a) == maxloc(qa)), failures)
    call check_true('findloc', all(findloc(a, a(2)) == findloc(qa, qa(2))), failures)

    m1 = reshape([float64x2(1.0_dp), float64x2(2.0_dp), float64x2(3.0_dp), float64x2(4.0_dp)], [2,2])
    m2 = reshape([float64x2(0.5_dp), float64x2(-1.0_dp), float64x2(1.5_dp), float64x2(2.0_dp)], [2,2])
    qm1 = reshape([1.0_qp, 2.0_qp, 3.0_qp, 4.0_qp], [2,2])
    qm2 = reshape([0.5_qp, -1.0_qp, 1.5_qp, 2.0_qp], [2,2])

    mm = matmul(m1, m2)
    qmm = matmul(qm1, qm2)
    call check_real_close('matmul(1,1)', mm(1,1), qmm(1,1), failures, scale=10.0_qp)
    call check_real_close('matmul(2,1)', mm(2,1), qmm(2,1), failures, scale=10.0_qp)
    call check_real_close('matmul(1,2)', mm(1,2), qmm(1,2), failures, scale=10.0_qp)
    call check_real_close('matmul(2,2)', mm(2,2), qmm(2,2), failures, scale=10.0_qp)

    mv = matmul(m1, a)
    qv = matmul(qm1, qa)
    call check_real_close('matmul mv 1', mv(1), qv(1), failures, scale=10.0_qp)
    call check_real_close('matmul mv 2', mv(2), qv(2), failures, scale=10.0_qp)

    cz = [complex64x2(1.0_dp, -0.5_dp), complex64x2(0.25_dp, 0.75_dp)]
    qcz = [cdd_to_cqp(cz(1)), cdd_to_cqp(cz(2))]
    call check_complex_close('complex sum', sum(cz), sum(qcz), failures)
    call check_complex_close('complex dot_product', dot_product(cz, cz), dot_product(qcz, qcz), failures)
  end subroutine

  subroutine test_rounding_and_random(failures)
    integer, intent(inout) :: failures
    type(float64x2) :: x, r

    x%limbs = [3.5_dp, 2.0e-30_dp]
    call check_true('dble', abs(dble(x) - real(dd_to_qp(x), dp)) <= epsilon(1.0_dp), failures)
    call check_true('int', int(x) == int(dd_to_qp(x)), failures)
    call check_true('nint', nint(x) == nint(dd_to_qp(x)), failures)
    call check_true('floor', floor(x) == floor(dd_to_qp(x)), failures)
    call check_true('ceiling', ceiling(x) == ceiling(dd_to_qp(x)), failures)

    call random_number(r)
    call check_true('random lower', dd_to_qp(r) >= 0.0_qp, failures)
    call check_true('random upper', dd_to_qp(r) < 1.0_qp, failures)
  end subroutine

  ! ----------------------------------------------------------------
  ! Edge-case sweeps
  ! ----------------------------------------------------------------

  subroutine test_signed_zero(failures)
    integer, intent(inout) :: failures
    type(float64x2) :: pz, nz, x, y
    pz%limbs = [0.0_dp, 0.0_dp]
    nz%limbs = [-0.0_dp, 0.0_dp]

    ! Equality of +0 and -0 (IEEE: +0 == -0)
    call check_true('+0 == -0', pz == nz, failures)
    call check_true('+0 not signbit', .not. (pz%limbs(1) < 0.0_dp), failures)

    ! Sign-bit semantics on lo limb (DD signbit uses limb(1) only)
    call check_true('-0 signbit dp', sign(1.0_dp, nz%limbs(1)) < 0.0_dp, failures)

    ! +0 + +0 = +0, -0 + -0 = -0, +0 + -0 = +0
    x = pz + pz
    call check_true('+0 + +0 = +0', sign(1.0_dp, x%limbs(1)) > 0.0_dp .and. x%limbs(1) == 0.0_dp, failures)
    x = nz + nz
    call check_true('-0 + -0 = -0', sign(1.0_dp, x%limbs(1)) < 0.0_dp .and. x%limbs(1) == 0.0_dp, failures)

    ! abs(-0) = +0
    x = abs(nz)
    call check_true('abs(-0) = +0', x%limbs(1) == 0.0_dp .and. sign(1.0_dp, x%limbs(1)) > 0.0_dp, failures)

    ! 1.0 / +0 = +inf, 1.0 / -0 = -inf
    y = float64x2(1.0_dp)
    x = y / pz
    call check_true('1/+0 = +inf', .not. ieee_is_finite(x%limbs(1)) .and. x%limbs(1) > 0.0_dp, failures)
    x = y / nz
    call check_true('1/-0 = -inf', .not. ieee_is_finite(x%limbs(1)) .and. x%limbs(1) < 0.0_dp, failures)

    ! sqrt(+0) = +0, sqrt(-0) = -0 (IEEE)
    x = sqrt(pz)
    call check_true('sqrt(+0) = +0', x%limbs(1) == 0.0_dp .and. sign(1.0_dp, x%limbs(1)) > 0.0_dp, failures)
    x = sqrt(nz)
    call check_true('sqrt(-0) = -0', x%limbs(1) == 0.0_dp .and. sign(1.0_dp, x%limbs(1)) < 0.0_dp, failures)
  end subroutine

  subroutine test_infinities(failures)
    integer, intent(inout) :: failures
    type(float64x2) :: pinf, ninf, one, zero, x
    real(dp) :: inf
    inf = ieee_value(inf, ieee_positive_inf)
    pinf%limbs = [inf, 0.0_dp]
    ninf%limbs = [-inf, 0.0_dp]
    one = float64x2(1.0_dp)
    zero%limbs = [0.0_dp, 0.0_dp]

    ! +inf + +inf = +inf
    x = pinf + pinf
    call check_true('+inf + +inf = +inf', .not. ieee_is_finite(x%limbs(1)) &
        .and. .not. ieee_is_nan(x%limbs(1)) .and. x%limbs(1) > 0.0_dp, failures)
    ! +inf - +inf = NaN
    x = pinf - pinf
    call check_true('+inf - +inf = NaN', ieee_is_nan(x%limbs(1)), failures)
    ! +inf + -inf = NaN
    x = pinf + ninf
    call check_true('+inf + -inf = NaN', ieee_is_nan(x%limbs(1)), failures)
    ! +inf * 2 = +inf
    x = pinf * float64x2(2.0_dp)
    call check_true('+inf * 2 = +inf', x%limbs(1) > huge(1.0_dp), failures)
    ! +inf * 0 = NaN
    x = pinf * zero
    call check_true('+inf * 0 = NaN', ieee_is_nan(x%limbs(1)), failures)
    ! +inf / +inf = NaN
    x = pinf / pinf
    call check_true('+inf / +inf = NaN', ieee_is_nan(x%limbs(1)), failures)
    ! 1 / +inf = 0
    x = one / pinf
    call check_true('1 / +inf = 0', x%limbs(1) == 0.0_dp, failures)
    ! sqrt(+inf) = +inf
    x = sqrt(pinf)
    call check_true('sqrt(+inf) = +inf', x%limbs(1) > huge(1.0_dp), failures)
    ! sqrt(-inf) = NaN
    x = sqrt(ninf)
    call check_true('sqrt(-inf) = NaN', ieee_is_nan(x%limbs(1)), failures)
    ! Comparison with inf
    call check_true('+inf > 1', pinf > one, failures)
    call check_true('-inf < 1', ninf < one, failures)
    call check_true('+inf == +inf', pinf == pinf, failures)
    ! abs(-inf) = +inf
    x = abs(ninf)
    call check_true('abs(-inf) = +inf', x%limbs(1) > huge(1.0_dp), failures)
  end subroutine

  subroutine test_nan_propagation(failures)
    integer, intent(inout) :: failures
    type(float64x2) :: nan, x, one
    real(dp) :: nan_d
    nan_d = ieee_value(nan_d, ieee_quiet_nan)
    nan%limbs = [nan_d, 0.0_dp]
    one = float64x2(1.0_dp)

    ! NaN propagates through arithmetic
    x = nan + one
    call check_true('NaN + 1 = NaN', ieee_is_nan(x%limbs(1)), failures)
    x = one - nan
    call check_true('1 - NaN = NaN', ieee_is_nan(x%limbs(1)), failures)
    x = nan * one
    call check_true('NaN * 1 = NaN', ieee_is_nan(x%limbs(1)), failures)
    x = one / nan
    call check_true('1 / NaN = NaN', ieee_is_nan(x%limbs(1)), failures)
    x = sqrt(nan)
    call check_true('sqrt(NaN) = NaN', ieee_is_nan(x%limbs(1)), failures)

    ! NaN comparisons
    call check_true('NaN /= NaN', nan /= nan, failures)
    call check_true('not NaN < 1', .not. (nan < one), failures)
    call check_true('not NaN > 1', .not. (nan > one), failures)
    call check_true('not NaN == 1', .not. (nan == one), failures)
  end subroutine

  subroutine test_subnormal_boundary(failures)
    integer, intent(inout) :: failures
    type(float64x2) :: smallest, x
    real(dp) :: t

    ! Smallest normal double in DD form
    t = tiny(1.0_dp)
    smallest%limbs = [t, 0.0_dp]

    ! tiny - tiny = 0
    x = smallest - smallest
    call check_true('tiny - tiny = 0', x%limbs(1) == 0.0_dp, failures)

    ! tiny + tiny = 2*tiny (still normal)
    x = smallest + smallest
    call check_true('tiny + tiny = 2*tiny', x%limbs(1) == 2.0_dp * t, failures)

    ! tiny * 0.5 → subnormal in dp
    x = smallest * float64x2(0.5_dp)
    call check_true('tiny * 0.5 finite', ieee_is_finite(x%limbs(1)), failures)
    call check_true('tiny * 0.5 < tiny', x%limbs(1) < t, failures)

    ! sqrt(tiny) — should be ~1.49e-154
    x = sqrt(smallest)
    call check_true('sqrt(tiny) > 0', x%limbs(1) > 0.0_dp, failures)
    call check_true('sqrt(tiny) finite', ieee_is_finite(x%limbs(1)), failures)
  end subroutine

  subroutine test_huge_boundary(failures)
    integer, intent(inout) :: failures
    type(float64x2) :: big, x
    real(dp) :: h
    h = huge(1.0_dp)
    big%limbs = [h, 0.0_dp]

    ! huge + huge = +inf
    x = big + big
    call check_true('huge + huge = inf', .not. ieee_is_finite(x%limbs(1)) &
        .and. .not. ieee_is_nan(x%limbs(1)) .and. x%limbs(1) > 0.0_dp, failures)

    ! huge * 2 = +inf
    x = big * float64x2(2.0_dp)
    call check_true('huge * 2 = inf', .not. ieee_is_finite(x%limbs(1)) &
        .and. .not. ieee_is_nan(x%limbs(1)), failures)

    ! -huge - huge = -inf
    x = -big - big
    call check_true('-huge - huge = -inf', x%limbs(1) < -huge(1.0_dp) * 0.0_dp, failures)
    call check_true('-huge - huge non-finite', .not. ieee_is_finite(x%limbs(1)), failures)

    ! huge / huge = 1
    x = big / big
    call check_true('huge / huge = 1', abs(x%limbs(1) - 1.0_dp) <= epsilon(1.0_dp), failures)
  end subroutine

  ! ----------------------------------------------------------------
  ! Reduction variants: mask, dim, back
  ! ----------------------------------------------------------------

  subroutine test_reduction_mask(failures)
    integer, intent(inout) :: failures
    type(float64x2) :: a(5), s
    real(qp) :: qa(5)
    logical :: mask(5)
    integer :: i

    do i = 1, 5
      a(i) = float64x2(real(i, dp))
      qa(i) = real(i, qp)
    end do
    mask = [.true., .false., .true., .false., .true.]  ! 1, 3, 5

    s = sum(a, mask=mask)
    call check_real_close('sum mask', s, sum(qa, mask=mask), failures)

    s = product(a, mask=mask)
    call check_real_close('product mask', s, product(qa, mask=mask), failures, &
        scale=abs(product(qa, mask=mask)))

    s = maxval(a, mask=mask)
    call check_real_close('maxval mask', s, maxval(qa, mask=mask), failures)

    s = minval(a, mask=mask)
    call check_real_close('minval mask', s, minval(qa, mask=mask), failures)

    call check_true('maxloc mask', all(maxloc(a, mask=mask) == maxloc(qa, mask=mask)), failures)
    call check_true('minloc mask', all(minloc(a, mask=mask) == minloc(qa, mask=mask)), failures)
    call check_true('findloc mask', &
        all(findloc(a, a(3), mask=mask) == findloc(qa, qa(3), mask=mask)), failures)
  end subroutine

  subroutine test_reduction_dim(failures)
    integer, intent(inout) :: failures
    type(float64x2) :: m(2,3), v1(3), v2(2)
    real(qp) :: qm(2,3), qv1(3), qv2(2)
    integer :: i, j

    do j = 1, 3
      do i = 1, 2
        m(i,j) = float64x2(real(10*i + j, dp))
        qm(i,j) = real(10*i + j, qp)
      end do
    end do

    ! sum dim=1 → 1D array of length 3
    v1 = sum(m, dim=1)
    qv1 = sum(qm, dim=1)
    do j = 1, 3
      call check_real_close('sum dim=1', v1(j), qv1(j), failures, scale=abs(qv1(j)))
    end do

    ! sum dim=2 → 1D array of length 2
    v2 = sum(m, dim=2)
    qv2 = sum(qm, dim=2)
    do i = 1, 2
      call check_real_close('sum dim=2', v2(i), qv2(i), failures, scale=abs(qv2(i)))
    end do

    ! product dim=1
    v1 = product(m, dim=1)
    qv1 = product(qm, dim=1)
    do j = 1, 3
      call check_real_close('product dim=1', v1(j), qv1(j), failures, scale=abs(qv1(j)))
    end do

    ! maxval / minval dim=1
    v1 = maxval(m, dim=1)
    qv1 = maxval(qm, dim=1)
    do j = 1, 3
      call check_real_close('maxval dim=1', v1(j), qv1(j), failures, scale=abs(qv1(j)))
    end do
    v1 = minval(m, dim=1)
    qv1 = minval(qm, dim=1)
    do j = 1, 3
      call check_real_close('minval dim=1', v1(j), qv1(j), failures, scale=abs(qv1(j)))
    end do

    ! maxloc / minloc dim=2 → 1D integer of length 2
    block
      integer :: idx2(2), qidx2(2)
      idx2 = maxloc(m, dim=2)
      qidx2 = maxloc(qm, dim=2)
      call check_true('maxloc dim=2', all(idx2 == qidx2), failures)
      idx2 = minloc(m, dim=2)
      qidx2 = minloc(qm, dim=2)
      call check_true('minloc dim=2', all(idx2 == qidx2), failures)
    end block

    ! findloc dim=1
    block
      integer :: idx3(3), qidx3(3)
      idx3 = findloc(m, m(1,2), dim=1)
      qidx3 = findloc(qm, qm(1,2), dim=1)
      call check_true('findloc dim=1', all(idx3 == qidx3), failures)
    end block
  end subroutine

  ! ----------------------------------------------------------------
  ! Assignment-operator coverage (every direction supported)
  ! ----------------------------------------------------------------

  subroutine test_assignments_to_dd(failures)
    integer, intent(inout) :: failures
    type(float64x2) :: x

    ! dd <- real / int / complex of every kind we expose.
    x = 3.5_dp
    call check_real_close('dd <- dp', x, 3.5_qp, failures)
    x = 3.5_sp
    call check_real_close('dd <- sp', x, real(3.5_sp, qp), failures)
    x = 7
    call check_real_close('dd <- int', x, 7.0_qp, failures)
    x = 7_int8
    call check_real_close('dd <- int8', x, 7.0_qp, failures)
    x = 1234_int16
    call check_real_close('dd <- int16', x, 1234.0_qp, failures)
    x = 123456789012_int64
    call check_real_close('dd <- int64', x, 123456789012.0_qp, failures)
    x = (2.5_dp, -8.0_dp)            ! complex(dp); imag dropped
    call check_real_close('dd <- cdp', x, 2.5_qp, failures)
    x = (1.25_sp, 4.0_sp)            ! complex(sp); imag dropped
    call check_real_close('dd <- csp', x, real(1.25_sp, qp), failures)
  end subroutine

  subroutine test_assignments_to_cdd(failures)
    integer, intent(inout) :: failures
    type(complex64x2) :: z

    z = 3.5_dp
    call check_complex_close('cdd <- dp', z, cmplx(3.5_qp, 0.0_qp, qp), failures)
    z = 3.5_sp
    call check_complex_close('cdd <- sp', z, cmplx(real(3.5_sp, qp), 0.0_qp, qp), failures)
    z = 7
    call check_complex_close('cdd <- int', z, cmplx(7.0_qp, 0.0_qp, qp), failures)
    z = 7_int8
    call check_complex_close('cdd <- int8', z, cmplx(7.0_qp, 0.0_qp, qp), failures)
    z = 1234_int16
    call check_complex_close('cdd <- int16', z, cmplx(1234.0_qp, 0.0_qp, qp), failures)
    z = 123456789012_int64
    call check_complex_close('cdd <- int64', z, cmplx(123456789012.0_qp, 0.0_qp, qp), failures)
    z = (2.5_dp, -8.0_dp)
    call check_complex_close('cdd <- cdp', z, cmplx(2.5_qp, -8.0_qp, qp), failures)
    z = (1.25_sp, 4.0_sp)
    call check_complex_close('cdd <- csp', z, &
        cmplx(real(1.25_sp, qp), real(4.0_sp, qp), qp), failures)
  end subroutine

  subroutine test_assignments_from_dd(failures)
    integer, intent(inout) :: failures
    type(float64x2) :: x
    real(dp) :: d
    real(sp) :: s
    integer :: i
    integer(int8) :: j8
    integer(int16) :: j16
    integer(int64) :: j64
    complex(dp) :: cd
    complex(sp) :: cs

    x = 3.14159265358979_dp
    d = x
    call check_true('dp <- dd', d == x%limbs(1), failures)
    s = x
    call check_true('sp <- dd', s == real(x%limbs(1), sp), failures)
    cd = x
    call check_true('cdp <- dd real', real(cd, dp) == x%limbs(1), failures)
    call check_true('cdp <- dd imag', aimag(cd) == 0.0_dp, failures)
    cs = x
    call check_true('csp <- dd real', real(cs, sp) == real(x%limbs(1), sp), failures)
    call check_true('csp <- dd imag', aimag(cs) == 0.0_sp, failures)

    ! Integer truncation toward zero (matches Fortran int(...) semantics).
    x = 42.7_dp
    i = x; j8 = x; j16 = x; j64 = x
    call check_true('int  <- dd 42.7', i == 42, failures)
    call check_true('int8 <- dd 42.7', j8 == 42_int8, failures)
    call check_true('int16<- dd 42.7', j16 == 42_int16, failures)
    call check_true('int64<- dd 42.7', j64 == 42_int64, failures)
    x = -42.7_dp
    i = x
    call check_true('int  <- dd -42.7', i == -42, failures)
    ! Boundary: hi == integer, lo < 0 → truncate toward zero by 1
    x%limbs = [4.0_dp, -1.0e-30_dp]
    i = x
    call check_true('int  <- dd (4, -tiny)', i == 3, failures)
    x%limbs = [-4.0_dp, 1.0e-30_dp]
    i = x
    call check_true('int  <- dd (-4, tiny)', i == -3, failures)
  end subroutine

  subroutine test_constructors(failures)
    integer, intent(inout) :: failures
    type(float64x2) :: x
    type(complex64x2) :: z

    ! float64x2(...) — every supported single-arg form.
    x = float64x2(3.5_dp)
    call check_real_close('ctor float64x2 dp', x, 3.5_qp, failures)
    x = float64x2(3.5_sp)
    call check_real_close('ctor float64x2 sp', x, real(3.5_sp, qp), failures)
    x = float64x2(7)
    call check_real_close('ctor float64x2 int', x, 7.0_qp, failures)
    x = float64x2(7_int8)
    call check_real_close('ctor float64x2 int8', x, 7.0_qp, failures)
    x = float64x2(1234_int16)
    call check_real_close('ctor float64x2 int16', x, 1234.0_qp, failures)
    x = float64x2(123456789012_int64)
    call check_real_close('ctor float64x2 int64', x, 123456789012.0_qp, failures)
    x = float64x2('1.2345678901234567890123456789')
    call check_real_close('ctor float64x2 char', x, &
        1.2345678901234567890123456789_qp, failures)
    x = float64x2((2.5_dp, -8.0_dp))      ! cdp → dd takes the real part
    call check_real_close('ctor float64x2 cdp', x, 2.5_qp, failures)
    x = float64x2((1.25_sp, 4.0_sp))      ! csp → dd takes the real part
    call check_real_close('ctor float64x2 csp', x, real(1.25_sp, qp), failures)
    z = complex64x2(float64x2(1.0_dp), float64x2(2.0_dp))
    x = float64x2(z)                      ! cdd → dd takes the real part
    call check_real_close('ctor float64x2 cdd', x, 1.0_qp, failures)

    ! complex64x2(...) — single-arg forms.
    z = complex64x2(float64x2(1.5_dp))
    call check_complex_close('ctor cdd dd', z, cmplx(1.5_qp, 0.0_qp, qp), failures)
    z = complex64x2(3.5_dp)
    call check_complex_close('ctor cdd dp', z, cmplx(3.5_qp, 0.0_qp, qp), failures)
    z = complex64x2(3.5_sp)
    call check_complex_close('ctor cdd sp', z, cmplx(real(3.5_sp, qp), 0.0_qp, qp), failures)
    z = complex64x2(7)
    call check_complex_close('ctor cdd int', z, cmplx(7.0_qp, 0.0_qp, qp), failures)
    z = complex64x2(7_int8)
    call check_complex_close('ctor cdd int8', z, cmplx(7.0_qp, 0.0_qp, qp), failures)
    z = complex64x2(1234_int16)
    call check_complex_close('ctor cdd int16', z, cmplx(1234.0_qp, 0.0_qp, qp), failures)
    z = complex64x2(123456789012_int64)
    call check_complex_close('ctor cdd int64', z, &
        cmplx(123456789012.0_qp, 0.0_qp, qp), failures)
    z = complex64x2((2.5_dp, -8.0_dp))
    call check_complex_close('ctor cdd cdp', z, cmplx(2.5_qp, -8.0_qp, qp), failures)
    z = complex64x2((1.25_sp, 4.0_sp))
    call check_complex_close('ctor cdd csp', z, &
        cmplx(real(1.25_sp, qp), real(4.0_sp, qp), qp), failures)

    ! complex64x2(re, im) — matching-kind two-arg forms.
    z = complex64x2(float64x2(1.5_dp), float64x2(-2.5_dp))
    call check_complex_close('ctor cdd dd,dd', z, cmplx(1.5_qp, -2.5_qp, qp), failures)
    z = complex64x2(1.5_dp, -2.5_dp)
    call check_complex_close('ctor cdd dp,dp', z, cmplx(1.5_qp, -2.5_qp, qp), failures)
    z = complex64x2(1.5_sp, -2.5_sp)
    call check_complex_close('ctor cdd sp,sp', z, &
        cmplx(real(1.5_sp, qp), real(-2.5_sp, qp), qp), failures)
    z = complex64x2(10, 20)
    call check_complex_close('ctor cdd int,int', z, cmplx(10.0_qp, 20.0_qp, qp), failures)
    z = complex64x2(1_int8, 2_int8)
    call check_complex_close('ctor cdd i8,i8', z, cmplx(1.0_qp, 2.0_qp, qp), failures)
    z = complex64x2(100_int16, 200_int16)
    call check_complex_close('ctor cdd i16,i16', z, cmplx(100.0_qp, 200.0_qp, qp), failures)
    z = complex64x2(10_int64, 20_int64)
    call check_complex_close('ctor cdd i64,i64', z, cmplx(10.0_qp, 20.0_qp, qp), failures)

    ! complex64x2 mixed-kind two-arg forms.
    z = complex64x2(float64x2(1.5_dp), 2.5_dp)
    call check_complex_close('ctor cdd dd,dp', z, cmplx(1.5_qp, 2.5_qp, qp), failures)
    z = complex64x2(1.5_dp, float64x2(2.5_dp))
    call check_complex_close('ctor cdd dp,dd', z, cmplx(1.5_qp, 2.5_qp, qp), failures)
  end subroutine

  subroutine test_assignments_from_cdd(failures)
    integer, intent(inout) :: failures
    type(complex64x2) :: z
    real(dp) :: d
    real(sp) :: s
    integer :: i
    integer(int8) :: j8
    integer(int16) :: j16
    integer(int64) :: j64
    complex(dp) :: cd
    complex(sp) :: cs

    z = (2.75_dp, -1.5_dp)
    d = z
    call check_true('dp <- cdd', d == z%re%limbs(1), failures)
    s = z
    call check_true('sp <- cdd', s == real(z%re%limbs(1), sp), failures)
    cd = z
    call check_true('cdp <- cdd re', real(cd, dp) == z%re%limbs(1), failures)
    call check_true('cdp <- cdd im', aimag(cd) == z%im%limbs(1), failures)
    cs = z
    call check_true('csp <- cdd re', real(cs, sp) == real(z%re%limbs(1), sp), failures)
    call check_true('csp <- cdd im', aimag(cs) == real(z%im%limbs(1), sp), failures)

    z = (10.99_dp, 0.0_dp)
    i = z; j8 = z; j16 = z; j64 = z
    call check_true('int  <- cdd 10.99', i == 10, failures)
    call check_true('int8 <- cdd 10.99', j8 == 10_int8, failures)
    call check_true('int16<- cdd 10.99', j16 == 10_int16, failures)
    call check_true('int64<- cdd 10.99', j64 == 10_int64, failures)
  end subroutine

  subroutine test_loc_back(failures)
    integer, intent(inout) :: failures
    type(float64x2) :: a(5)
    real(qp) :: qa(5)
    integer :: i

    ! Build an array with a duplicated max so back= matters
    a = [float64x2(1.0_dp), float64x2(3.0_dp), float64x2(2.0_dp), &
         float64x2(3.0_dp), float64x2(0.0_dp)]
    do i = 1, 5
      qa(i) = real(a(i)%limbs(1), qp)
    end do

    call check_true('maxloc forward', all(maxloc(a) == maxloc(qa)), failures)
    call check_true('maxloc back', &
        all(maxloc(a, back=.true.) == maxloc(qa, back=.true.)), failures)
    call check_true('minloc forward', all(minloc(a) == minloc(qa)), failures)
    call check_true('minloc back', &
        all(minloc(a, back=.true.) == minloc(qa, back=.true.)), failures)

    call check_true('findloc forward', &
        all(findloc(a, float64x2(3.0_dp)) == findloc(qa, 3.0_qp)), failures)
    call check_true('findloc back', &
        all(findloc(a, float64x2(3.0_dp), back=.true.) == findloc(qa, 3.0_qp, back=.true.)), &
        failures)
  end subroutine

  subroutine test_ulp_boundary(failures)
    integer, intent(inout) :: failures
    type(float64x2) :: one, eps_dd, x
    real(dp) :: e
    one = float64x2(1.0_dp)
    e = epsilon(1.0_dp)
    eps_dd%limbs = [e, 0.0_dp]

    ! 1 + epsilon - 1 = epsilon (no precision loss in DD)
    x = one + eps_dd - one
    call check_true('(1+eps)-1 = eps', x%limbs(1) == e .and. x%limbs(2) == 0.0_dp, failures)

    ! Lo limb representable: 1 + eps^2/2 should be exact in DD
    x = one
    x%limbs(2) = 0.5_dp * e * e
    call check_true('lo limb stays', x%limbs(2) == 0.5_dp * e * e, failures)

    ! 1 + ulp(lo) should round into lo limb of DD
    x = one + float64x2(0.5_dp * e * e)
    call check_true('1 + 0.5 eps^2 finite', ieee_is_finite(x%limbs(1)) &
        .and. ieee_is_finite(x%limbs(2)), failures)
  end subroutine

end program test_multifloats_precision
