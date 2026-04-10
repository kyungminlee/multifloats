program test_multifloats_precision
  use multifloats
  implicit none

  integer, parameter :: qp = 16
  integer, parameter :: dp = 8
  integer, parameter :: sp = 4
  real(qp), parameter :: tol = 1.0e5_qp * epsilon(1.0_qp)
  integer :: failures

  failures = 0

  call test_scalar_constructors(failures)
  call test_scalar_arithmetic(failures)
  call test_unary_math(failures)
  call test_complex_arithmetic(failures)
  call test_array_ops(failures)
  call test_rounding_and_random(failures)

  if (failures /= 0) then
    write(*,'(a,i0)') 'FAIL: ', failures
    error stop 1
  end if

  write(*,'(a)') 'PASS'

contains

  pure elemental function mf_to_qp(x) result(y)
    type(float64x2), intent(in) :: x
    real(qp) :: y
    y = real(x%limbs(1), qp) + real(x%limbs(2), qp)
  end function

  pure elemental function cx_to_cqp(x) result(y)
    type(complex128x2), intent(in) :: x
    complex(qp) :: y
    y = cmplx(mf_to_qp(x%re), mf_to_qp(x%im), kind=qp)
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
    real(qp) :: err, bound, mag
    if (present(scale)) then
      mag = scale
    else
      mag = max(1.0_qp, abs(expect))
    end if
    bound = tol * mag
    err = abs(mf_to_qp(got) - expect)
    if (err > bound) then
      failures = failures + 1
      write(*,'(a,1x,es24.16,1x,es24.16,1x,es24.16)') 'FAIL '//trim(label), err, bound, expect
    end if
  end subroutine

  subroutine check_complex_close(label, got, expect, failures, scale)
    character(*), intent(in) :: label
    type(complex128x2), intent(in) :: got
    complex(qp), intent(in) :: expect
    integer, intent(inout) :: failures
    real(qp), intent(in), optional :: scale
    real(qp) :: err, bound, mag
    if (present(scale)) then
      mag = scale
    else
      mag = max(1.0_qp, abs(expect))
    end if
    bound = tol * mag
    err = abs(cx_to_cqp(got) - expect)
    if (err > bound) then
      failures = failures + 1
      write(*,'(a,1x,es24.16,1x,es24.16)') 'FAIL '//trim(label), err, bound
    end if
  end subroutine

  subroutine test_scalar_constructors(failures)
    integer, intent(inout) :: failures
    type(float64x2) :: x
    type(complex128x2) :: z

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

    z = complex128x2(float64x2(1.5_dp), float64x2(-0.25_dp))
    call check_complex_close('ctor complex mf', z, cmplx(1.5_qp, -0.25_qp, qp), failures)

    z = complex128x2(3.0_dp, -2.0_dp)
    call check_complex_close('ctor complex dp', z, cmplx(3.0_qp, -2.0_qp, qp), failures)
  end subroutine

  subroutine test_scalar_arithmetic(failures)
    integer, intent(inout) :: failures
    type(float64x2) :: x, y
    real(qp) :: qx, qy

    x%limbs = [1.0_dp, 2.0e-30_dp]
    y%limbs = [2.0_dp, -3.0e-30_dp]
    qx = mf_to_qp(x)
    qy = mf_to_qp(y)

    call check_real_close('add mf-mf', x + y, qx + qy, failures)
    call check_real_close('sub mf-mf', x - y, qx - qy, failures)
    call check_real_close('mul mf-mf', x * y, qx * qy, failures)
    call check_real_close('div mf-mf', x / y, qx / qy, failures)

    call check_real_close('add mf-dp', x + 2.0_dp, qx + 2.0_qp, failures)
    call check_real_close('add dp-mf', 2.0_dp + x, 2.0_qp + qx, failures)
    call check_real_close('sub mf-int', x - 2, qx - 2.0_qp, failures)
    call check_real_close('sub int-mf', 2 - x, 2.0_qp - qx, failures)
    call check_real_close('mul mf-sp', x * real(0.5_sp, sp), qx * 0.5_qp, failures)
    call check_real_close('mul sp-mf', real(0.5_sp, sp) * x, 0.5_qp * qx, failures)
    call check_real_close('div mf-dp', x / 2.0_dp, qx / 2.0_qp, failures)
    call check_real_close('div dp-mf', 2.0_dp / x, 2.0_qp / qx, failures)
    call check_real_close('pow int', x ** 3, qx ** 3, failures, scale=abs(qx ** 3))

    call check_true('comparison lt', x < y, failures)
    call check_true('comparison gt', y > x, failures)
    call check_true('comparison eq', float64x2(2) == 2, failures)
    call check_true('comparison ne', x /= y, failures)
  end subroutine

  subroutine test_unary_math(failures)
    integer, intent(inout) :: failures
    type(float64x2) :: x
    real(qp) :: qx

    x%limbs = [1.25_dp, 1.0e-30_dp]
    qx = mf_to_qp(x)

    call check_real_close('sqrt', sqrt(x), sqrt(qx), failures)
    call check_real_close('exp', exp(x), exp(qx), failures, scale=abs(exp(qx)))
    call check_real_close('log', log(x), log(qx), failures)
    call check_real_close('sin', sin(x), sin(qx), failures)
    call check_real_close('cos', cos(x), cos(qx), failures)
    call check_real_close('abs', abs(-x), abs(-qx), failures)
  end subroutine

  subroutine test_complex_arithmetic(failures)
    integer, intent(inout) :: failures
    type(complex128x2) :: z1, z2
    complex(qp) :: q1, q2

    z1 = complex128x2(float64x2(1.25_dp), float64x2(-0.5_dp))
    z2 = complex128x2(0.75_dp, 0.25_dp)
    q1 = cx_to_cqp(z1)
    q2 = cx_to_cqp(z2)

    call check_complex_close('cx add', z1 + z2, q1 + q2, failures)
    call check_complex_close('cx sub mixed', z1 - cmplx(0.5_dp, -0.125_dp, dp), q1 - cmplx(0.5_qp, -0.125_qp, qp), failures)
    call check_complex_close('cx mul', z1 * z2, q1 * q2, failures)
    call check_complex_close('cx div', z1 / z2, q1 / q2, failures)
    call check_complex_close('cx sqrt', sqrt(z1), sqrt(q1), failures)
    call check_complex_close('cx exp', exp(z1), exp(q1), failures, scale=abs(exp(q1)))

    call check_real_close('cx abs', abs(z1), abs(q1), failures)
    call check_real_close('cx aimag', aimag(z1), aimag(q1), failures)
  end subroutine

  subroutine test_array_ops(failures)
    integer, intent(inout) :: failures
    type(float64x2) :: a(2), m1(2,2), m2(2,2), mv(2), mm(2,2)
    type(complex128x2) :: cz(2)
    real(qp) :: qa(2), qm1(2,2), qm2(2,2), qv(2), qmm(2,2)
    complex(qp) :: qcz(2)

    a(1)%limbs = [1.0_dp, 1.0e-30_dp]
    a(2)%limbs = [2.0_dp, -2.0e-30_dp]
    qa = [mf_to_qp(a(1)), mf_to_qp(a(2))]

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

    cz = [complex128x2(1.0_dp, -0.5_dp), complex128x2(0.25_dp, 0.75_dp)]
    qcz = [cx_to_cqp(cz(1)), cx_to_cqp(cz(2))]
    call check_complex_close('complex sum', sum(cz), sum(qcz), failures)
    call check_complex_close('complex dot_product', dot_product(cz, cz), dot_product(qcz, qcz), failures)
  end subroutine

  subroutine test_rounding_and_random(failures)
    integer, intent(inout) :: failures
    type(float64x2) :: x, r

    x%limbs = [3.5_dp, 2.0e-30_dp]
    call check_true('dble', abs(dble(x) - real(mf_to_qp(x), dp)) <= epsilon(1.0_dp), failures)
    call check_true('int', int(x) == int(mf_to_qp(x)), failures)
    call check_true('nint', nint(x) == nint(mf_to_qp(x)), failures)
    call check_true('floor', floor(x) == floor(mf_to_qp(x)), failures)
    call check_true('ceiling', ceiling(x) == ceiling(mf_to_qp(x)), failures)

    call random_number(r)
    call check_true('random lower', mf_to_qp(r) >= 0.0_qp, failures)
    call check_true('random upper', mf_to_qp(r) < 1.0_qp, failures)
  end subroutine

end program test_multifloats_precision
