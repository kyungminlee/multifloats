program bench_c_wrappers
  ! Three-way benchmark: quad precision (qp) vs Fortran multifloats (mf)
  ! vs C-ABI wrappers around C++ multifloats (cmf).
  !
  ! The C wrappers return float64x2_t (16 bytes) by value via the platform C
  ! ABI, which on ARM64 uses d0/d1 registers — no hidden pointer. This
  ! isolates the cost of gfortran's derived-type ABI vs the C register-
  ! return convention for the same underlying DD arithmetic.
  use multifloats
  use, intrinsic :: iso_fortran_env, only: int64, real64
  use, intrinsic :: iso_c_binding,   only: c_double, c_int
  implicit none

  integer, parameter :: qp = 16
  integer, parameter :: dp = 8

  ! C-interoperable mirror of float64x2_t.
  type, bind(c) :: float64x2_t
    real(c_double) :: hi, lo
  end type float64x2_t

  ! ================================================================
  ! C wrapper interfaces (pass & return by value)
  ! ================================================================
  interface
    ! Arithmetic
    type(float64x2_t) function c_dd_add(a, b) bind(c, name='adddd')
      import :: float64x2_t
      type(float64x2_t), value :: a, b
    end function
    type(float64x2_t) function c_dd_sub(a, b) bind(c, name='subdd')
      import :: float64x2_t
      type(float64x2_t), value :: a, b
    end function
    type(float64x2_t) function c_dd_mul(a, b) bind(c, name='muldd')
      import :: float64x2_t
      type(float64x2_t), value :: a, b
    end function
    type(float64x2_t) function c_dd_div(a, b) bind(c, name='divdd')
      import :: float64x2_t
      type(float64x2_t), value :: a, b
    end function
    ! Unary
    type(float64x2_t) function c_dd_neg(a) bind(c, name='negdd')
      import :: float64x2_t
      type(float64x2_t), value :: a
    end function
    type(float64x2_t) function c_dd_abs(a) bind(c, name='fabsdd')
      import :: float64x2_t
      type(float64x2_t), value :: a
    end function
    type(float64x2_t) function c_dd_sqrt(a) bind(c, name='sqrtdd')
      import :: float64x2_t
      type(float64x2_t), value :: a
    end function
    ! Binary
    type(float64x2_t) function c_dd_fmin(a, b) bind(c, name='fmindd')
      import :: float64x2_t
      type(float64x2_t), value :: a, b
    end function
    type(float64x2_t) function c_dd_fmax(a, b) bind(c, name='fmaxdd')
      import :: float64x2_t
      type(float64x2_t), value :: a, b
    end function
    type(float64x2_t) function c_dd_hypot(a, b) bind(c, name='hypotdd')
      import :: float64x2_t
      type(float64x2_t), value :: a, b
    end function
    type(float64x2_t) function c_dd_fdim(a, b) bind(c, name='fdimdd')
      import :: float64x2_t
      type(float64x2_t), value :: a, b
    end function
    type(float64x2_t) function c_dd_fmod(a, b) bind(c, name='fmoddd')
      import :: float64x2_t
      type(float64x2_t), value :: a, b
    end function
    type(float64x2_t) function c_dd_copysign(a, b) bind(c, name='copysigndd')
      import :: float64x2_t
      type(float64x2_t), value :: a, b
    end function
    ! Transcendental
    type(float64x2_t) function c_dd_exp(a) bind(c, name='expdd')
      import :: float64x2_t
      type(float64x2_t), value :: a
    end function
    type(float64x2_t) function c_dd_log(a) bind(c, name='logdd')
      import :: float64x2_t
      type(float64x2_t), value :: a
    end function
    type(float64x2_t) function c_dd_log10(a) bind(c, name='log10dd')
      import :: float64x2_t
      type(float64x2_t), value :: a
    end function
    type(float64x2_t) function c_dd_pow(a, b) bind(c, name='powdd')
      import :: float64x2_t
      type(float64x2_t), value :: a, b
    end function
    type(float64x2_t) function c_dd_sin(a) bind(c, name='sindd')
      import :: float64x2_t
      type(float64x2_t), value :: a
    end function
    type(float64x2_t) function c_dd_cos(a) bind(c, name='cosdd')
      import :: float64x2_t
      type(float64x2_t), value :: a
    end function
    type(float64x2_t) function c_dd_tan(a) bind(c, name='tandd')
      import :: float64x2_t
      type(float64x2_t), value :: a
    end function
    type(float64x2_t) function c_dd_asin(a) bind(c, name='asindd')
      import :: float64x2_t
      type(float64x2_t), value :: a
    end function
    type(float64x2_t) function c_dd_acos(a) bind(c, name='acosdd')
      import :: float64x2_t
      type(float64x2_t), value :: a
    end function
    type(float64x2_t) function c_dd_atan(a) bind(c, name='atandd')
      import :: float64x2_t
      type(float64x2_t), value :: a
    end function
    type(float64x2_t) function c_dd_atan2(a, b) bind(c, name='atan2dd')
      import :: float64x2_t
      type(float64x2_t), value :: a, b
    end function
    type(float64x2_t) function c_dd_sinh(a) bind(c, name='sinhdd')
      import :: float64x2_t
      type(float64x2_t), value :: a
    end function
    type(float64x2_t) function c_dd_cosh(a) bind(c, name='coshdd')
      import :: float64x2_t
      type(float64x2_t), value :: a
    end function
    type(float64x2_t) function c_dd_tanh(a) bind(c, name='tanhdd')
      import :: float64x2_t
      type(float64x2_t), value :: a
    end function
    type(float64x2_t) function c_dd_asinh(a) bind(c, name='asinhdd')
      import :: float64x2_t
      type(float64x2_t), value :: a
    end function
    type(float64x2_t) function c_dd_acosh(a) bind(c, name='acoshdd')
      import :: float64x2_t
      type(float64x2_t), value :: a
    end function
    type(float64x2_t) function c_dd_atanh(a) bind(c, name='atanhdd')
      import :: float64x2_t
      type(float64x2_t), value :: a
    end function
    type(float64x2_t) function c_dd_erf(a) bind(c, name='erfdd')
      import :: float64x2_t
      type(float64x2_t), value :: a
    end function
    type(float64x2_t) function c_dd_erfc(a) bind(c, name='erfcdd')
      import :: float64x2_t
      type(float64x2_t), value :: a
    end function
    ! Comparison
    integer(c_int) function c_dd_lt(a, b) bind(c, name='ltdd')
      import :: float64x2_t, c_int
      type(float64x2_t), value :: a, b
    end function
    integer(c_int) function c_dd_le(a, b) bind(c, name='ledd')
      import :: float64x2_t, c_int
      type(float64x2_t), value :: a, b
    end function
    integer(c_int) function c_dd_eq(a, b) bind(c, name='eqdd')
      import :: float64x2_t, c_int
      type(float64x2_t), value :: a, b
    end function
  end interface

  ! ================================================================
  ! Data
  ! ================================================================
  integer, parameter :: N = 1024
  integer, parameter :: REPS_FAST = 400
  integer, parameter :: REPS_TRIG = 40
  integer, parameter :: REPS_VERY_SLOW = 4

  real(qp), allocatable :: q1(:), q2(:), qpos(:), qsmall(:), qbnd(:), qres(:)
  type(float64x2), allocatable :: f1(:), f2(:), fpos(:), fsmall(:), fbnd(:), fres(:)
  type(float64x2_t), allocatable :: c1(:), c2(:), cpos(:), csmall(:), cbnd(:), cres(:)

  real(qp) :: q_sink = 0.0_qp
  real(dp) :: f_sink = 0.0_dp
  real(dp) :: c_sink = 0.0_dp

  call init_data()

  print '(a)', "================================================================"
  print '(a)', " Three-way benchmark: qp vs Fortran-mf vs C-wrapped-mf"
  print '(a,i0)', " N=", N
  print '(a)', "================================================================"
  print '(a)', ""
  print '(a)', &
    " op                      n_ops     qp [s]     mf [s]    cmf [s]  qp/mf   qp/cmf  mf/cmf"
  print '(a)', &
    " --------------------------------------------------------------------------------------"

  call bench_arith()
  call bench_unary()
  call bench_binary()
  call bench_exp_log()
  call bench_trig()
  call bench_inv_trig()
  call bench_hyperbolic()
  call bench_inv_hyperbolic()
  call bench_erf()

  print '(a)', &
    " --------------------------------------------------------------------------------------"
  print '(a,es11.3,a,es11.3,a,es11.3)', &
    " sinks: qp=", real(q_sink, dp), " mf=", f_sink, " cmf=", c_sink

contains

  ! ================================================================
  ! Plumbing
  ! ================================================================

  type(float64x2) function to_f64x2(q)
    real(qp), intent(in) :: q
    real(dp) :: h, l, s, b
    h = real(q, dp)
    l = real(q - real(h, qp), dp)
    s = h + l
    b = s - h
    to_f64x2%limbs(1) = s
    to_f64x2%limbs(2) = l - b
  end function

  type(float64x2_t) function mf_to_dd(f)
    type(float64x2), intent(in) :: f
    mf_to_dd%hi = f%limbs(1)
    mf_to_dd%lo = f%limbs(2)
  end function

  subroutine init_data()
    integer :: i
    real(dp) :: r(3)
    allocate(q1(N), q2(N), qpos(N), qsmall(N), qbnd(N), qres(N))
    allocate(f1(N), f2(N), fpos(N), fsmall(N), fbnd(N), fres(N))
    allocate(c1(N), c2(N), cpos(N), csmall(N), cbnd(N), cres(N))
    call random_seed()
    do i = 1, N
      call random_number(r)
      q1(i) = real((r(1) - 0.5_dp) * 8.0_dp + sign(0.25_dp, r(1) - 0.5_dp), qp)
      q2(i) = real((r(2) - 0.5_dp) * 8.0_dp + sign(0.25_dp, r(2) - 0.5_dp), qp)
      qpos(i) = real(0.1_dp + 9.9_dp * r(1), qp)
      qsmall(i) = real((r(2) - 0.5_dp) * 6.0_dp, qp)
      qbnd(i) = real((r(3) - 0.5_dp) * 1.8_dp, qp)
    end do
    do i = 1, N
      f1(i) = to_f64x2(q1(i));     c1(i) = mf_to_dd(f1(i))
      f2(i) = to_f64x2(q2(i));     c2(i) = mf_to_dd(f2(i))
      fpos(i) = to_f64x2(qpos(i)); cpos(i) = mf_to_dd(fpos(i))
      fsmall(i) = to_f64x2(qsmall(i)); csmall(i) = mf_to_dd(fsmall(i))
      fbnd(i) = to_f64x2(qbnd(i)); cbnd(i) = mf_to_dd(fbnd(i))
    end do
  end subroutine

  subroutine tick(t)
    integer(int64), intent(out) :: t
    call system_clock(t)
  end subroutine

  real(real64) function elapsed(t0)
    integer(int64), intent(in) :: t0
    integer(int64) :: t1, rate
    call system_clock(t1, rate)
    elapsed = real(t1 - t0, real64) / real(rate, real64)
  end function

  ! NOINLINE drains — one per leg.
  subroutine qfeed(q_in)
    !GCC$ ATTRIBUTES NOINLINE :: qfeed
    real(qp), intent(inout) :: q_in(:)
    real(dp) :: s
    integer :: i
    s = 0.0_dp
    do i = 1, size(qres)
      s = s + real(qres(i), dp)
    end do
    q_in(1) = q_in(1) + real(s, qp) * 1.0e-30_qp
    q_sink = q_sink + real(s, qp)
  end subroutine

  subroutine ffeed(f_in)
    !GCC$ ATTRIBUTES NOINLINE :: ffeed
    type(float64x2), intent(inout) :: f_in(:)
    real(dp) :: s
    integer :: i
    s = 0.0_dp
    do i = 1, size(fres)
      s = s + fres(i)%limbs(1)
    end do
    f_in(1)%limbs(1) = f_in(1)%limbs(1) + s * 1.0e-30_dp
    f_sink = f_sink + s
  end subroutine

  subroutine cfeed(c_in)
    !GCC$ ATTRIBUTES NOINLINE :: cfeed
    type(float64x2_t), intent(inout) :: c_in(:)
    real(dp) :: s
    integer :: i
    s = 0.0_dp
    do i = 1, size(cres)
      s = s + cres(i)%hi
    end do
    c_in(1)%hi = c_in(1)%hi + s * 1.0e-30_dp
    c_sink = c_sink + s
  end subroutine

  subroutine report3(name, n_ops, tq, tf, tc)
    character(*), intent(in) :: name
    integer(int64), intent(in) :: n_ops
    real(real64), intent(in) :: tq, tf, tc
    real(real64) :: sq, sc, mc
    sq = merge(tq / tf, 0.0_real64, tf > 0)
    sc = merge(tq / tc, 0.0_real64, tc > 0)
    mc = merge(tf / tc, 0.0_real64, tc > 0)
    write(*, '(1x,a22,1x,i10,1x,f9.4,1x,f9.4,1x,f9.4,1x,f7.2,"x",1x,f7.2,"x",1x,f7.2,"x")') &
        name, n_ops, tq, tf, tc, sq, sc, mc
  end subroutine

  ! ================================================================
  ! Benchmark sections
  ! ================================================================

  subroutine bench_arith()
    integer :: i, r
    integer(int64) :: t0, n_ops
    real(real64) :: tq, tf, tc
    n_ops = int(N, int64) * int(REPS_FAST, int64)

    ! add
    call tick(t0)
    do r = 1, REPS_FAST; do i = 1, N; qres(i) = q1(i) + q2(i); end do; call qfeed(q1); end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST; do i = 1, N; fres(i) = f1(i) + f2(i); end do; call ffeed(f1); end do
    tf = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST; do i = 1, N; cres(i) = c_dd_add(c1(i), c2(i)); end do; call cfeed(c1); end do
    tc = elapsed(t0)
    call report3("add", n_ops, tq, tf, tc)

    ! sub
    call tick(t0)
    do r = 1, REPS_FAST; do i = 1, N; qres(i) = q1(i) - q2(i); end do; call qfeed(q1); end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST; do i = 1, N; fres(i) = f1(i) - f2(i); end do; call ffeed(f1); end do
    tf = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST; do i = 1, N; cres(i) = c_dd_sub(c1(i), c2(i)); end do; call cfeed(c1); end do
    tc = elapsed(t0)
    call report3("sub", n_ops, tq, tf, tc)

    ! mul
    call tick(t0)
    do r = 1, REPS_FAST; do i = 1, N; qres(i) = q1(i) * q2(i); end do; call qfeed(q1); end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST; do i = 1, N; fres(i) = f1(i) * f2(i); end do; call ffeed(f1); end do
    tf = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST; do i = 1, N; cres(i) = c_dd_mul(c1(i), c2(i)); end do; call cfeed(c1); end do
    tc = elapsed(t0)
    call report3("mul", n_ops, tq, tf, tc)

    ! div
    call tick(t0)
    do r = 1, REPS_FAST; do i = 1, N; qres(i) = q1(i) / q2(i); end do; call qfeed(q1); end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST; do i = 1, N; fres(i) = f1(i) / f2(i); end do; call ffeed(f1); end do
    tf = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST; do i = 1, N; cres(i) = c_dd_div(c1(i), c2(i)); end do; call cfeed(c1); end do
    tc = elapsed(t0)
    call report3("div", n_ops, tq, tf, tc)
  end subroutine

  subroutine bench_unary()
    integer :: i, r
    integer(int64) :: t0, n_ops
    real(real64) :: tq, tf, tc
    n_ops = int(N, int64) * int(REPS_FAST, int64)

    ! sqrt
    call tick(t0)
    do r = 1, REPS_FAST; do i = 1, N; qres(i) = sqrt(qpos(i)); end do; call qfeed(qpos); end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST; do i = 1, N; fres(i) = sqrt(fpos(i)); end do; call ffeed(fpos); end do
    tf = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST; do i = 1, N; cres(i) = c_dd_sqrt(cpos(i)); end do; call cfeed(cpos); end do
    tc = elapsed(t0)
    call report3("sqrt", n_ops, tq, tf, tc)

    ! abs
    call tick(t0)
    do r = 1, REPS_FAST; do i = 1, N; qres(i) = abs(q1(i)); end do; call qfeed(q1); end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST; do i = 1, N; fres(i) = abs(f1(i)); end do; call ffeed(f1); end do
    tf = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST; do i = 1, N; cres(i) = c_dd_abs(c1(i)); end do; call cfeed(c1); end do
    tc = elapsed(t0)
    call report3("abs", n_ops, tq, tf, tc)

    ! neg
    call tick(t0)
    do r = 1, REPS_FAST; do i = 1, N; qres(i) = -q1(i); end do; call qfeed(q1); end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST; do i = 1, N; fres(i) = -f1(i); end do; call ffeed(f1); end do
    tf = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST; do i = 1, N; cres(i) = c_dd_neg(c1(i)); end do; call cfeed(c1); end do
    tc = elapsed(t0)
    call report3("neg", n_ops, tq, tf, tc)
  end subroutine

  subroutine bench_binary()
    integer :: i, r
    integer(int64) :: t0, n_ops
    real(real64) :: tq, tf, tc
    n_ops = int(N, int64) * int(REPS_FAST, int64)

    call tick(t0)
    do r = 1, REPS_FAST; do i = 1, N; qres(i) = min(q1(i),q2(i)); end do; call qfeed(q1); end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST; do i = 1, N; fres(i) = min(f1(i),f2(i)); end do; call ffeed(f1); end do
    tf = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST; do i = 1, N; cres(i) = c_dd_fmin(c1(i),c2(i)); end do; call cfeed(c1); end do
    tc = elapsed(t0)
    call report3("fmin", n_ops, tq, tf, tc)

    call tick(t0)
    do r = 1, REPS_FAST; do i = 1, N; qres(i) = max(q1(i),q2(i)); end do; call qfeed(q1); end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST; do i = 1, N; fres(i) = max(f1(i),f2(i)); end do; call ffeed(f1); end do
    tf = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST; do i = 1, N; cres(i) = c_dd_fmax(c1(i),c2(i)); end do; call cfeed(c1); end do
    tc = elapsed(t0)
    call report3("fmax", n_ops, tq, tf, tc)

    call tick(t0)
    do r = 1, REPS_FAST; do i = 1, N; qres(i) = hypot(q1(i),q2(i)); end do; call qfeed(q1); end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST; do i = 1, N; fres(i) = hypot(f1(i),f2(i)); end do; call ffeed(f1); end do
    tf = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST; do i = 1, N; cres(i) = c_dd_hypot(c1(i),c2(i)); end do; call cfeed(c1); end do
    tc = elapsed(t0)
    call report3("hypot", n_ops, tq, tf, tc)

    call tick(t0)
    do r = 1, REPS_FAST; do i = 1, N; qres(i) = dim(q1(i),q2(i)); end do; call qfeed(q1); end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST; do i = 1, N; fres(i) = dim(f1(i),f2(i)); end do; call ffeed(f1); end do
    tf = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST; do i = 1, N; cres(i) = c_dd_fdim(c1(i),c2(i)); end do; call cfeed(c1); end do
    tc = elapsed(t0)
    call report3("fdim", n_ops, tq, tf, tc)

    call tick(t0)
    do r = 1, REPS_FAST; do i = 1, N; qres(i) = mod(q1(i),q2(i)); end do; call qfeed(q1); end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST; do i = 1, N; fres(i) = mod(f1(i),f2(i)); end do; call ffeed(f1); end do
    tf = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST; do i = 1, N; cres(i) = c_dd_fmod(c1(i),c2(i)); end do; call cfeed(c1); end do
    tc = elapsed(t0)
    call report3("fmod", n_ops, tq, tf, tc)
  end subroutine

  subroutine bench_exp_log()
    integer :: i, r
    integer(int64) :: t0, n_ops
    real(real64) :: tq, tf, tc
    n_ops = int(N, int64) * int(REPS_TRIG, int64)

    call tick(t0)
    do r = 1, REPS_TRIG; do i = 1, N; qres(i) = exp(qsmall(i)); end do; call qfeed(qsmall); end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_TRIG; do i = 1, N; fres(i) = exp(fsmall(i)); end do; call ffeed(fsmall); end do
    tf = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_TRIG; do i = 1, N; cres(i) = c_dd_exp(csmall(i)); end do; call cfeed(csmall); end do
    tc = elapsed(t0)
    call report3("exp", n_ops, tq, tf, tc)

    call tick(t0)
    do r = 1, REPS_TRIG; do i = 1, N; qres(i) = log(qpos(i)); end do; call qfeed(qpos); end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_TRIG; do i = 1, N; fres(i) = log(fpos(i)); end do; call ffeed(fpos); end do
    tf = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_TRIG; do i = 1, N; cres(i) = c_dd_log(cpos(i)); end do; call cfeed(cpos); end do
    tc = elapsed(t0)
    call report3("log", n_ops, tq, tf, tc)

    call tick(t0)
    do r = 1, REPS_TRIG; do i = 1, N; qres(i) = log10(qpos(i)); end do; call qfeed(qpos); end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_TRIG; do i = 1, N; fres(i) = log10(fpos(i)); end do; call ffeed(fpos); end do
    tf = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_TRIG; do i = 1, N; cres(i) = c_dd_log10(cpos(i)); end do; call cfeed(cpos); end do
    tc = elapsed(t0)
    call report3("log10", n_ops, tq, tf, tc)

    call tick(t0)
    do r = 1, REPS_TRIG; do i = 1, N; qres(i) = qpos(i)**qbnd(i); end do; call qfeed(qpos); end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_TRIG; do i = 1, N; fres(i) = fpos(i)**fbnd(i); end do; call ffeed(fpos); end do
    tf = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_TRIG; do i = 1, N; cres(i) = c_dd_pow(cpos(i),cbnd(i)); end do; call cfeed(cpos); end do
    tc = elapsed(t0)
    call report3("pow", n_ops, tq, tf, tc)
  end subroutine

  subroutine bench_trig()
    integer :: i, r
    integer(int64) :: t0, n_ops
    real(real64) :: tq, tf, tc
    n_ops = int(N, int64) * int(REPS_TRIG, int64)

    call tick(t0)
    do r=1,REPS_TRIG; do i=1,N; qres(i)=sin(qsmall(i)); end do; call qfeed(qsmall); end do
    tq = elapsed(t0)
    call tick(t0)
    do r=1,REPS_TRIG; do i=1,N; fres(i)=sin(fsmall(i)); end do; call ffeed(fsmall); end do
    tf = elapsed(t0)
    call tick(t0)
    do r=1,REPS_TRIG; do i=1,N; cres(i)=c_dd_sin(csmall(i)); end do; call cfeed(csmall); end do
    tc = elapsed(t0)
    call report3("sin", n_ops, tq, tf, tc)

    call tick(t0)
    do r=1,REPS_TRIG; do i=1,N; qres(i)=cos(qsmall(i)); end do; call qfeed(qsmall); end do
    tq = elapsed(t0)
    call tick(t0)
    do r=1,REPS_TRIG; do i=1,N; fres(i)=cos(fsmall(i)); end do; call ffeed(fsmall); end do
    tf = elapsed(t0)
    call tick(t0)
    do r=1,REPS_TRIG; do i=1,N; cres(i)=c_dd_cos(csmall(i)); end do; call cfeed(csmall); end do
    tc = elapsed(t0)
    call report3("cos", n_ops, tq, tf, tc)

    call tick(t0)
    do r=1,REPS_TRIG; do i=1,N; qres(i)=tan(qbnd(i)); end do; call qfeed(qbnd); end do
    tq = elapsed(t0)
    call tick(t0)
    do r=1,REPS_TRIG; do i=1,N; fres(i)=tan(fbnd(i)); end do; call ffeed(fbnd); end do
    tf = elapsed(t0)
    call tick(t0)
    do r=1,REPS_TRIG; do i=1,N; cres(i)=c_dd_tan(cbnd(i)); end do; call cfeed(cbnd); end do
    tc = elapsed(t0)
    call report3("tan", n_ops, tq, tf, tc)
  end subroutine

  subroutine bench_inv_trig()
    integer :: i, r
    integer(int64) :: t0, n_ops
    real(real64) :: tq, tf, tc
    n_ops = int(N, int64) * int(REPS_TRIG, int64)

    call tick(t0)
    do r=1,REPS_TRIG; do i=1,N; qres(i)=asin(qbnd(i)); end do; call qfeed(qbnd); end do
    tq = elapsed(t0)
    call tick(t0)
    do r=1,REPS_TRIG; do i=1,N; fres(i)=asin(fbnd(i)); end do; call ffeed(fbnd); end do
    tf = elapsed(t0)
    call tick(t0)
    do r=1,REPS_TRIG; do i=1,N; cres(i)=c_dd_asin(cbnd(i)); end do; call cfeed(cbnd); end do
    tc = elapsed(t0)
    call report3("asin", n_ops, tq, tf, tc)

    call tick(t0)
    do r=1,REPS_TRIG; do i=1,N; qres(i)=acos(qbnd(i)); end do; call qfeed(qbnd); end do
    tq = elapsed(t0)
    call tick(t0)
    do r=1,REPS_TRIG; do i=1,N; fres(i)=acos(fbnd(i)); end do; call ffeed(fbnd); end do
    tf = elapsed(t0)
    call tick(t0)
    do r=1,REPS_TRIG; do i=1,N; cres(i)=c_dd_acos(cbnd(i)); end do; call cfeed(cbnd); end do
    tc = elapsed(t0)
    call report3("acos", n_ops, tq, tf, tc)

    call tick(t0)
    do r=1,REPS_TRIG; do i=1,N; qres(i)=atan(q1(i)); end do; call qfeed(q1); end do
    tq = elapsed(t0)
    call tick(t0)
    do r=1,REPS_TRIG; do i=1,N; fres(i)=atan(f1(i)); end do; call ffeed(f1); end do
    tf = elapsed(t0)
    call tick(t0)
    do r=1,REPS_TRIG; do i=1,N; cres(i)=c_dd_atan(c1(i)); end do; call cfeed(c1); end do
    tc = elapsed(t0)
    call report3("atan", n_ops, tq, tf, tc)

    call tick(t0)
    do r=1,REPS_TRIG; do i=1,N; qres(i)=atan2(q1(i),q2(i)); end do; call qfeed(q1); end do
    tq = elapsed(t0)
    call tick(t0)
    do r=1,REPS_TRIG; do i=1,N; fres(i)=atan2(f1(i),f2(i)); end do; call ffeed(f1); end do
    tf = elapsed(t0)
    call tick(t0)
    do r=1,REPS_TRIG; do i=1,N; cres(i)=c_dd_atan2(c1(i),c2(i)); end do; call cfeed(c1); end do
    tc = elapsed(t0)
    call report3("atan2", n_ops, tq, tf, tc)
  end subroutine

  subroutine bench_hyperbolic()
    integer :: i, r
    integer(int64) :: t0, n_ops
    real(real64) :: tq, tf, tc
    n_ops = int(N, int64) * int(REPS_TRIG, int64)

    call tick(t0)
    do r=1,REPS_TRIG; do i=1,N; qres(i)=sinh(qsmall(i)); end do; call qfeed(qsmall); end do
    tq = elapsed(t0)
    call tick(t0)
    do r=1,REPS_TRIG; do i=1,N; fres(i)=sinh(fsmall(i)); end do; call ffeed(fsmall); end do
    tf = elapsed(t0)
    call tick(t0)
    do r=1,REPS_TRIG; do i=1,N; cres(i)=c_dd_sinh(csmall(i)); end do; call cfeed(csmall); end do
    tc = elapsed(t0)
    call report3("sinh", n_ops, tq, tf, tc)

    call tick(t0)
    do r=1,REPS_TRIG; do i=1,N; qres(i)=cosh(qsmall(i)); end do; call qfeed(qsmall); end do
    tq = elapsed(t0)
    call tick(t0)
    do r=1,REPS_TRIG; do i=1,N; fres(i)=cosh(fsmall(i)); end do; call ffeed(fsmall); end do
    tf = elapsed(t0)
    call tick(t0)
    do r=1,REPS_TRIG; do i=1,N; cres(i)=c_dd_cosh(csmall(i)); end do; call cfeed(csmall); end do
    tc = elapsed(t0)
    call report3("cosh", n_ops, tq, tf, tc)

    call tick(t0)
    do r=1,REPS_TRIG; do i=1,N; qres(i)=tanh(qsmall(i)); end do; call qfeed(qsmall); end do
    tq = elapsed(t0)
    call tick(t0)
    do r=1,REPS_TRIG; do i=1,N; fres(i)=tanh(fsmall(i)); end do; call ffeed(fsmall); end do
    tf = elapsed(t0)
    call tick(t0)
    do r=1,REPS_TRIG; do i=1,N; cres(i)=c_dd_tanh(csmall(i)); end do; call cfeed(csmall); end do
    tc = elapsed(t0)
    call report3("tanh", n_ops, tq, tf, tc)
  end subroutine

  subroutine bench_inv_hyperbolic()
    integer :: i, r
    integer(int64) :: t0, n_ops
    real(real64) :: tq, tf, tc
    n_ops = int(N, int64) * int(REPS_TRIG, int64)

    call tick(t0)
    do r=1,REPS_TRIG; do i=1,N; qres(i)=asinh(q1(i)); end do; call qfeed(q1); end do
    tq = elapsed(t0)
    call tick(t0)
    do r=1,REPS_TRIG; do i=1,N; fres(i)=asinh(f1(i)); end do; call ffeed(f1); end do
    tf = elapsed(t0)
    call tick(t0)
    do r=1,REPS_TRIG; do i=1,N; cres(i)=c_dd_asinh(c1(i)); end do; call cfeed(c1); end do
    tc = elapsed(t0)
    call report3("asinh", n_ops, tq, tf, tc)

    call tick(t0)
    do r=1,REPS_TRIG; do i=1,N; qres(i)=acosh(1.0_qp+qpos(i)); end do; call qfeed(qpos); end do
    tq = elapsed(t0)
    call tick(t0)
    do r=1,REPS_TRIG; do i=1,N; fres(i)=acosh(1.0_dp+fpos(i)); end do; call ffeed(fpos); end do
    tf = elapsed(t0)
    call tick(t0)
    do r=1,REPS_TRIG; do i=1,N
      block; type(float64x2_t) :: one; one%hi=1.0_dp; one%lo=0.0_dp
        cres(i)=c_dd_acosh(c_dd_add(one,cpos(i))); end block
    end do; call cfeed(cpos); end do
    tc = elapsed(t0)
    call report3("acosh", n_ops, tq, tf, tc)

    call tick(t0)
    do r=1,REPS_TRIG; do i=1,N; qres(i)=atanh(qbnd(i)); end do; call qfeed(qbnd); end do
    tq = elapsed(t0)
    call tick(t0)
    do r=1,REPS_TRIG; do i=1,N; fres(i)=atanh(fbnd(i)); end do; call ffeed(fbnd); end do
    tf = elapsed(t0)
    call tick(t0)
    do r=1,REPS_TRIG; do i=1,N; cres(i)=c_dd_atanh(cbnd(i)); end do; call cfeed(cbnd); end do
    tc = elapsed(t0)
    call report3("atanh", n_ops, tq, tf, tc)
  end subroutine

  subroutine bench_erf()
    integer :: i, r
    integer(int64) :: t0, n_ops
    real(real64) :: tq, tf, tc
    n_ops = int(N, int64) * int(REPS_VERY_SLOW, int64)

    call tick(t0)
    do r=1,REPS_VERY_SLOW; do i=1,N; qres(i)=erf(q1(i)); end do; call qfeed(q1); end do
    tq = elapsed(t0)
    call tick(t0)
    do r=1,REPS_VERY_SLOW; do i=1,N; fres(i)=erf(f1(i)); end do; call ffeed(f1); end do
    tf = elapsed(t0)
    call tick(t0)
    do r=1,REPS_VERY_SLOW; do i=1,N; cres(i)=c_dd_erf(c1(i)); end do; call cfeed(c1); end do
    tc = elapsed(t0)
    call report3("erf", n_ops, tq, tf, tc)

    call tick(t0)
    do r=1,REPS_VERY_SLOW; do i=1,N; qres(i)=erfc(q1(i)); end do; call qfeed(q1); end do
    tq = elapsed(t0)
    call tick(t0)
    do r=1,REPS_VERY_SLOW; do i=1,N; fres(i)=erfc(f1(i)); end do; call ffeed(f1); end do
    tf = elapsed(t0)
    call tick(t0)
    do r=1,REPS_VERY_SLOW; do i=1,N; cres(i)=c_dd_erfc(c1(i)); end do; call cfeed(c1); end do
    tc = elapsed(t0)
    call report3("erfc", n_ops, tq, tf, tc)
  end subroutine

end program
