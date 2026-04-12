program multifloat_bench
  ! Benchmark every kernel category exercised by test/fuzz.f90, comparing
  ! the multifloats float64x2 / complex128x2 implementations against the
  ! REAL(KIND=16) / COMPLEX(KIND=16) reference. The multifloats version is
  ! expected to be substantially faster — this program quantifies that.
  !
  ! Anti–dead-code-elimination strategy:
  !   - Each rep computes the entire output array.
  !   - At the end of each rep we read the FULL output via sum() into a
  !     scalar accumulator (forces all elements live), then feed the
  !     accumulator back into one element of an input array (forces a
  !     true cross-rep dependency). Without this, gfortran -O3 hoists or
  !     elides the entire rep loop because all reps would otherwise
  !     produce the same overwritten output.
  !   - The drain (sum + feedback) is INSIDE the timed region. It costs
  !     ≈N adds per rep, the same overhead in both the qp and mf legs,
  !     so the qp/mf speedup ratio is still meaningful even when the
  !     "real" op being benchmarked is itself a single add.
  use multifloats
  use, intrinsic :: iso_fortran_env, only: int64, real64
  implicit none

  integer, parameter :: qp = 16
  integer, parameter :: dp = 8

  ! Batch size and rep counts. Per-op total work = N * REPS, chosen so
  ! the fastest ops still take ~tens of ms in either leg.
  integer, parameter :: N = 1024
  integer, parameter :: REPS_FAST       = 400   ! +, -, *, /, abs, neg, ...
  integer, parameter :: REPS_TRIG       = 40    ! exp, log, sin, cos, ...
  integer, parameter :: REPS_VERY_SLOW  = 4     ! gamma, bessel, complex
  integer, parameter :: REPS_ARR        = 200000  ! sum, dot, matmul, ...

  ! Workspaces. Inputs are deliberately mutable: every rep nudges one
  ! element via the feedback path described above.
  real(qp), allocatable :: q1(:), q2(:), q3(:), qpos(:), qsmall(:), qbnd(:), qres(:)
  type(float64x2), allocatable :: f1(:), f2(:), f3(:), fpos(:), fsmall(:), fbnd(:), fres(:)
  complex(qp), allocatable :: cq1(:), cq2(:), cqres(:)
  type(complex128x2), allocatable :: cf1(:), cf2(:), cfres(:)

  ! Sinks: each op funnels its accumulator here, then the program prints
  ! both totals at the end so the optimizer cannot dead-strip the whole
  ! benchmark.
  real(qp) :: q_sink = 0.0_qp
  type(float64x2) :: f_sink
  complex(qp) :: cq_sink = (0.0_qp, 0.0_qp)
  type(complex128x2) :: cf_sink

  call init_data()
  f_sink     = float64x2(0.0_dp)
  cf_sink%re = float64x2(0.0_dp)
  cf_sink%im = float64x2(0.0_dp)

  print '(a)',         "================================================================"
  print '(a,i0,a,i0)', " multifloats benchmark — N=", N, " elements per rep, batched"
  print '(a)',         "================================================================"
  print '(a)',         ""
  print '(a)', " op                       n_ops      qp [s]      mf [s]    speedup"
  print '(a)', " ----------------------------------------------------------------"

  call bench_arith()
  call bench_unary_basic()
  call bench_minmax_binary()
  call bench_exp_log()
  call bench_trig()
  call bench_inv_trig()
  call bench_hyperbolic()
  call bench_inv_hyperbolic()
  call bench_erf_gamma()
  call bench_bessel()
  call bench_misc()
  call bench_complex()
  call bench_arrays()

  print '(a)', " ----------------------------------------------------------------"
  print '(a,es12.4,a,es12.4)', " sinks (must remain live): qp_sink=", &
      real(q_sink, dp), "  mf_sink=", f_sink%limbs(1) + f_sink%limbs(2)
  print '(a,es12.4,a,es12.4)', "                           cqp_re= ", &
      real(real(cq_sink, dp), dp), "  cmf_re=", cf_sink%re%limbs(1) + cf_sink%re%limbs(2)

contains

  ! ----------------------------------------------------------------
  ! Plumbing
  ! ----------------------------------------------------------------

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

  subroutine init_data()
    integer :: i
    real(dp) :: r(3)
    allocate(q1(N), q2(N), q3(N), qpos(N), qsmall(N), qbnd(N), qres(N))
    allocate(f1(N), f2(N), f3(N), fpos(N), fsmall(N), fbnd(N), fres(N))
    allocate(cq1(N), cq2(N), cqres(N))
    allocate(cf1(N), cf2(N), cfres(N))

    call random_seed()
    do i = 1, N
      call random_number(r)
      ! General arithmetic inputs in (-10, 10), nonzero (avoid div-by-0).
      q1(i) = real((r(1) - 0.5_dp) * 8.0_dp + sign(0.25_dp, r(1) - 0.5_dp), qp)
      q2(i) = real((r(2) - 0.5_dp) * 8.0_dp + sign(0.25_dp, r(2) - 0.5_dp), qp)
      q3(i) = real((r(3) - 0.5_dp) * 8.0_dp + sign(0.25_dp, r(3) - 0.5_dp), qp)
      ! Positive in (0.1, 10) for log / sqrt / pow base / gamma.
      qpos(i)   = real(0.1_dp + 9.9_dp * r(1), qp)
      ! Small in (-3, 3) for trig / hyperbolic.
      qsmall(i) = real((r(2) - 0.5_dp) * 6.0_dp, qp)
      ! Bounded in (-0.9, 0.9) for asin / acos / atanh / pow exponent.
      qbnd(i)   = real((r(3) - 0.5_dp) * 1.8_dp, qp)
    end do
    do i = 1, N
      f1(i)    = to_f64x2(q1(i))
      f2(i)    = to_f64x2(q2(i))
      f3(i)    = to_f64x2(q3(i))
      fpos(i)  = to_f64x2(qpos(i))
      fsmall(i)= to_f64x2(qsmall(i))
      fbnd(i)  = to_f64x2(qbnd(i))
      cq1(i)   = cmplx(qbnd(i), qsmall(i), qp)
      cq2(i)   = cmplx(qsmall(i), qbnd(i), qp)
      cf1(i)   = complex128x2(fbnd(i), fsmall(i))
      cf2(i)   = complex128x2(fsmall(i), fbnd(i))
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

  subroutine report(name, n_ops, tq, tf)
    character(*), intent(in) :: name
    integer(int64), intent(in) :: n_ops
    real(real64), intent(in) :: tq, tf
    real(real64) :: speed
    if (tf > 0.0_real64) then
      speed = tq / tf
    else
      speed = 0.0_real64
    end if
    write(*, '(1x, a22, 1x, i12, 1x, f10.4, 1x, f10.4, 1x, f9.2, "x")') &
        name, n_ops, tq, tf, speed
  end subroutine

  ! qp feedback: read every qres element into a dp accumulator (cheap,
  ! no libquadmath calls), then feed back into q_in(1) so the next rep
  ! sees a (slightly) different input. NOINLINE so the optimizer cannot
  ! fuse this with the rep-loop body and elide reps after the first.
  ! Draining via dp instead of qp/mf sum keeps drain cost negligible
  ! relative to the op being measured.
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

  ! mf feedback: drain via the leading dp limbs only.
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
    f_sink%limbs(1) = f_sink%limbs(1) + s
  end subroutine

  ! Complex feedback variants. Drain via dp on the real part only.
  subroutine cqfeed(cq_in)
    !GCC$ ATTRIBUTES NOINLINE :: cqfeed
    complex(qp), intent(inout) :: cq_in(:)
    real(dp) :: sr, si
    integer :: i
    sr = 0.0_dp; si = 0.0_dp
    do i = 1, size(cqres)
      sr = sr + real(real(cqres(i),  dp), dp)
      si = si + real(aimag(cqres(i)), dp)
    end do
    cq_in(1) = cq_in(1) + cmplx(real(sr, qp) * 1.0e-30_qp, 0.0_qp, qp)
    cq_sink = cq_sink + cmplx(real(sr, qp), real(si, qp), qp)
  end subroutine

  subroutine cffeed(cf_in)
    !GCC$ ATTRIBUTES NOINLINE :: cffeed
    type(complex128x2), intent(inout) :: cf_in(:)
    real(dp) :: sr, si
    integer :: i
    sr = 0.0_dp; si = 0.0_dp
    do i = 1, size(cfres)
      sr = sr + cfres(i)%re%limbs(1)
      si = si + cfres(i)%im%limbs(1)
    end do
    cf_in(1)%re%limbs(1) = cf_in(1)%re%limbs(1) + sr * 1.0e-30_dp
    cf_sink%re%limbs(1) = cf_sink%re%limbs(1) + sr
    cf_sink%im%limbs(1) = cf_sink%im%limbs(1) + si
  end subroutine

  ! ----------------------------------------------------------------
  ! Arithmetic kernels
  ! ----------------------------------------------------------------

  subroutine bench_arith()
    integer :: i, r
    integer(int64) :: t0, n_ops
    real(real64) :: tq, tf

    n_ops = int(N, int64) * int(REPS_FAST, int64)

    ! add
    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; qres(i) = q1(i) + q2(i); end do
      call qfeed(q1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; fres(i) = f1(i) + f2(i); end do
      call ffeed(f1)
    end do
    tf = elapsed(t0)
    call report("add", n_ops, tq, tf)

    ! sub
    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; qres(i) = q1(i) - q2(i); end do
      call qfeed(q1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; fres(i) = f1(i) - f2(i); end do
      call ffeed(f1)
    end do
    tf = elapsed(t0)
    call report("sub", n_ops, tq, tf)

    ! mul
    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; qres(i) = q1(i) * q2(i); end do
      call qfeed(q1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; fres(i) = f1(i) * f2(i); end do
      call ffeed(f1)
    end do
    tf = elapsed(t0)
    call report("mul", n_ops, tq, tf)

    ! div
    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; qres(i) = q1(i) / q2(i); end do
      call qfeed(q1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; fres(i) = f1(i) / f2(i); end do
      call ffeed(f1)
    end do
    tf = elapsed(t0)
    call report("div", n_ops, tq, tf)

    ! sqrt
    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; qres(i) = sqrt(qpos(i)); end do
      call qfeed(qpos)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; fres(i) = sqrt(fpos(i)); end do
      call ffeed(fpos)
    end do
    tf = elapsed(t0)
    call report("sqrt", n_ops, tq, tf)

    ! mixed-mode add (mf + dp)
    block
      real(dp) :: d
      call tick(t0)
      do r = 1, REPS_FAST
        d = real(q2(1), dp)
        do i = 1, N; qres(i) = q1(i) + real(d, qp); end do
        call qfeed(q1)
      end do
      tq = elapsed(t0)
      call tick(t0)
      do r = 1, REPS_FAST
        d = f2(1)%limbs(1)
        do i = 1, N; fres(i) = f1(i) + d; end do
        call ffeed(f1)
      end do
      tf = elapsed(t0)
      call report("add (mf+dp)", n_ops, tq, tf)
    end block

    ! mixed-mode mul (dp * mf)
    block
      real(dp) :: d
      call tick(t0)
      do r = 1, REPS_FAST
        d = real(q1(1), dp)
        do i = 1, N; qres(i) = real(d, qp) * q2(i); end do
        call qfeed(q2)
      end do
      tq = elapsed(t0)
      call tick(t0)
      do r = 1, REPS_FAST
        d = f1(1)%limbs(1)
        do i = 1, N; fres(i) = d * f2(i); end do
        call ffeed(f2)
      end do
      tf = elapsed(t0)
      call report("mul (dp*mf)", n_ops, tq, tf)
    end block
  end subroutine

  ! ----------------------------------------------------------------
  ! Unary basics: abs, neg, aint, anint, fraction
  ! ----------------------------------------------------------------

  subroutine bench_unary_basic()
    integer :: i, r
    integer(int64) :: t0, n_ops
    real(real64) :: tq, tf
    n_ops = int(N, int64) * int(REPS_FAST, int64)

    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; qres(i) = abs(q1(i)); end do
      call qfeed(q1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; fres(i) = abs(f1(i)); end do
      call ffeed(f1)
    end do
    tf = elapsed(t0)
    call report("abs", n_ops, tq, tf)

    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; qres(i) = -q1(i); end do
      call qfeed(q1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; fres(i) = -f1(i); end do
      call ffeed(f1)
    end do
    tf = elapsed(t0)
    call report("neg", n_ops, tq, tf)

    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; qres(i) = aint(q1(i)); end do
      call qfeed(q1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; fres(i) = aint(f1(i)); end do
      call ffeed(f1)
    end do
    tf = elapsed(t0)
    call report("aint", n_ops, tq, tf)

    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; qres(i) = anint(q1(i)); end do
      call qfeed(q1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; fres(i) = anint(f1(i)); end do
      call ffeed(f1)
    end do
    tf = elapsed(t0)
    call report("anint", n_ops, tq, tf)

    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; qres(i) = fraction(q1(i)); end do
      call qfeed(q1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; fres(i) = fraction(f1(i)); end do
      call ffeed(f1)
    end do
    tf = elapsed(t0)
    call report("fraction", n_ops, tq, tf)
  end subroutine

  ! ----------------------------------------------------------------
  ! min, max, sign, dim, mod, modulo, hypot, min3, max3
  ! ----------------------------------------------------------------

  subroutine bench_minmax_binary()
    integer :: i, r
    integer(int64) :: t0, n_ops
    real(real64) :: tq, tf
    n_ops = int(N, int64) * int(REPS_FAST, int64)

    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; qres(i) = min(q1(i), q2(i)); end do
      call qfeed(q1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; fres(i) = min(f1(i), f2(i)); end do
      call ffeed(f1)
    end do
    tf = elapsed(t0)
    call report("min", n_ops, tq, tf)

    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; qres(i) = max(q1(i), q2(i)); end do
      call qfeed(q1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; fres(i) = max(f1(i), f2(i)); end do
      call ffeed(f1)
    end do
    tf = elapsed(t0)
    call report("max", n_ops, tq, tf)

    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; qres(i) = sign(q1(i), q2(i)); end do
      call qfeed(q1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; fres(i) = sign(f1(i), f2(i)); end do
      call ffeed(f1)
    end do
    tf = elapsed(t0)
    call report("sign", n_ops, tq, tf)

    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; qres(i) = dim(q1(i), q2(i)); end do
      call qfeed(q1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; fres(i) = dim(f1(i), f2(i)); end do
      call ffeed(f1)
    end do
    tf = elapsed(t0)
    call report("dim", n_ops, tq, tf)

    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; qres(i) = mod(q1(i), q2(i)); end do
      call qfeed(q1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; fres(i) = mod(f1(i), f2(i)); end do
      call ffeed(f1)
    end do
    tf = elapsed(t0)
    call report("mod", n_ops, tq, tf)

    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; qres(i) = modulo(q1(i), q2(i)); end do
      call qfeed(q1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; fres(i) = modulo(f1(i), f2(i)); end do
      call ffeed(f1)
    end do
    tf = elapsed(t0)
    call report("modulo", n_ops, tq, tf)

    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; qres(i) = hypot(q1(i), q2(i)); end do
      call qfeed(q1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; fres(i) = hypot(f1(i), f2(i)); end do
      call ffeed(f1)
    end do
    tf = elapsed(t0)
    call report("hypot", n_ops, tq, tf)

    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; qres(i) = min(q1(i), q2(i), q3(i)); end do
      call qfeed(q1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; fres(i) = min(f1(i), f2(i), f3(i)); end do
      call ffeed(f1)
    end do
    tf = elapsed(t0)
    call report("min3", n_ops, tq, tf)

    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; qres(i) = max(q1(i), q2(i), q3(i)); end do
      call qfeed(q1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; fres(i) = max(f1(i), f2(i), f3(i)); end do
      call ffeed(f1)
    end do
    tf = elapsed(t0)
    call report("max3", n_ops, tq, tf)
  end subroutine

  ! ----------------------------------------------------------------
  ! exp / log / log10 / pow
  ! ----------------------------------------------------------------

  subroutine bench_exp_log()
    integer :: i, r
    integer(int64) :: t0, n_ops
    real(real64) :: tq, tf
    n_ops = int(N, int64) * int(REPS_TRIG, int64)

    call tick(t0)
    do r = 1, REPS_TRIG
      do i = 1, N; qres(i) = exp(qsmall(i)); end do
      call qfeed(qsmall)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_TRIG
      do i = 1, N; fres(i) = exp(fsmall(i)); end do
      call ffeed(fsmall)
    end do
    tf = elapsed(t0)
    call report("exp", n_ops, tq, tf)

    call tick(t0)
    do r = 1, REPS_TRIG
      do i = 1, N; qres(i) = log(qpos(i)); end do
      call qfeed(qpos)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_TRIG
      do i = 1, N; fres(i) = log(fpos(i)); end do
      call ffeed(fpos)
    end do
    tf = elapsed(t0)
    call report("log", n_ops, tq, tf)

    call tick(t0)
    do r = 1, REPS_TRIG
      do i = 1, N; qres(i) = log10(qpos(i)); end do
      call qfeed(qpos)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_TRIG
      do i = 1, N; fres(i) = log10(fpos(i)); end do
      call ffeed(fpos)
    end do
    tf = elapsed(t0)
    call report("log10", n_ops, tq, tf)

    call tick(t0)
    do r = 1, REPS_TRIG
      do i = 1, N; qres(i) = qpos(i) ** qbnd(i); end do
      call qfeed(qpos)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_TRIG
      do i = 1, N; fres(i) = fpos(i) ** fbnd(i); end do
      call ffeed(fpos)
    end do
    tf = elapsed(t0)
    call report("pow", n_ops, tq, tf)

    call tick(t0)
    do r = 1, REPS_TRIG
      do i = 1, N; qres(i) = q1(i) ** 3; end do
      call qfeed(q1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_TRIG
      do i = 1, N; fres(i) = f1(i) ** 3; end do
      call ffeed(f1)
    end do
    tf = elapsed(t0)
    call report("pow_int", n_ops, tq, tf)
  end subroutine

  ! ----------------------------------------------------------------
  ! sin, cos, tan
  ! ----------------------------------------------------------------

  subroutine bench_trig()
    integer :: i, r
    integer(int64) :: t0, n_ops
    real(real64) :: tq, tf
    n_ops = int(N, int64) * int(REPS_TRIG, int64)

    call tick(t0)
    do r = 1, REPS_TRIG
      do i = 1, N; qres(i) = sin(qsmall(i)); end do
      call qfeed(qsmall)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_TRIG
      do i = 1, N; fres(i) = sin(fsmall(i)); end do
      call ffeed(fsmall)
    end do
    tf = elapsed(t0)
    call report("sin", n_ops, tq, tf)

    call tick(t0)
    do r = 1, REPS_TRIG
      do i = 1, N; qres(i) = cos(qsmall(i)); end do
      call qfeed(qsmall)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_TRIG
      do i = 1, N; fres(i) = cos(fsmall(i)); end do
      call ffeed(fsmall)
    end do
    tf = elapsed(t0)
    call report("cos", n_ops, tq, tf)

    call tick(t0)
    do r = 1, REPS_TRIG
      do i = 1, N; qres(i) = sin(qsmall(i) * acos(-1.0_qp)); end do
      call qfeed(qsmall)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_TRIG
      do i = 1, N; fres(i) = sinpi(fsmall(i)); end do
      call ffeed(fsmall)
    end do
    tf = elapsed(t0)
    call report("sinpi", n_ops, tq, tf)

    call tick(t0)
    do r = 1, REPS_TRIG
      do i = 1, N; qres(i) = cos(qsmall(i) * acos(-1.0_qp)); end do
      call qfeed(qsmall)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_TRIG
      do i = 1, N; fres(i) = cospi(fsmall(i)); end do
      call ffeed(fsmall)
    end do
    tf = elapsed(t0)
    call report("cospi", n_ops, tq, tf)

    call tick(t0)
    do r = 1, REPS_TRIG
      do i = 1, N; qres(i) = tan(qbnd(i)); end do
      call qfeed(qbnd)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_TRIG
      do i = 1, N; fres(i) = tan(fbnd(i)); end do
      call ffeed(fbnd)
    end do
    tf = elapsed(t0)
    call report("tan", n_ops, tq, tf)
  end subroutine

  ! ----------------------------------------------------------------
  ! asin, acos, atan, atan2
  ! ----------------------------------------------------------------

  subroutine bench_inv_trig()
    integer :: i, r
    integer(int64) :: t0, n_ops
    real(real64) :: tq, tf
    n_ops = int(N, int64) * int(REPS_TRIG, int64)

    call tick(t0)
    do r = 1, REPS_TRIG
      do i = 1, N; qres(i) = asin(qbnd(i)); end do
      call qfeed(qbnd)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_TRIG
      do i = 1, N; fres(i) = asin(fbnd(i)); end do
      call ffeed(fbnd)
    end do
    tf = elapsed(t0)
    call report("asin", n_ops, tq, tf)

    call tick(t0)
    do r = 1, REPS_TRIG
      do i = 1, N; qres(i) = acos(qbnd(i)); end do
      call qfeed(qbnd)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_TRIG
      do i = 1, N; fres(i) = acos(fbnd(i)); end do
      call ffeed(fbnd)
    end do
    tf = elapsed(t0)
    call report("acos", n_ops, tq, tf)

    call tick(t0)
    do r = 1, REPS_TRIG
      do i = 1, N; qres(i) = atan(q1(i)); end do
      call qfeed(q1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_TRIG
      do i = 1, N; fres(i) = atan(f1(i)); end do
      call ffeed(f1)
    end do
    tf = elapsed(t0)
    call report("atan", n_ops, tq, tf)

    call tick(t0)
    do r = 1, REPS_TRIG
      do i = 1, N; qres(i) = atan2(q1(i), q2(i)); end do
      call qfeed(q1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_TRIG
      do i = 1, N; fres(i) = atan2(f1(i), f2(i)); end do
      call ffeed(f1)
    end do
    tf = elapsed(t0)
    call report("atan2", n_ops, tq, tf)
  end subroutine

  ! ----------------------------------------------------------------
  ! sinh, cosh, tanh
  ! ----------------------------------------------------------------

  subroutine bench_hyperbolic()
    integer :: i, r
    integer(int64) :: t0, n_ops
    real(real64) :: tq, tf
    n_ops = int(N, int64) * int(REPS_TRIG, int64)

    call tick(t0)
    do r = 1, REPS_TRIG
      do i = 1, N; qres(i) = sinh(qsmall(i)); end do
      call qfeed(qsmall)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_TRIG
      do i = 1, N; fres(i) = sinh(fsmall(i)); end do
      call ffeed(fsmall)
    end do
    tf = elapsed(t0)
    call report("sinh", n_ops, tq, tf)

    call tick(t0)
    do r = 1, REPS_TRIG
      do i = 1, N; qres(i) = cosh(qsmall(i)); end do
      call qfeed(qsmall)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_TRIG
      do i = 1, N; fres(i) = cosh(fsmall(i)); end do
      call ffeed(fsmall)
    end do
    tf = elapsed(t0)
    call report("cosh", n_ops, tq, tf)

    call tick(t0)
    do r = 1, REPS_TRIG
      do i = 1, N; qres(i) = tanh(qsmall(i)); end do
      call qfeed(qsmall)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_TRIG
      do i = 1, N; fres(i) = tanh(fsmall(i)); end do
      call ffeed(fsmall)
    end do
    tf = elapsed(t0)
    call report("tanh", n_ops, tq, tf)
  end subroutine

  ! ----------------------------------------------------------------
  ! asinh, acosh, atanh
  ! ----------------------------------------------------------------

  subroutine bench_inv_hyperbolic()
    integer :: i, r
    integer(int64) :: t0, n_ops
    real(real64) :: tq, tf
    n_ops = int(N, int64) * int(REPS_TRIG, int64)

    call tick(t0)
    do r = 1, REPS_TRIG
      do i = 1, N; qres(i) = asinh(q1(i)); end do
      call qfeed(q1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_TRIG
      do i = 1, N; fres(i) = asinh(f1(i)); end do
      call ffeed(f1)
    end do
    tf = elapsed(t0)
    call report("asinh", n_ops, tq, tf)

    ! acosh requires x >= 1: use 1 + qpos
    call tick(t0)
    do r = 1, REPS_TRIG
      do i = 1, N; qres(i) = acosh(1.0_qp + qpos(i)); end do
      call qfeed(qpos)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_TRIG
      do i = 1, N; fres(i) = acosh(1.0_dp + fpos(i)); end do
      call ffeed(fpos)
    end do
    tf = elapsed(t0)
    call report("acosh", n_ops, tq, tf)

    call tick(t0)
    do r = 1, REPS_TRIG
      do i = 1, N; qres(i) = atanh(qbnd(i)); end do
      call qfeed(qbnd)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_TRIG
      do i = 1, N; fres(i) = atanh(fbnd(i)); end do
      call ffeed(fbnd)
    end do
    tf = elapsed(t0)
    call report("atanh", n_ops, tq, tf)
  end subroutine

  ! ----------------------------------------------------------------
  ! erf, erfc, erfc_scaled, gamma, log_gamma
  ! ----------------------------------------------------------------

  subroutine bench_erf_gamma()
    integer :: i, r
    integer(int64) :: t0, n_ops
    real(real64) :: tq, tf
    n_ops = int(N, int64) * int(REPS_VERY_SLOW, int64)

    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; qres(i) = erf(q1(i)); end do
      call qfeed(q1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; fres(i) = erf(f1(i)); end do
      call ffeed(f1)
    end do
    tf = elapsed(t0)
    call report("erf", n_ops, tq, tf)

    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; qres(i) = erfc(q1(i)); end do
      call qfeed(q1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; fres(i) = erfc(f1(i)); end do
      call ffeed(f1)
    end do
    tf = elapsed(t0)
    call report("erfc", n_ops, tq, tf)

    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; qres(i) = erfc_scaled(q1(i)); end do
      call qfeed(q1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; fres(i) = erfc_scaled(f1(i)); end do
      call ffeed(f1)
    end do
    tf = elapsed(t0)
    call report("erfc_scaled", n_ops, tq, tf)

    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; qres(i) = gamma(qpos(i)); end do
      call qfeed(qpos)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; fres(i) = gamma(fpos(i)); end do
      call ffeed(fpos)
    end do
    tf = elapsed(t0)
    call report("gamma", n_ops, tq, tf)

    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; qres(i) = log_gamma(qpos(i)); end do
      call qfeed(qpos)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; fres(i) = log_gamma(fpos(i)); end do
      call ffeed(fpos)
    end do
    tf = elapsed(t0)
    call report("log_gamma", n_ops, tq, tf)
  end subroutine

  ! ----------------------------------------------------------------
  ! Bessel J0, J1, Jn, Y0, Y1, Yn
  ! ----------------------------------------------------------------

  subroutine bench_bessel()
    integer :: i, r
    integer(int64) :: t0, n_ops
    real(real64) :: tq, tf
    n_ops = int(N, int64) * int(REPS_VERY_SLOW, int64)

    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; qres(i) = bessel_j0(q1(i)); end do
      call qfeed(q1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; fres(i) = bessel_j0(f1(i)); end do
      call ffeed(f1)
    end do
    tf = elapsed(t0)
    call report("bessel_j0", n_ops, tq, tf)

    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; qres(i) = bessel_j1(q1(i)); end do
      call qfeed(q1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; fres(i) = bessel_j1(f1(i)); end do
      call ffeed(f1)
    end do
    tf = elapsed(t0)
    call report("bessel_j1", n_ops, tq, tf)

    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; qres(i) = bessel_jn(3, q1(i)); end do
      call qfeed(q1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; fres(i) = bessel_jn(3, f1(i)); end do
      call ffeed(f1)
    end do
    tf = elapsed(t0)
    call report("bessel_jn(3,.)", n_ops, tq, tf)

    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; qres(i) = bessel_y0(qpos(i)); end do
      call qfeed(qpos)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; fres(i) = bessel_y0(fpos(i)); end do
      call ffeed(fpos)
    end do
    tf = elapsed(t0)
    call report("bessel_y0", n_ops, tq, tf)

    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; qres(i) = bessel_y1(qpos(i)); end do
      call qfeed(qpos)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; fres(i) = bessel_y1(fpos(i)); end do
      call ffeed(fpos)
    end do
    tf = elapsed(t0)
    call report("bessel_y1", n_ops, tq, tf)

    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; qres(i) = bessel_yn(3, qpos(i)); end do
      call qfeed(qpos)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; fres(i) = bessel_yn(3, fpos(i)); end do
      call ffeed(fpos)
    end do
    tf = elapsed(t0)
    call report("bessel_yn(3,.)", n_ops, tq, tf)
  end subroutine

  ! ----------------------------------------------------------------
  ! scale, set_exponent
  ! ----------------------------------------------------------------

  subroutine bench_misc()
    integer :: i, r
    integer(int64) :: t0, n_ops
    real(real64) :: tq, tf
    n_ops = int(N, int64) * int(REPS_FAST, int64)

    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; qres(i) = scale(q1(i), 5); end do
      call qfeed(q1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; fres(i) = scale(f1(i), 5); end do
      call ffeed(f1)
    end do
    tf = elapsed(t0)
    call report("scale", n_ops, tq, tf)

    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; qres(i) = set_exponent(q1(i), 5); end do
      call qfeed(q1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; fres(i) = set_exponent(f1(i), 5); end do
      call ffeed(f1)
    end do
    tf = elapsed(t0)
    call report("set_exponent", n_ops, tq, tf)
  end subroutine

  ! ----------------------------------------------------------------
  ! Complex arithmetic and transcendentals
  ! ----------------------------------------------------------------

  subroutine bench_complex()
    integer :: i, r
    integer(int64) :: t0, n_ops_fast, n_ops_slow
    real(real64) :: tq, tf
    n_ops_fast = int(N, int64) * int(REPS_FAST, int64)
    n_ops_slow = int(N, int64) * int(REPS_VERY_SLOW, int64)

    ! cx_add
    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; cqres(i) = cq1(i) + cq2(i); end do
      call cqfeed(cq1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; cfres(i) = cf1(i) + cf2(i); end do
      call cffeed(cf1)
    end do
    tf = elapsed(t0)
    call report("cx_add", n_ops_fast, tq, tf)

    ! cx_sub
    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; cqres(i) = cq1(i) - cq2(i); end do
      call cqfeed(cq1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; cfres(i) = cf1(i) - cf2(i); end do
      call cffeed(cf1)
    end do
    tf = elapsed(t0)
    call report("cx_sub", n_ops_fast, tq, tf)

    ! cx_mul
    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; cqres(i) = cq1(i) * cq2(i); end do
      call cqfeed(cq1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; cfres(i) = cf1(i) * cf2(i); end do
      call cffeed(cf1)
    end do
    tf = elapsed(t0)
    call report("cx_mul", n_ops_fast, tq, tf)

    ! cx_div
    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; cqres(i) = cq1(i) / cq2(i); end do
      call cqfeed(cq1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; cfres(i) = cf1(i) / cf2(i); end do
      call cffeed(cf1)
    end do
    tf = elapsed(t0)
    call report("cx_div", n_ops_fast, tq, tf)

    ! cx_conjg
    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; cqres(i) = conjg(cq1(i)); end do
      call cqfeed(cq1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_FAST
      do i = 1, N; cfres(i) = conjg(cf1(i)); end do
      call cffeed(cf1)
    end do
    tf = elapsed(t0)
    call report("cx_conjg", n_ops_fast, tq, tf)

    ! cx_abs (complex → real)
    block
      real(qp) :: qs
      type(float64x2) :: fs
      call tick(t0)
      do r = 1, REPS_FAST
        do i = 1, N; qres(i) = abs(cq1(i)); end do
        qs = sum(qres)
        q_sink = q_sink + qs
        cq1(1) = cq1(1) + cmplx(qs * 1.0e-30_qp, 0.0_qp, qp)
      end do
      tq = elapsed(t0)
      call tick(t0)
      do r = 1, REPS_FAST
        do i = 1, N; fres(i) = abs(cf1(i)); end do
        fs = sum(fres)
        f_sink = f_sink + fs
        cf1(1)%re%limbs(1) = cf1(1)%re%limbs(1) + fs%limbs(1) * 1.0e-30_dp
      end do
      tf = elapsed(t0)
      call report("cx_abs", n_ops_fast, tq, tf)
    end block

    ! cx_sqrt
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; cqres(i) = sqrt(cq1(i)); end do
      call cqfeed(cq1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; cfres(i) = sqrt(cf1(i)); end do
      call cffeed(cf1)
    end do
    tf = elapsed(t0)
    call report("cx_sqrt", n_ops_slow, tq, tf)

    ! cx_exp
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; cqres(i) = exp(cq1(i)); end do
      call cqfeed(cq1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; cfres(i) = exp(cf1(i)); end do
      call cffeed(cf1)
    end do
    tf = elapsed(t0)
    call report("cx_exp", n_ops_slow, tq, tf)

    ! cx_log
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; cqres(i) = log(cq1(i)); end do
      call cqfeed(cq1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; cfres(i) = log(cf1(i)); end do
      call cffeed(cf1)
    end do
    tf = elapsed(t0)
    call report("cx_log", n_ops_slow, tq, tf)

    ! cx_sin
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; cqres(i) = sin(cq1(i)); end do
      call cqfeed(cq1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; cfres(i) = sin(cf1(i)); end do
      call cffeed(cf1)
    end do
    tf = elapsed(t0)
    call report("cx_sin", n_ops_slow, tq, tf)

    ! cx_cos
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; cqres(i) = cos(cq1(i)); end do
      call cqfeed(cq1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; cfres(i) = cos(cf1(i)); end do
      call cffeed(cf1)
    end do
    tf = elapsed(t0)
    call report("cx_cos", n_ops_slow, tq, tf)

    ! cx_tan
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; cqres(i) = tan(cq1(i)); end do
      call cqfeed(cq1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; cfres(i) = tan(cf1(i)); end do
      call cffeed(cf1)
    end do
    tf = elapsed(t0)
    call report("cx_tan", n_ops_slow, tq, tf)

    ! cx_sinh
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; cqres(i) = sinh(cq1(i)); end do
      call cqfeed(cq1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; cfres(i) = sinh(cf1(i)); end do
      call cffeed(cf1)
    end do
    tf = elapsed(t0)
    call report("cx_sinh", n_ops_slow, tq, tf)

    ! cx_cosh
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; cqres(i) = cosh(cq1(i)); end do
      call cqfeed(cq1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; cfres(i) = cosh(cf1(i)); end do
      call cffeed(cf1)
    end do
    tf = elapsed(t0)
    call report("cx_cosh", n_ops_slow, tq, tf)

    ! cx_tanh
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; cqres(i) = tanh(cq1(i)); end do
      call cqfeed(cq1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; cfres(i) = tanh(cf1(i)); end do
      call cffeed(cf1)
    end do
    tf = elapsed(t0)
    call report("cx_tanh", n_ops_slow, tq, tf)

    ! cx_asin
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; cqres(i) = asin(cq1(i)); end do
      call cqfeed(cq1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; cfres(i) = asin(cf1(i)); end do
      call cffeed(cf1)
    end do
    tf = elapsed(t0)
    call report("cx_asin", n_ops_slow, tq, tf)

    ! cx_acos
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; cqres(i) = acos(cq1(i)); end do
      call cqfeed(cq1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; cfres(i) = acos(cf1(i)); end do
      call cffeed(cf1)
    end do
    tf = elapsed(t0)
    call report("cx_acos", n_ops_slow, tq, tf)

    ! cx_atan
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; cqres(i) = atan(cq1(i)); end do
      call cqfeed(cq1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; cfres(i) = atan(cf1(i)); end do
      call cffeed(cf1)
    end do
    tf = elapsed(t0)
    call report("cx_atan", n_ops_slow, tq, tf)

    ! cx_asinh
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; cqres(i) = asinh(cq1(i)); end do
      call cqfeed(cq1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; cfres(i) = asinh(cf1(i)); end do
      call cffeed(cf1)
    end do
    tf = elapsed(t0)
    call report("cx_asinh", n_ops_slow, tq, tf)

    ! cx_acosh
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; cqres(i) = acosh(cq1(i)); end do
      call cqfeed(cq1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; cfres(i) = acosh(cf1(i)); end do
      call cffeed(cf1)
    end do
    tf = elapsed(t0)
    call report("cx_acosh", n_ops_slow, tq, tf)

    ! cx_atanh
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; cqres(i) = atanh(cq1(i)); end do
      call cqfeed(cq1)
    end do
    tq = elapsed(t0)
    call tick(t0)
    do r = 1, REPS_VERY_SLOW
      do i = 1, N; cfres(i) = atanh(cf1(i)); end do
      call cffeed(cf1)
    end do
    tf = elapsed(t0)
    call report("cx_atanh", n_ops_slow, tq, tf)
  end subroutine

  ! ----------------------------------------------------------------
  ! Array reductions: sum, product, maxval, minval, dot_product, norm2,
  ! matmul. Each call processes a length-NN slice; the rep loop accumulates
  ! into a sink so each iteration's reduction stays live.
  ! ----------------------------------------------------------------

  subroutine bench_arrays()
    integer, parameter :: NN = 8
    integer :: r, k
    integer(int64) :: t0, n_ops
    real(real64) :: tq, tf
    real(qp) :: qa(NN), qb(NN), qm(NN, NN), qmv(NN), qacc
    type(float64x2) :: a(NN), b(NN), m(NN, NN), mv(NN), facc

    n_ops = int(NN, int64) * int(REPS_ARR, int64)

    do k = 1, NN
      qa(k) = q1(k); qb(k) = q2(k)
      a(k) = f1(k); b(k) = f2(k)
    end do
    do k = 1, NN
      qm(:, k) = q1(k * NN + 1 : k * NN + NN)
      m(:, k)  = f1(k * NN + 1 : k * NN + NN)
    end do

    ! sum
    call tick(t0)
    qacc = 0.0_qp
    do r = 1, REPS_ARR
      qacc = qacc + sum(qa)
      qa(1) = qa(1) + qacc * 1.0e-30_qp  ! prevent hoisting
    end do
    tq = elapsed(t0)
    q_sink = q_sink + qacc
    call tick(t0)
    facc = float64x2(0.0_dp)
    do r = 1, REPS_ARR
      facc = facc + sum(a)
      a(1)%limbs(1) = a(1)%limbs(1) + facc%limbs(1) * 1.0e-30_dp
    end do
    tf = elapsed(t0)
    f_sink = f_sink + facc
    call report("arr_sum (n=8)", n_ops, tq, tf)

    ! product
    call tick(t0)
    qacc = 0.0_qp
    do r = 1, REPS_ARR
      qacc = qacc + product(qa)
      qa(1) = qa(1) + qacc * 1.0e-30_qp
    end do
    tq = elapsed(t0)
    q_sink = q_sink + qacc
    call tick(t0)
    facc = float64x2(0.0_dp)
    do r = 1, REPS_ARR
      facc = facc + product(a)
      a(1)%limbs(1) = a(1)%limbs(1) + facc%limbs(1) * 1.0e-30_dp
    end do
    tf = elapsed(t0)
    f_sink = f_sink + facc
    call report("arr_product (n=8)", n_ops, tq, tf)

    ! maxval
    call tick(t0)
    qacc = 0.0_qp
    do r = 1, REPS_ARR
      qacc = qacc + maxval(qa)
      qa(1) = qa(1) + qacc * 1.0e-30_qp
    end do
    tq = elapsed(t0)
    q_sink = q_sink + qacc
    call tick(t0)
    facc = float64x2(0.0_dp)
    do r = 1, REPS_ARR
      facc = facc + maxval(a)
      a(1)%limbs(1) = a(1)%limbs(1) + facc%limbs(1) * 1.0e-30_dp
    end do
    tf = elapsed(t0)
    f_sink = f_sink + facc
    call report("arr_maxval (n=8)", n_ops, tq, tf)

    ! minval
    call tick(t0)
    qacc = 0.0_qp
    do r = 1, REPS_ARR
      qacc = qacc + minval(qa)
      qa(1) = qa(1) + qacc * 1.0e-30_qp
    end do
    tq = elapsed(t0)
    q_sink = q_sink + qacc
    call tick(t0)
    facc = float64x2(0.0_dp)
    do r = 1, REPS_ARR
      facc = facc + minval(a)
      a(1)%limbs(1) = a(1)%limbs(1) + facc%limbs(1) * 1.0e-30_dp
    end do
    tf = elapsed(t0)
    f_sink = f_sink + facc
    call report("arr_minval (n=8)", n_ops, tq, tf)

    ! dot_product
    call tick(t0)
    qacc = 0.0_qp
    do r = 1, REPS_ARR
      qacc = qacc + dot_product(qa, qb)
      qa(1) = qa(1) + qacc * 1.0e-30_qp
    end do
    tq = elapsed(t0)
    q_sink = q_sink + qacc
    call tick(t0)
    facc = float64x2(0.0_dp)
    do r = 1, REPS_ARR
      facc = facc + dot_product(a, b)
      a(1)%limbs(1) = a(1)%limbs(1) + facc%limbs(1) * 1.0e-30_dp
    end do
    tf = elapsed(t0)
    f_sink = f_sink + facc
    call report("arr_dot (n=8)", n_ops, tq, tf)

    ! norm2
    call tick(t0)
    qacc = 0.0_qp
    do r = 1, REPS_ARR
      qacc = qacc + norm2(qa)
      qa(1) = qa(1) + qacc * 1.0e-30_qp
    end do
    tq = elapsed(t0)
    q_sink = q_sink + qacc
    call tick(t0)
    facc = float64x2(0.0_dp)
    do r = 1, REPS_ARR
      facc = facc + norm2(a)
      a(1)%limbs(1) = a(1)%limbs(1) + facc%limbs(1) * 1.0e-30_dp
    end do
    tf = elapsed(t0)
    f_sink = f_sink + facc
    call report("arr_norm2 (n=8)", n_ops, tq, tf)

    ! matmul (8x8 * 8): roughly NN^2 = 64 scalar muls per call.
    n_ops = int(NN, int64) * int(NN, int64) * int(REPS_ARR, int64)
    call tick(t0)
    do r = 1, REPS_ARR
      qmv = matmul(qm, qa)
      qa(1) = qa(1) + qmv(1) * 1.0e-30_qp
    end do
    tq = elapsed(t0)
    q_sink = q_sink + sum(qmv)
    call tick(t0)
    do r = 1, REPS_ARR
      mv = matmul(m, a)
      a(1)%limbs(1) = a(1)%limbs(1) + mv(1)%limbs(1) * 1.0e-30_dp
    end do
    tf = elapsed(t0)
    f_sink = f_sink + sum(mv)
    call report("arr_matmul (8x8*8)", n_ops, tq, tf)
  end subroutine

end program
