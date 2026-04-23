program bench_abi
  ! Phase-1 overhead baseline: Fortran-native operators/intrinsics vs the
  ! C-ABI shim path (extern "C" from src/multifloats_math.cc).
  !
  ! For each scalar DD op we measure two legs:
  !   native : whatever the Fortran module currently does (inline-able)
  !   c_abi  : bind(c) call to the corresponding kernel in multifloats.h
  !
  ! Ratio dd/cw:
  !   ~ 1.0x  : parity (no overhead from the boundary)
  !   > 1.0x  : the C kernel wins despite the call
  !   < 1.0x  : the native inline wins; overhead cost of delegating
  !
  ! Iteration counts are sized so each leg runs ~0.1-2 s; variance per op
  ! is printed so the reader can judge precision of the ratio.
  use multifloats
  use, intrinsic :: iso_fortran_env, only: int64, real64
  use, intrinsic :: iso_c_binding,   only: c_double, c_int
  implicit none

  integer, parameter :: dp = 8
  integer, parameter :: N    = 4096      ! inner array length
  integer, parameter :: REPS = 16000     ! outer repeats (N*REPS ~= 6.5e7 ops)
  integer, parameter :: TRIALS = 5       ! trials per leg (report best, show spread)

  type, bind(c) :: dd_c
    real(c_double) :: hi, lo
  end type dd_c

  interface
    pure type(dd_c) function c_dd_add(a, b) bind(c, name='adddd')
      import :: dd_c; type(dd_c), value :: a, b
    end function
    pure type(dd_c) function c_dd_sub(a, b) bind(c, name='subdd')
      import :: dd_c; type(dd_c), value :: a, b
    end function
    pure type(dd_c) function c_dd_mul(a, b) bind(c, name='muldd')
      import :: dd_c; type(dd_c), value :: a, b
    end function
    pure type(dd_c) function c_dd_div(a, b) bind(c, name='divdd')
      import :: dd_c; type(dd_c), value :: a, b
    end function
    pure type(dd_c) function c_dd_neg(a) bind(c, name='negdd')
      import :: dd_c; type(dd_c), value :: a
    end function
    pure type(dd_c) function c_dd_abs(a) bind(c, name='fabsdd')
      import :: dd_c; type(dd_c), value :: a
    end function
    pure type(dd_c) function c_dd_trunc(a) bind(c, name='truncdd')
      import :: dd_c; type(dd_c), value :: a
    end function
    pure type(dd_c) function c_dd_round(a) bind(c, name='rounddd')
      import :: dd_c; type(dd_c), value :: a
    end function
    pure type(dd_c) function c_dd_sqrt(a) bind(c, name='sqrtdd')
      import :: dd_c; type(dd_c), value :: a
    end function
    pure type(dd_c) function c_dd_min(a, b) bind(c, name='fmindd')
      import :: dd_c; type(dd_c), value :: a, b
    end function
    pure type(dd_c) function c_dd_max(a, b) bind(c, name='fmaxdd')
      import :: dd_c; type(dd_c), value :: a, b
    end function
    pure type(dd_c) function c_dd_copysign(a, b) bind(c, name='copysigndd')
      import :: dd_c; type(dd_c), value :: a, b
    end function
    pure type(dd_c) function c_dd_fdim(a, b) bind(c, name='fdimdd')
      import :: dd_c; type(dd_c), value :: a, b
    end function
    pure integer(c_int) function c_dd_eq(a, b) bind(c, name='eqdd')
      import :: dd_c, c_int; type(dd_c), value :: a, b
    end function
    pure integer(c_int) function c_dd_ne(a, b) bind(c, name='nedd')
      import :: dd_c, c_int; type(dd_c), value :: a, b
    end function
    pure integer(c_int) function c_dd_lt(a, b) bind(c, name='ltdd')
      import :: dd_c, c_int; type(dd_c), value :: a, b
    end function
    pure integer(c_int) function c_dd_le(a, b) bind(c, name='ledd')
      import :: dd_c, c_int; type(dd_c), value :: a, b
    end function
    pure integer(c_int) function c_dd_gt(a, b) bind(c, name='gtdd')
      import :: dd_c, c_int; type(dd_c), value :: a, b
    end function
    pure integer(c_int) function c_dd_ge(a, b) bind(c, name='gedd')
      import :: dd_c, c_int; type(dd_c), value :: a, b
    end function
  end interface

  ! Data arrays (native Fortran type)
  type(float64x2), allocatable :: f1(:), f2(:), fres(:)
  ! Mirror arrays (C ABI struct)
  type(dd_c),      allocatable :: d1(:), d2(:), dres(:)
  ! Int/logical result buffers
  integer, allocatable :: ires_n(:), ires_c(:)

  real(dp) :: dd_sink = 0.0_dp, cw_sink = 0.0_dp
  integer  :: i

  allocate(f1(N), f2(N), fres(N), d1(N), d2(N), dres(N), &
           ires_n(N), ires_c(N))

  call seed_and_fill()

  print '(a)', "=============================================================="
  print '(a)', " Phase-1 ABI overhead: Fortran-native vs C-ABI (bind(c))"
  print '(a,i0,a,i0,a,i0,a,i0,a)', " N=", N, "  REPS=", REPS, &
        "  ops/leg=", N*REPS, "  trials=", TRIALS, " (best-of)"
  print '(a)', "=============================================================="
  print '(a)', ""
  print '(a)', " op              n_ops      dd [s]     cw [s]    dd/cw   ns/op(cw)"
  print '(a)', " --------------------------------------------------------------------"

  ! Unary: dd -> dd
  call bench("neg")
  call bench("abs")
  call bench("aint")
  call bench("anint")
  call bench("sqrt")

  ! Binary: dd,dd -> dd
  call bench("add")
  call bench("sub")
  call bench("mul")
  call bench("div")
  call bench("min")
  call bench("max")
  call bench("sign")
  call bench("dim")

  ! Binary predicates: dd,dd -> int
  call bench("eq")
  call bench("ne")
  call bench("lt")
  call bench("le")
  call bench("gt")
  call bench("ge")

  ! Phase-2 candidates (native only until the C ABI kernels land)
  call bench_native("pow7")
  call bench_native("pow20")
  call bench_native("modulo")

  print '(a)', " --------------------------------------------------------------------"
  print '(a,es11.3,a,es11.3)', " sinks: dd=", dd_sink, " cw=", cw_sink

contains

  subroutine seed_and_fill()
    integer :: seed_size
    integer, allocatable :: seed(:)
    real(dp) :: r, rlo
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed = 42
    call random_seed(put=seed)
    deallocate(seed)
    do i = 1, N
      call random_number(r); call random_number(rlo)
      f1(i)%limbs(1) = (r - 0.5_dp) * 8.0_dp + sign(0.25_dp, r - 0.5_dp)
      f1(i)%limbs(2) = f1(i)%limbs(1) * (rlo - 0.5_dp) * 2.0_dp**(-52)
      call random_number(r); call random_number(rlo)
      f2(i)%limbs(1) = (r - 0.5_dp) * 8.0_dp + sign(0.25_dp, r - 0.5_dp)
      f2(i)%limbs(2) = f2(i)%limbs(1) * (rlo - 0.5_dp) * 2.0_dp**(-52)
      d1(i)%hi = f1(i)%limbs(1); d1(i)%lo = f1(i)%limbs(2)
      d2(i)%hi = f2(i)%limbs(1); d2(i)%lo = f2(i)%limbs(2)
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

  subroutine dd_drain(arr)
    !GCC$ ATTRIBUTES NOINLINE :: dd_drain
    type(float64x2), intent(inout) :: arr(:)
    real(dp) :: s
    integer :: j
    s = 0.0_dp
    do j = 1, size(fres); s = s + fres(j)%limbs(1) + fres(j)%limbs(2); end do
    arr(1)%limbs(1) = arr(1)%limbs(1) + s * 1e-300_dp
    dd_sink = dd_sink + s
  end subroutine

  subroutine cw_drain(arr)
    !GCC$ ATTRIBUTES NOINLINE :: cw_drain
    type(dd_c), intent(inout) :: arr(:)
    real(dp) :: s
    integer :: j
    s = 0.0_dp
    do j = 1, size(dres); s = s + dres(j)%hi + dres(j)%lo; end do
    arr(1)%hi = arr(1)%hi + s * 1e-300_dp
    cw_sink = cw_sink + s
  end subroutine

  subroutine dd_drain_int(arr)
    !GCC$ ATTRIBUTES NOINLINE :: dd_drain_int
    integer, intent(inout) :: arr(:)
    integer :: s, j
    s = 0
    do j = 1, size(ires_n); s = s + ires_n(j); end do
    arr(1) = arr(1) + s
    dd_sink = dd_sink + real(s, dp)
  end subroutine

  subroutine cw_drain_int(arr)
    !GCC$ ATTRIBUTES NOINLINE :: cw_drain_int
    integer, intent(inout) :: arr(:)
    integer :: s, j
    s = 0
    do j = 1, size(ires_c); s = s + ires_c(j); end do
    arr(1) = arr(1) + s
    cw_sink = cw_sink + real(s, dp)
  end subroutine

  subroutine run_native(op, tsec)
    character(*), intent(in) :: op
    real(real64), intent(out) :: tsec
    integer(int64) :: t0
    integer :: r
    call tick(t0)
    select case (op)
    case ("neg");   do r=1,REPS; do i=1,N; fres(i) = -f1(i);                end do; call dd_drain(f1); end do
    case ("abs");   do r=1,REPS; do i=1,N; fres(i) = abs(f1(i));            end do; call dd_drain(f1); end do
    case ("aint");  do r=1,REPS; do i=1,N; fres(i) = aint(f1(i));           end do; call dd_drain(f1); end do
    case ("anint"); do r=1,REPS; do i=1,N; fres(i) = anint(f1(i));          end do; call dd_drain(f1); end do
    case ("sqrt");  do r=1,REPS; do i=1,N; fres(i) = sqrt(abs(f1(i)));      end do; call dd_drain(f1); end do
    case ("add");   do r=1,REPS; do i=1,N; fres(i) = f1(i) + f2(i);         end do; call dd_drain(f1); end do
    case ("sub");   do r=1,REPS; do i=1,N; fres(i) = f1(i) - f2(i);         end do; call dd_drain(f1); end do
    case ("mul");   do r=1,REPS; do i=1,N; fres(i) = f1(i) * f2(i);         end do; call dd_drain(f1); end do
    case ("div");   do r=1,REPS; do i=1,N; fres(i) = f1(i) / f2(i);         end do; call dd_drain(f1); end do
    case ("min");   do r=1,REPS; do i=1,N; fres(i) = min(f1(i), f2(i));     end do; call dd_drain(f1); end do
    case ("max");   do r=1,REPS; do i=1,N; fres(i) = max(f1(i), f2(i));     end do; call dd_drain(f1); end do
    case ("sign");  do r=1,REPS; do i=1,N; fres(i) = sign(f1(i), f2(i));    end do; call dd_drain(f1); end do
    case ("dim");   do r=1,REPS; do i=1,N; fres(i) = dim(f1(i), f2(i));     end do; call dd_drain(f1); end do
    case ("eq");    do r=1,REPS; do i=1,N; ires_n(i) = merge(1,0, f1(i) == f2(i)); end do; call dd_drain_int(ires_n); end do
    case ("ne");    do r=1,REPS; do i=1,N; ires_n(i) = merge(1,0, f1(i) /= f2(i)); end do; call dd_drain_int(ires_n); end do
    case ("lt");    do r=1,REPS; do i=1,N; ires_n(i) = merge(1,0, f1(i) <  f2(i)); end do; call dd_drain_int(ires_n); end do
    case ("le");    do r=1,REPS; do i=1,N; ires_n(i) = merge(1,0, f1(i) <= f2(i)); end do; call dd_drain_int(ires_n); end do
    case ("gt");    do r=1,REPS; do i=1,N; ires_n(i) = merge(1,0, f1(i) >  f2(i)); end do; call dd_drain_int(ires_n); end do
    case ("ge");    do r=1,REPS; do i=1,N; ires_n(i) = merge(1,0, f1(i) >= f2(i)); end do; call dd_drain_int(ires_n); end do
    case ("pow7");   do r=1,REPS; do i=1,N; fres(i) = f1(i) ** 7;   end do; call dd_drain(f1); end do
    case ("pow20");  do r=1,REPS; do i=1,N; fres(i) = f1(i) ** 20;  end do; call dd_drain(f1); end do
    case ("modulo"); do r=1,REPS; do i=1,N; fres(i) = modulo(f1(i), f2(i)); end do; call dd_drain(f1); end do
    end select
    tsec = elapsed(t0)
  end subroutine

  subroutine run_cabi(op, tsec)
    character(*), intent(in) :: op
    real(real64), intent(out) :: tsec
    integer(int64) :: t0
    integer :: r
    call tick(t0)
    select case (op)
    case ("neg");   do r=1,REPS; do i=1,N; dres(i) = c_dd_neg(d1(i));               end do; call cw_drain(d1); end do
    case ("abs");   do r=1,REPS; do i=1,N; dres(i) = c_dd_abs(d1(i));               end do; call cw_drain(d1); end do
    case ("aint");  do r=1,REPS; do i=1,N; dres(i) = c_dd_trunc(d1(i));             end do; call cw_drain(d1); end do
    case ("anint"); do r=1,REPS; do i=1,N; dres(i) = c_dd_round(d1(i));             end do; call cw_drain(d1); end do
    case ("sqrt");  do r=1,REPS; do i=1,N; dres(i) = c_dd_sqrt(c_dd_abs(d1(i)));    end do; call cw_drain(d1); end do
    case ("add");   do r=1,REPS; do i=1,N; dres(i) = c_dd_add(d1(i), d2(i));        end do; call cw_drain(d1); end do
    case ("sub");   do r=1,REPS; do i=1,N; dres(i) = c_dd_sub(d1(i), d2(i));        end do; call cw_drain(d1); end do
    case ("mul");   do r=1,REPS; do i=1,N; dres(i) = c_dd_mul(d1(i), d2(i));        end do; call cw_drain(d1); end do
    case ("div");   do r=1,REPS; do i=1,N; dres(i) = c_dd_div(d1(i), d2(i));        end do; call cw_drain(d1); end do
    case ("min");   do r=1,REPS; do i=1,N; dres(i) = c_dd_min(d1(i), d2(i));        end do; call cw_drain(d1); end do
    case ("max");   do r=1,REPS; do i=1,N; dres(i) = c_dd_max(d1(i), d2(i));        end do; call cw_drain(d1); end do
    case ("sign");  do r=1,REPS; do i=1,N; dres(i) = c_dd_copysign(d1(i), d2(i));   end do; call cw_drain(d1); end do
    case ("dim");   do r=1,REPS; do i=1,N; dres(i) = c_dd_fdim(d1(i), d2(i));       end do; call cw_drain(d1); end do
    case ("eq");    do r=1,REPS; do i=1,N; ires_c(i) = int(c_dd_eq(d1(i), d2(i)));  end do; call cw_drain_int(ires_c); end do
    case ("ne");    do r=1,REPS; do i=1,N; ires_c(i) = int(c_dd_ne(d1(i), d2(i)));  end do; call cw_drain_int(ires_c); end do
    case ("lt");    do r=1,REPS; do i=1,N; ires_c(i) = int(c_dd_lt(d1(i), d2(i)));  end do; call cw_drain_int(ires_c); end do
    case ("le");    do r=1,REPS; do i=1,N; ires_c(i) = int(c_dd_le(d1(i), d2(i)));  end do; call cw_drain_int(ires_c); end do
    case ("gt");    do r=1,REPS; do i=1,N; ires_c(i) = int(c_dd_gt(d1(i), d2(i)));  end do; call cw_drain_int(ires_c); end do
    case ("ge");    do r=1,REPS; do i=1,N; ires_c(i) = int(c_dd_ge(d1(i), d2(i)));  end do; call cw_drain_int(ires_c); end do
    end select
    tsec = elapsed(t0)
  end subroutine

  subroutine bench(op)
    character(*), intent(in) :: op
    integer(int64) :: n_ops
    real(real64) :: tn, tc, tn_best, tc_best, ratio, ns_per
    integer :: k

    n_ops = int(N, int64) * int(REPS, int64)

    ! warmup
    call run_native(op, tn)
    call run_cabi(op, tc)

    tn_best = huge(1.0_real64); tc_best = huge(1.0_real64)
    do k = 1, TRIALS
      call run_native(op, tn); if (tn < tn_best) tn_best = tn
      call run_cabi(op, tc);   if (tc < tc_best) tc_best = tc
    end do

    ratio = 0.0_real64; if (tc_best > 0.0_real64) ratio = tn_best / tc_best
    ns_per = 0.0_real64
    if (n_ops > 0) ns_per = tc_best * 1.0e9_real64 / real(n_ops, real64)

    write(*, '(1x,a10,1x,i12,2x,f10.4,1x,f10.4,2x,f6.2,"x",2x,f8.2)') &
        op, n_ops, tn_best, tc_best, ratio, ns_per
  end subroutine

  ! Native-only variant for ops that do not yet have a C-ABI shim. Prints
  ! the Fortran-native ns/op so a later stage can compare against it.
  subroutine bench_native(op)
    character(*), intent(in) :: op
    integer(int64) :: n_ops
    real(real64) :: tn, tn_best, ns_per
    integer :: k

    n_ops = int(N, int64) * int(REPS, int64)
    call run_native(op, tn)              ! warmup
    tn_best = huge(1.0_real64)
    do k = 1, TRIALS
      call run_native(op, tn); if (tn < tn_best) tn_best = tn
    end do
    ns_per = 0.0_real64
    if (n_ops > 0) ns_per = tn_best * 1.0e9_real64 / real(n_ops, real64)
    write(*, '(1x,a10,1x,i12,2x,f10.4,1x,"(native only)",6x,f8.2)') &
        op, n_ops, tn_best, ns_per
  end subroutine

end program
