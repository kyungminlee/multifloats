program bench_abi
  ! Test whether bind(c) + value on Fortran-native DD kernels matches
  ! the C wrapper speed. Three legs:
  !   fortran_mf : current elemental float64x2 operator (Fortran ABI)
  !   fortran_c  : same algorithm, but bind(c) with value args (C ABI)
  !   c_wrapper  : extern "C" wrapper calling C++ multifloats.hh (C ABI)
  !
  ! If fortran_c ≈ c_wrapper, the remaining gap is purely ABI, not codegen.
  use multifloats
  use dd_bindc
  use, intrinsic :: iso_fortran_env, only: int64, real64
  use, intrinsic :: iso_c_binding,   only: c_double
  implicit none

  integer, parameter :: dp = 8
  integer, parameter :: N = 1024
  integer, parameter :: REPS = 400

  ! C wrapper interfaces (from multifloats_c.cc)
  interface
    type(dd_c) function c_dd_add(a, b) bind(c, name='adddd')
      import :: dd_c; type(dd_c), value :: a, b
    end function
    type(dd_c) function c_dd_sub(a, b) bind(c, name='subdd')
      import :: dd_c; type(dd_c), value :: a, b
    end function
    type(dd_c) function c_dd_mul(a, b) bind(c, name='muldd')
      import :: dd_c; type(dd_c), value :: a, b
    end function
    type(dd_c) function c_dd_div(a, b) bind(c, name='divdd')
      import :: dd_c; type(dd_c), value :: a, b
    end function
    type(dd_c) function c_dd_sqrt(a) bind(c, name='sqrtdd')
      import :: dd_c; type(dd_c), value :: a
    end function
  end interface

  ! Data
  type(float64x2), allocatable :: f1(:), f2(:), fres(:)
  type(dd_c), allocatable :: d1(:), d2(:), dres(:)
  real(dp) :: mf_sink = 0.0_dp, fc_sink = 0.0_dp, cw_sink = 0.0_dp
  integer :: i

  allocate(f1(N), f2(N), fres(N), d1(N), d2(N), dres(N))
  block
    real(dp) :: r
    do i = 1, N
      call random_number(r)
      f1(i)%limbs(1) = (r - 0.5_dp) * 8.0_dp + sign(0.25_dp, r - 0.5_dp)
      f1(i)%limbs(2) = 0.0_dp
      call random_number(r)
      f2(i)%limbs(1) = (r - 0.5_dp) * 8.0_dp + sign(0.25_dp, r - 0.5_dp)
      f2(i)%limbs(2) = 0.0_dp
      d1(i)%hi = f1(i)%limbs(1); d1(i)%lo = f1(i)%limbs(2)
      d2(i)%hi = f2(i)%limbs(1); d2(i)%lo = f2(i)%limbs(2)
    end do
  end block

  print '(a)', "================================================================"
  print '(a)', " ABI test: Fortran-mf vs Fortran-bind(c) vs C-wrapper"
  print '(a)', "================================================================"
  print '(a)', ""
  print '(a)', &
    " op              n_ops    mf [s]    fc [s]    cw [s]  mf/fc   mf/cw   fc/cw"
  print '(a)', &
    " ---------------------------------------------------------------------------"

  call bench("add")
  call bench("sub")
  call bench("mul")
  call bench("div")
  call bench("sqrt")

  print '(a)', &
    " ---------------------------------------------------------------------------"
  print '(a,es11.3,a,es11.3,a,es11.3)', &
    " sinks: mf=", mf_sink, " fc=", fc_sink, " cw=", cw_sink

contains

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

  subroutine mf_drain(arr)
    !GCC$ ATTRIBUTES NOINLINE :: mf_drain
    type(float64x2), intent(inout) :: arr(:)
    real(dp) :: s
    integer :: j
    s = 0.0_dp
    do j = 1, size(fres); s = s + fres(j)%limbs(1); end do
    arr(1)%limbs(1) = arr(1)%limbs(1) + s * 1e-30_dp
    mf_sink = mf_sink + s
  end subroutine

  subroutine fc_drain(arr)
    !GCC$ ATTRIBUTES NOINLINE :: fc_drain
    type(dd_c), intent(inout) :: arr(:)
    real(dp) :: s
    integer :: j
    s = 0.0_dp
    do j = 1, size(dres); s = s + dres(j)%hi; end do
    arr(1)%hi = arr(1)%hi + s * 1e-30_dp
    fc_sink = fc_sink + s
  end subroutine

  subroutine cw_drain(arr)
    !GCC$ ATTRIBUTES NOINLINE :: cw_drain
    type(dd_c), intent(inout) :: arr(:)
    real(dp) :: s
    integer :: j
    s = 0.0_dp
    do j = 1, size(dres); s = s + dres(j)%hi; end do
    arr(1)%hi = arr(1)%hi + s * 1e-30_dp
    cw_sink = cw_sink + s
  end subroutine

  subroutine bench(op)
    character(*), intent(in) :: op
    integer :: r
    integer(int64) :: t0, n_ops
    real(real64) :: tmf, tfc, tcw, r1, r2, r3

    n_ops = int(N, int64) * int(REPS, int64)

    ! --- Fortran multifloats (Fortran ABI) ---
    call tick(t0)
    select case (op)
    case ("add");  do r=1,REPS; do i=1,N; fres(i)=f1(i)+f2(i);      end do; call mf_drain(f1); end do
    case ("sub");  do r=1,REPS; do i=1,N; fres(i)=f1(i)-f2(i);      end do; call mf_drain(f1); end do
    case ("mul");  do r=1,REPS; do i=1,N; fres(i)=f1(i)*f2(i);      end do; call mf_drain(f1); end do
    case ("div");  do r=1,REPS; do i=1,N; fres(i)=f1(i)/f2(i);      end do; call mf_drain(f1); end do
    case ("sqrt"); do r=1,REPS; do i=1,N; fres(i)=sqrt(f1(i));       end do; call mf_drain(f1); end do
    end select
    tmf = elapsed(t0)

    ! --- Fortran bind(c) kernels (C ABI, gfortran codegen) ---
    call tick(t0)
    select case (op)
    case ("add");  do r=1,REPS; do i=1,N; dres(i)=fc_add(d1(i),d2(i));  end do; call fc_drain(d1); end do
    case ("sub");  do r=1,REPS; do i=1,N; dres(i)=fc_sub(d1(i),d2(i));  end do; call fc_drain(d1); end do
    case ("mul");  do r=1,REPS; do i=1,N; dres(i)=fc_mul(d1(i),d2(i));  end do; call fc_drain(d1); end do
    case ("div");  do r=1,REPS; do i=1,N; dres(i)=fc_div(d1(i),d2(i));  end do; call fc_drain(d1); end do
    case ("sqrt"); do r=1,REPS; do i=1,N; dres(i)=fc_sqrt(d1(i));       end do; call fc_drain(d1); end do
    end select
    tfc = elapsed(t0)

    ! --- C wrapper (C ABI, g++ codegen) ---
    call tick(t0)
    select case (op)
    case ("add");  do r=1,REPS; do i=1,N; dres(i)=c_dd_add(d1(i),d2(i)); end do; call cw_drain(d1); end do
    case ("sub");  do r=1,REPS; do i=1,N; dres(i)=c_dd_sub(d1(i),d2(i)); end do; call cw_drain(d1); end do
    case ("mul");  do r=1,REPS; do i=1,N; dres(i)=c_dd_mul(d1(i),d2(i)); end do; call cw_drain(d1); end do
    case ("div");  do r=1,REPS; do i=1,N; dres(i)=c_dd_div(d1(i),d2(i)); end do; call cw_drain(d1); end do
    case ("sqrt"); do r=1,REPS; do i=1,N; dres(i)=c_dd_sqrt(d1(i));      end do; call cw_drain(d1); end do
    end select
    tcw = elapsed(t0)

    r1 = merge(tmf/tfc, 0.0_real64, tfc > 0)
    r2 = merge(tmf/tcw, 0.0_real64, tcw > 0)
    r3 = merge(tfc/tcw, 0.0_real64, tcw > 0)
    write(*, '(1x,a14,1x,i10,1x,f9.4,1x,f9.4,1x,f9.4,1x,f6.2,"x",1x,f6.2,"x",1x,f6.2,"x")') &
        op, n_ops, tmf, tfc, tcw, r1, r2, r3
  end subroutine

end program
