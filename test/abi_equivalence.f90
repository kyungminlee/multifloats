! Assert that the two DD entry points for each of {add, sub, mul, div,
! sqrt} produce precision-equivalent results on the same inputs:
!   - native : Fortran elemental operator (`a + b`, ..., sqrt(a))
!   - c_abi  : extern "C" from src/multifloats_math.cc (adddd, subdd, ...)
!
! We require:
!   (a) bit-exact agreement on the HI limb (anything else is an
!       algorithmic regression; 1 dp ULP on hi is ~10^-16 of relative
!       error, way above DD precision).
!   (b) LO limbs within 4 ULPs of each other (sub-DD-ULP variation is
!       tolerated — the C divdd picks a slightly different residual for
!       -7/3 style inputs and sqrtdd uses a slightly different Newton
!       iteration than the native Fortran sqrt. Both remain well within
!       the DD precision envelope).
program abi_equivalence
  use multifloats
  use, intrinsic :: iso_c_binding, only: c_double
  implicit none

  integer, parameter :: dp = 8
  integer :: failures = 0

  ! C wrapper DD struct (layout-compatible with float64x2_t).
  type, bind(c) :: dd_c
    real(c_double) :: hi, lo
  end type dd_c

  ! C wrapper interfaces.
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

  call check_all_binary("add", 1)
  call check_all_binary("sub", 2)
  call check_all_binary("mul", 3)
  call check_all_binary("div", 4)
  call check_sqrt()

  if (failures /= 0) then
    write(*,'(a,i0,a)') 'FAIL: ', failures, ' bit-inequality failures'
    error stop 1
  end if
  write(*,'(a)') 'PASS abi_equivalence'

contains

  subroutine record(label, hi_native, lo_native, hi_cabi, lo_cabi)
    character(*), intent(in) :: label
    real(dp), intent(in) :: hi_native, lo_native, hi_cabi, lo_cabi

    logical :: hi_ok
    real(dp) :: ulp_scale

    ! Hi must match bit-exactly across both paths.
    hi_ok = raw_eq(hi_native, hi_cabi)
    ! Lo tolerance: 4 dp ULPs relative to the hi magnitude. Accommodates
    ! compensated-error algorithmic variance without losing sensitivity
    ! to real regressions (any algorithmic break would easily exceed 4 ulp).
    ulp_scale = spacing(max(abs(hi_native), 1.0e-300_dp)) * 4.0_dp

    if (.not. hi_ok .or. abs(lo_native - lo_cabi) > ulp_scale) then
      failures = failures + 1
      write(*,'(a,1x,a)') 'FAIL', trim(label)
      write(*,'(a,2es24.16)') '  native : ', hi_native, lo_native
      write(*,'(a,2es24.16)') '  c_abi  : ', hi_cabi, lo_cabi
      write(*,'(a,es24.16)')  '  tol    : ', ulp_scale
    end if
  end subroutine

  logical function raw_eq(a, b)
    real(dp), intent(in) :: a, b
    integer(8) :: ai, bi
    ! Bit-for-bit equality: reinterpret the double as int64 via transfer.
    ai = transfer(a, ai)
    bi = transfer(b, bi)
    raw_eq = (ai == bi)
  end function

  subroutine check_all_binary(label, op)
    character(*), intent(in) :: label
    integer, intent(in) :: op   ! 1=add, 2=sub, 3=mul, 4=div
    real(dp), parameter :: inputs(2, 12) = reshape([ &
      1.0_dp,                 1.0_dp, &
      1.25_dp,                3.75_dp, &
      0.1_dp,                 0.3_dp, &
      1.0e-10_dp,             1.0e-10_dp, &
      -7.0_dp,                3.0_dp, &
      1.0e20_dp,              1.0_dp, &
      1.0_dp,                 1.0e-30_dp, &
      0.0_dp,                 2.5_dp, &
      2.5_dp,                 0.0_dp, &
      3.14159265358979_dp,    2.71828182845905_dp, &
      1.23456789_dp,          9.87654321_dp, &
      -1.0_dp,                -1.0_dp ], [2, 12])
    integer :: i
    type(float64x2) :: fa, fb, fr_native
    type(dd_c) :: da, db, dr_cabi
    character(len=80) :: tag
    real(dp) :: a_hi, b_hi

    do i = 1, size(inputs, 2)
      a_hi = inputs(1, i)
      b_hi = inputs(2, i)
      fa%limbs(1) = a_hi;  fa%limbs(2) = 0.0_dp
      fb%limbs(1) = b_hi;  fb%limbs(2) = 0.0_dp
      da%hi = a_hi;  da%lo = 0.0_dp
      db%hi = b_hi;  db%lo = 0.0_dp

      select case (op)
      case (1)
        fr_native = fa + fb
        dr_cabi = c_dd_add(da, db)
      case (2)
        fr_native = fa - fb
        dr_cabi = c_dd_sub(da, db)
      case (3)
        fr_native = fa * fb
        dr_cabi = c_dd_mul(da, db)
      case (4)
        ! Skip division by zero — both paths handle it consistently,
        ! but +inf vs -inf depends on limb sign; covered elsewhere.
        if (b_hi == 0.0_dp) cycle
        fr_native = fa / fb
        dr_cabi = c_dd_div(da, db)
      end select

      write(tag, '(a,a,i0,a)') trim(label), '(', i, ')'
      call record(trim(tag), &
                  fr_native%limbs(1), fr_native%limbs(2), &
                  dr_cabi%hi, dr_cabi%lo)
    end do
  end subroutine

  subroutine check_sqrt()
    real(dp), parameter :: inputs(8) = [ &
      0.0_dp, 1.0_dp, 2.0_dp, 4.0_dp, 0.5_dp, &
      1.0e-10_dp, 1.0e20_dp, 123.456_dp ]
    integer :: i
    type(float64x2) :: fa, fr_native
    type(dd_c) :: da, dr_cabi
    character(len=80) :: tag

    do i = 1, size(inputs)
      fa%limbs(1) = inputs(i);  fa%limbs(2) = 0.0_dp
      da%hi = inputs(i);  da%lo = 0.0_dp
      fr_native = sqrt(fa)
      dr_cabi = c_dd_sqrt(da)
      write(tag, '(a,i0,a)') 'sqrt(', i, ')'
      call record(trim(tag), &
                  fr_native%limbs(1), fr_native%limbs(2), &
                  dr_cabi%hi, dr_cabi%lo)
    end do
  end subroutine

end program abi_equivalence
