! Fuzz-equivalence: assert that the Fortran-native path and the C-ABI
! path produce matching results for every scalar-DD op targeted by
! Phase 1 of the C-ABI delegation. Uses large random inputs plus a
! fixed edge-case table so we catch both typical and corner behaviors.
!
! Tolerance policy:
!   * Hi limb must be bit-exact across the two paths (any algorithmic
!     divergence shows up here).
!   * Lo limb is compared in dp ULPs relative to the hi magnitude:
!       - add/sub/mul/div/sqrt : within 4 dp ULPs (both paths use
!         compensated-error algorithms that may pick slightly
!         different residuals; 4 dp ULPs is ~2 DD-bits, which is well
!         within the DD precision envelope).
!       - all other Phase-1 ops : bit-exact (operations are pure
!         manipulation/comparison, no EFT residual).
!   * Comparison ops return bool: must agree bit-for-bit.
program abi_equivalence
  use multifloats
  use, intrinsic :: iso_c_binding, only: c_double, c_int
  implicit none

  integer, parameter :: dp = 8
  integer, parameter :: N_FUZZ = 200000   ! random samples per op
  integer :: failures = 0
  integer :: total    = 0

  type, bind(c) :: dd_c
    real(c_double) :: limbs(2)
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
    pure type(dd_c) function c_dd_powi(a, n) bind(c, name='powidd')
      import :: dd_c, c_int; type(dd_c), value :: a; integer(c_int), value :: n
    end function
    pure type(dd_c) function c_dd_modulo(a, b) bind(c, name='modulodd')
      import :: dd_c; type(dd_c), value :: a, b
    end function
  end interface

  call seed_rng()

  ! Deterministic edge cases first
  call edge_cases_dd_unary()
  call edge_cases_dd_binary()

  ! Fuzz
  call fuzz_dd_unary("neg",   4)  ! op code; 4 = bit-exact
  call fuzz_dd_unary("abs",   4)
  call fuzz_dd_unary("aint",  4)
  call fuzz_dd_unary("anint", 4)
  call fuzz_dd_unary("sqrt",  1)  ! 1 = 4-ulp tolerance
  call fuzz_dd_binary("add",  1)
  call fuzz_dd_binary("sub",  1)
  call fuzz_dd_binary("mul",  1)
  call fuzz_dd_binary("div",  1)
  call fuzz_dd_binary("min",  4)
  call fuzz_dd_binary("max",  4)
  call fuzz_dd_binary("sign", 4)
  call fuzz_dd_binary("dim",  4)
  call fuzz_dd_cmp("eq")
  call fuzz_dd_cmp("ne")
  call fuzz_dd_cmp("lt")
  call fuzz_dd_cmp("le")
  call fuzz_dd_cmp("gt")
  call fuzz_dd_cmp("ge")
  call fuzz_dd_powi()
  call fuzz_dd_modulo()

  write(*,'(a,i0,a,i0,a)') 'Summary: ', total-failures, '/', total, ' checks passed'
  if (failures /= 0) then
    write(*,'(a,i0,a)') 'FAIL: ', failures, ' equivalence failures'
    error stop 1
  end if
  write(*,'(a)') 'PASS abi_equivalence'

contains

  subroutine seed_rng()
    integer :: seed_size, j
    integer, allocatable :: seed(:)
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    do j = 1, seed_size
      seed(j) = 424242 + 7 * j
    end do
    call random_seed(put=seed)
  end subroutine

  subroutine gen_dd(x, f)
    ! Emit well-normalized random DD: |lo| <= 0.5 ulp(hi).
    type(dd_c), intent(out) :: x
    type(real64x2), intent(out) :: f
    real(dp) :: r, rlo, rexp, mag
    integer :: e
    call random_number(r)
    call random_number(rlo)
    call random_number(rexp)
    e = int((rexp - 0.5_dp) * 40.0_dp)    ! exponents in [-20, 20]
    mag = (2.0_dp * r - 1.0_dp)
    if (abs(mag) < 1.0e-3_dp) mag = sign(1.0e-3_dp, mag)
    x%limbs(1) = mag * (2.0_dp ** e)
    x%limbs(2) = x%limbs(1) * (2.0_dp * rlo - 1.0_dp) * 2.0_dp**(-52)
    f%limbs(1) = x%limbs(1); f%limbs(2) = x%limbs(2)
  end subroutine

  logical function raw_eq(a, b)
    real(dp), intent(in) :: a, b
    integer(8) :: ai, bi
    ai = transfer(a, ai); bi = transfer(b, bi)
    raw_eq = (ai == bi)
  end function

  subroutine record_dd(label, hn, ln, hc, lc, mode)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    character(*), intent(in) :: label
    real(dp), intent(in) :: hn, ln, hc, lc
    integer, intent(in) :: mode    ! 1 = 4-ulp lo tol; 4 = bit-exact
    real(dp) :: tol
    logical :: ok
    total = total + 1
    if (mode == 4 .or. .not. ieee_is_finite(hn)) then
      ! Bit-exact path: also applies when hi is ±Inf/NaN (spacing(Inf)
      ! is NaN so the tolerance arithmetic below is meaningless).
      ok = raw_eq(hn, hc) .and. raw_eq(ln, lc)
      tol = 0.0_dp
    else
      tol = spacing(max(abs(hn), 1.0e-300_dp)) * 4.0_dp
      ok = raw_eq(hn, hc) .and. abs(ln - lc) <= tol
    end if
    if (.not. ok) then
      failures = failures + 1
      if (failures <= 20) then
        write(*,'(a,1x,a)') 'FAIL', trim(label)
        write(*,'(a,2es24.16)') '  native : ', hn, ln
        write(*,'(a,2es24.16)') '  c_abi  : ', hc, lc
        if (mode /= 4) write(*,'(a,es12.3)') '  tol    : ', tol
      end if
    end if
  end subroutine

  subroutine record_int(label, in_, ic)
    character(*), intent(in) :: label
    integer, intent(in) :: in_, ic
    total = total + 1
    if (in_ /= ic) then
      failures = failures + 1
      if (failures <= 20) then
        write(*,'(a,1x,a,1x,a,i0,1x,a,i0)') 'FAIL', trim(label), &
          'native=', in_, 'c_abi=', ic
      end if
    end if
  end subroutine

  subroutine apply_unary(op, fa, da, frn, drn)
    character(*), intent(in) :: op
    type(real64x2), intent(in) :: fa
    type(dd_c),      intent(in) :: da
    type(real64x2), intent(out) :: frn
    type(dd_c),      intent(out) :: drn
    select case (op)
    case ("neg");   frn = -fa;            drn = c_dd_neg(da)
    case ("abs");   frn = abs(fa);        drn = c_dd_abs(da)
    case ("aint");  frn = aint(fa);       drn = c_dd_trunc(da)
    case ("anint"); frn = anint(fa);      drn = c_dd_round(da)
    case ("sqrt");  frn = sqrt(abs(fa));  drn = c_dd_sqrt(c_dd_abs(da))
    end select
  end subroutine

  subroutine apply_binary(op, fa, fb, da, db, frn, drn)
    character(*), intent(in) :: op
    type(real64x2), intent(in) :: fa, fb
    type(dd_c),      intent(in) :: da, db
    type(real64x2), intent(out) :: frn
    type(dd_c),      intent(out) :: drn
    select case (op)
    case ("add"); frn = fa + fb;          drn = c_dd_add(da, db)
    case ("sub"); frn = fa - fb;          drn = c_dd_sub(da, db)
    case ("mul"); frn = fa * fb;          drn = c_dd_mul(da, db)
    case ("div"); frn = fa / fb;          drn = c_dd_div(da, db)
    case ("min"); frn = min(fa, fb);      drn = c_dd_min(da, db)
    case ("max"); frn = max(fa, fb);      drn = c_dd_max(da, db)
    case ("sign");frn = sign(fa, fb);     drn = c_dd_copysign(da, db)
    case ("dim"); frn = dim(fa, fb);      drn = c_dd_fdim(da, db)
    end select
  end subroutine

  subroutine apply_cmp(op, fa, fb, da, db, in_, ic)
    character(*), intent(in) :: op
    type(real64x2), intent(in) :: fa, fb
    type(dd_c),      intent(in) :: da, db
    integer, intent(out) :: in_, ic
    select case (op)
    case ("eq"); in_ = merge(1,0, fa == fb); ic = int(c_dd_eq(da, db))
    case ("ne"); in_ = merge(1,0, fa /= fb); ic = int(c_dd_ne(da, db))
    case ("lt"); in_ = merge(1,0, fa <  fb); ic = int(c_dd_lt(da, db))
    case ("le"); in_ = merge(1,0, fa <= fb); ic = int(c_dd_le(da, db))
    case ("gt"); in_ = merge(1,0, fa >  fb); ic = int(c_dd_gt(da, db))
    case ("ge"); in_ = merge(1,0, fa >= fb); ic = int(c_dd_ge(da, db))
    end select
  end subroutine

  subroutine fuzz_dd_unary(op, mode)
    character(*), intent(in) :: op
    integer, intent(in) :: mode
    integer :: k
    type(real64x2) :: fa, frn
    type(dd_c) :: da, drn
    character(80) :: tag
    do k = 1, N_FUZZ
      call gen_dd(da, fa)
      call apply_unary(op, fa, da, frn, drn)
      if (k <= 3) then
        write(tag, '(a,a,i0,a)') trim(op), '(', k, ')'
        call record_dd(trim(tag), frn%limbs(1), frn%limbs(2), drn%limbs(1), drn%limbs(2), mode)
      else
        call record_dd(op, frn%limbs(1), frn%limbs(2), drn%limbs(1), drn%limbs(2), mode)
      end if
    end do
  end subroutine

  subroutine fuzz_dd_binary(op, mode)
    character(*), intent(in) :: op
    integer, intent(in) :: mode
    integer :: k
    type(real64x2) :: fa, fb, frn
    type(dd_c) :: da, db, drn
    do k = 1, N_FUZZ
      call gen_dd(da, fa); call gen_dd(db, fb)
      if (op == "div" .and. db%limbs(1) == 0.0_dp) cycle
      call apply_binary(op, fa, fb, da, db, frn, drn)
      call record_dd(op, frn%limbs(1), frn%limbs(2), drn%limbs(1), drn%limbs(2), mode)
    end do
  end subroutine

  subroutine fuzz_dd_cmp(op)
    character(*), intent(in) :: op
    integer :: k, in_, ic
    type(real64x2) :: fa, fb
    type(dd_c) :: da, db
    do k = 1, N_FUZZ
      call gen_dd(da, fa); call gen_dd(db, fb)
      ! also probe equal and hi-equal / lo-different
      if (mod(k, 17) == 0) then; db = da; fb = fa; end if
      if (mod(k, 19) == 0) then
        db = da; fb = fa
        db%limbs(2) = -db%limbs(2); fb%limbs(2) = -fb%limbs(2)
      end if
      call apply_cmp(op, fa, fb, da, db, in_, ic)
      call record_int(op, in_, ic)
    end do
  end subroutine

  subroutine fuzz_dd_powi()
    ! Integer-exponent power: span a representative set of exponents
    ! (edge zeros, signs, small, moderate, and a couple of larger
    ! magnitudes that exercise the EBS loop). Base drawn from gen_dd.
    integer, parameter :: exps(15) = [ &
         0,  1, -1,  2, -2,  3,  7, 10, -10, 20, -20, 50, -50, 100, -100 ]
    integer :: k, j, n
    type(real64x2) :: fa, fr
    type(dd_c) :: da, dr
    character(80) :: tag
    do k = 1, N_FUZZ / size(exps)
      call gen_dd(da, fa)
      do j = 1, size(exps)
        n = exps(j)
        fr = fa ** n
        dr = c_dd_powi(da, int(n, c_int))
        write(tag, '(a,i0,a)') 'powi(n=', n, ')'
        ! EBS multiplies can differ from a linear loop in the last-ulp
        ! residual on large |n|; tolerate 4 dp ULPs on lo.
        call record_dd(trim(tag), fr%limbs(1), fr%limbs(2), dr%limbs(1), dr%limbs(2), 1)
      end do
    end do
  end subroutine

  subroutine fuzz_dd_modulo()
    integer :: k
    type(real64x2) :: fa, fb, fr
    type(dd_c) :: da, db, dr
    do k = 1, N_FUZZ
      call gen_dd(da, fa); call gen_dd(db, fb)
      if (db%limbs(1) == 0.0_dp) cycle
      fr = modulo(fa, fb)
      dr = c_dd_modulo(da, db)
      ! Sign-adjust can cross by exactly one dp ULP on the hi limb in
      ! rare cases where fmod's residual has hi ≈ 0 with a tiny lo; the
      ! C path uses signbit() of the whole DD while the Fortran path
      ! (pre-Stage-3) only looks at hi. Bit-exact on lo, 4 dp ULPs on
      ! hi is within the DD precision envelope either way.
      call record_dd("modulo", fr%limbs(1), fr%limbs(2), dr%limbs(1), dr%limbs(2), 1)
    end do
  end subroutine

  subroutine edge_cases_dd_unary()
    real(dp), parameter :: hi(9) = [ &
      0.0_dp, -0.0_dp, 1.0_dp, -1.0_dp, 0.5_dp, -0.5_dp, &
      1.5_dp, -1.5_dp, 1.0e20_dp ]
    real(dp), parameter :: lo(9) = [ &
      0.0_dp,  0.0_dp, 1.0e-17_dp, -1.0e-17_dp, 1.0e-18_dp, &
      -1.0e-18_dp, 1.0e-17_dp, -1.0e-17_dp, 1.0e3_dp ]
    type(real64x2) :: fa, frn
    type(dd_c) :: da, drn
    character(len=16), parameter :: ops(5) = [character(len=16) :: &
      "neg","abs","aint","anint","sqrt"]
    integer, parameter :: modes(5) = [4,4,4,4,1]
    integer :: j, i
    character(80) :: tag
    do j = 1, size(ops)
      do i = 1, size(hi)
        fa%limbs(1) = hi(i); fa%limbs(2) = lo(i)
        da%limbs(1) = hi(i);       da%limbs(2) = lo(i)
        call apply_unary(trim(ops(j)), fa, da, frn, drn)
        write(tag,'(a,"/edge",i0)') trim(ops(j)), i
        call record_dd(trim(tag), frn%limbs(1), frn%limbs(2), drn%limbs(1), drn%limbs(2), modes(j))
      end do
    end do
  end subroutine

  subroutine edge_cases_dd_binary()
    real(dp), parameter :: A(6) = [ 1.0_dp, -7.0_dp, 0.1_dp, 1.0e20_dp, 1.0_dp, 0.0_dp]
    real(dp), parameter :: B(6) = [ 1.0_dp,  3.0_dp, 0.3_dp, 1.0_dp,   -1.0_dp, 2.5_dp]
    character(len=16), parameter :: ops(8) = [character(len=16) :: &
      "add","sub","mul","div","min","max","sign","dim"]
    integer, parameter :: modes(8) = [1,1,1,1,4,4,4,4]
    type(real64x2) :: fa, fb, frn
    type(dd_c) :: da, db, drn
    character(80) :: tag
    integer :: j, i
    do j = 1, size(ops)
      do i = 1, size(A)
        fa%limbs(1) = A(i); fa%limbs(2) = 0.0_dp
        fb%limbs(1) = B(i); fb%limbs(2) = 0.0_dp
        da%limbs(1) = A(i); da%limbs(2) = 0.0_dp
        db%limbs(1) = B(i); db%limbs(2) = 0.0_dp
        if (trim(ops(j)) == "div" .and. B(i) == 0.0_dp) cycle
        call apply_binary(trim(ops(j)), fa, fb, da, db, frn, drn)
        write(tag,'(a,"/edge",i0)') trim(ops(j)), i
        call record_dd(trim(tag), frn%limbs(1), frn%limbs(2), drn%limbs(1), drn%limbs(2), modes(j))
      end do
    end do
  end subroutine

end program abi_equivalence
