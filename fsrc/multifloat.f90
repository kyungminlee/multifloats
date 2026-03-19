module multifloat
  use iso_c_binding
  use, intrinsic :: ieee_arithmetic

  public :: float64x2

  type :: float64x2
    real*8 :: limbs(2)
  contains
    generic :: operator(*) => mul, mulf, fmul
    procedure :: mul, mulf
    procedure, pass(b) :: fmul
    generic :: operator(/) => div, divf, fdiv
    procedure :: div, divf
    procedure, pass(b) :: fdiv
    procedure :: isfinite => isfinite_f
  end type

  interface operator (+)
    module procedure add, add_fd, add_df
  end interface

  interface operator (-)
    module procedure sub, sub_fd, sub_df
  end interface

  interface operator (.lt.)
    module procedure lt, lt_fd, lt_df
  end interface

  interface operator (.gt.)
    module procedure gt, gt_fd, gt_df
  end interface

  interface operator (.le.)
    module procedure le, le_fd, le_df
  end interface

  interface operator (.ge.)
    module procedure ge, ge_fd, ge_df
  end interface

  interface assignment(=)
    module procedure assign_from_double
  end interface

  interface operator (.eq.)
    module procedure eq_ff, eq_fd, eq_df
  end interface

  interface operator (.ne.)
    module procedure ne_ff, ne_fd, ne_df
  end interface

contains

  elemental function isfinite_f(x) result(z)
    class(float64x2), intent(in) :: x
    logical :: z
    z = ieee_is_finite(x%limbs(1))
  end function

  elemental subroutine assign_from_double(lhs, rhs)
    type(float64x2), intent(out) :: lhs
    double precision, intent(in) :: rhs
    lhs%limbs(1) = rhs
    lhs%limbs(2) = 0.0d0
  end subroutine

  pure function eq_ff(x, y) result(z)
    type(float64x2), intent(in) :: x, y
    logical :: z
    type(float64x2) :: x2, y2
    x2 = x
    y2 = y
    call renormalize(x2)
    call renormalize(y2)
    z = (x2%limbs(1) .eq. y2%limbs(1)) .and. (x2%limbs(2) .eq. y2%limbs(2))
  end function

  pure function ne_ff(x, y) result(z)
    type(float64x2), intent(in) :: x, y
    logical :: z
    z = .not. eq_ff(x, y)
  end function

  pure function eq_fd(x, d) result(z)
    type(float64x2), intent(in) :: x
    double precision, intent(in) :: d
    logical :: z
    type(float64x2) :: x2
    x2 = x
    call renormalize(x2)
    z = (x2%limbs(1) .eq. d) .and. (x2%limbs(2) .eq. 0.0d0)
  end function

  pure function eq_df(d, x) result(z)
    double precision, intent(in) :: d
    type(float64x2), intent(in) :: x
    logical :: z
    z = eq_fd(x, d)
  end function

  pure function ne_fd(x, d) result(z)
    type(float64x2), intent(in) :: x
    double precision, intent(in) :: d
    logical :: z
    z = .not. eq_fd(x, d)
  end function

  pure function ne_df(d, x) result(z)
    double precision, intent(in) :: d
    type(float64x2), intent(in) :: x
    logical :: z
    z = .not. eq_fd(x, d)
  end function

  elemental function div(a, b) result(c)
    class(float64x2), intent(in) :: a, b
    type(float64x2) :: c
    double precision :: q1, q2, r
    if (b%limbs(1) == 0.0d0) then
      c%limbs(1) = a%limbs(1) / b%limbs(1)
      c%limbs(2) = 0.0d0
      return
    end if
    if (.not. ieee_is_finite(a%limbs(1)) .or. .not. ieee_is_finite(b%limbs(1))) then
      c%limbs(1) = a%limbs(1) / b%limbs(1)
      c%limbs(2) = 0.0d0
      return
    end if
    q1 = a%limbs(1) / b%limbs(1)
    r = ieee_fma(-q1, b%limbs(1), a%limbs(1))
    r = r + a%limbs(2) - q1 * b%limbs(2)
    q2 = r / b%limbs(1)
    c%limbs(1) = q1
    c%limbs(2) = q2
    call fast_two_sum(c%limbs(1), c%limbs(2))
  end function

  elemental function fdiv(a, b) result(c)
    double precision, intent(in) :: a
    class(float64x2), intent(in) :: b
    type(float64x2) :: c
    double precision :: q1, q2, r
    if (b%limbs(1) == 0.0d0) then
      c%limbs(1) = a / b%limbs(1)
      c%limbs(2) = 0.0d0
      return
    end if
    if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(b%limbs(1))) then
      c%limbs(1) = a / b%limbs(1)
      c%limbs(2) = 0.0d0
      return
    end if
    q1 = a / b%limbs(1)
    r = ieee_fma(-q1, b%limbs(1), a)
    r = r - q1 * b%limbs(2)
    q2 = r / b%limbs(1)
    c%limbs(1) = q1
    c%limbs(2) = q2
    call fast_two_sum(c%limbs(1), c%limbs(2))
  end function

  elemental function divf(a, b) result(c)
    class(float64x2), intent(in) :: a
    double precision, intent(in) :: b
    type(float64x2) :: c
    double precision :: q1, q2, r
    if (b == 0.0d0) then
      c%limbs(1) = a%limbs(1) / b
      c%limbs(2) = 0.0d0
      return
    end if
    if (.not. ieee_is_finite(a%limbs(1)) .or. .not. ieee_is_finite(b)) then
      c%limbs(1) = a%limbs(1) / b
      c%limbs(2) = 0.0d0
      return
    end if
    q1 = a%limbs(1) / b
    r = ieee_fma(-q1, b, a%limbs(1))
    r = r + a%limbs(2)
    q2 = r / b
    c%limbs(1) = q1
    c%limbs(2) = q2
    call fast_two_sum(c%limbs(1), c%limbs(2))
  end function

  elemental function mul(a, b) result(c)
    class(float64x2), intent(in) :: a, b
    type(float64x2) :: c
    double precision :: p, e
    if (.not. ieee_is_finite(a%limbs(1)) .or. .not. ieee_is_finite(b%limbs(1))) then
      c%limbs(1) = a%limbs(1) * b%limbs(1)
      c%limbs(2) = 0.0d0
      return
    end if
    call two_prod(a%limbs(1), b%limbs(1), p, e)
    e = e + a%limbs(1) * b%limbs(2)
    e = e + a%limbs(2) * b%limbs(1)
    e = e + a%limbs(2) * b%limbs(2)
    c%limbs(1) = p
    c%limbs(2) = e
    call fast_two_sum(c%limbs(1), c%limbs(2))
  end function

  elemental function fmul(a, b) result(c)
    double precision, intent(in) :: a
    class(float64x2), intent(in) :: b
    type(float64x2) :: c
    double precision :: p, e
    if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(b%limbs(1))) then
      c%limbs(1) = a * b%limbs(1)
      c%limbs(2) = 0.0d0
      return
    end if
    call two_prod(a, b%limbs(1), p, e)
    e = e + a * b%limbs(2)
    c%limbs(1) = p
    c%limbs(2) = e
    call fast_two_sum(c%limbs(1), c%limbs(2))
  end function

  elemental function mulf(a, b) result(c)
    class(float64x2), intent(in) :: a
    double precision, intent(in) :: b
    type(float64x2) :: c
    double precision :: p, e
    if (.not. ieee_is_finite(a%limbs(1)) .or. .not. ieee_is_finite(b)) then
      c%limbs(1) = a%limbs(1) * b
      c%limbs(2) = 0.0d0
      return
    end if
    call two_prod(a%limbs(1), b, p, e)
    e = e + a%limbs(2) * b
    c%limbs(1) = p
    c%limbs(2) = e
    call fast_two_sum(c%limbs(1), c%limbs(2))
  end function

  elemental subroutine fast_two_sum(a, b)
    double precision, intent(inout):: a, b
    double precision :: s, b_prime, b_err
    s = a + b
    b_prime = s - a
    b_err = b - b_prime
    a = s
    b = b_err
  end subroutine

  elemental subroutine two_sum(a, b)
    double precision, intent(inout):: a, b
    double precision :: s, a_prime, b_prime, a_err, b_err
    s = a + b
    a_prime = s - b
    b_prime = s - a_prime
    a_err = a - a_prime
    b_err = b - b_prime
    a = s
    b = a_err + b_err
  end subroutine

  elemental subroutine two_prod(a, b, p, e)
    double precision, intent(in) :: a, b
    double precision, intent(out) :: p, e
    p = a * b
    e = ieee_fma(a, b, -p)
  end subroutine

  elemental subroutine renormalize(x)
    type(float64x2), intent(inout) :: x
    call two_sum(x%limbs(1), x%limbs(2))
  end subroutine

  elemental function add(x, y) result(z)
    type(float64x2), intent(in) :: x, y
    type(float64x2) :: z
    double precision :: s, e
    if (.not. ieee_is_finite(x%limbs(1)) .or. .not. ieee_is_finite(y%limbs(1))) then
      z%limbs(1) = x%limbs(1) + y%limbs(1)
      z%limbs(2) = 0.0d0
      return
    end if
    s = x%limbs(1)
    e = y%limbs(1)
    call two_sum(s, e)
    e = e + x%limbs(2) + y%limbs(2)
    z%limbs(1) = s
    z%limbs(2) = e
    call fast_two_sum(z%limbs(1), z%limbs(2))
  end function

  elemental function add_fd(x, d) result(z)
    type(float64x2), intent(in) :: x
    double precision, intent(in) :: d
    type(float64x2) :: z
    double precision :: s, e
    if (.not. ieee_is_finite(x%limbs(1)) .or. .not. ieee_is_finite(d)) then
      z%limbs(1) = x%limbs(1) + d
      z%limbs(2) = 0.0d0
      return
    end if
    s = x%limbs(1)
    e = d
    call two_sum(s, e)
    e = e + x%limbs(2)
    z%limbs(1) = s
    z%limbs(2) = e
    call fast_two_sum(z%limbs(1), z%limbs(2))
  end function

  elemental function add_df(d, x) result(z)
    double precision, intent(in) :: d
    type(float64x2), intent(in) :: x
    type(float64x2) :: z
    z = add_fd(x, d)
  end function

  elemental function sub(x, y) result(z)
    type(float64x2), intent(in) :: x, y
    type(float64x2) :: z
    double precision :: s, e
    if (.not. ieee_is_finite(x%limbs(1)) .or. .not. ieee_is_finite(y%limbs(1))) then
      z%limbs(1) = x%limbs(1) - y%limbs(1)
      z%limbs(2) = 0.0d0
      return
    end if
    s = x%limbs(1)
    e = -y%limbs(1)
    call two_sum(s, e)
    e = e + x%limbs(2) - y%limbs(2)
    z%limbs(1) = s
    z%limbs(2) = e
    call fast_two_sum(z%limbs(1), z%limbs(2))
  end function

  elemental function sub_fd(x, d) result(z)
    type(float64x2), intent(in) :: x
    double precision, intent(in) :: d
    type(float64x2) :: z
    double precision :: s, e
    if (.not. ieee_is_finite(x%limbs(1)) .or. .not. ieee_is_finite(d)) then
      z%limbs(1) = x%limbs(1) - d
      z%limbs(2) = 0.0d0
      return
    end if
    s = x%limbs(1)
    e = -d
    call two_sum(s, e)
    e = e + x%limbs(2)
    z%limbs(1) = s
    z%limbs(2) = e
    call fast_two_sum(z%limbs(1), z%limbs(2))
  end function

  elemental function sub_df(d, x) result(z)
    double precision, intent(in) :: d
    type(float64x2), intent(in) :: x
    type(float64x2) :: z
    double precision :: s, e
    if (.not. ieee_is_finite(d) .or. .not. ieee_is_finite(x%limbs(1))) then
      z%limbs(1) = d - x%limbs(1)
      z%limbs(2) = 0.0d0
      return
    end if
    s = d
    e = -x%limbs(1)
    call two_sum(s, e)
    e = e - x%limbs(2)
    z%limbs(1) = s
    z%limbs(2) = e
    call fast_two_sum(z%limbs(1), z%limbs(2))
  end function

  pure function lt(x, y) result(z)
    type(float64x2), intent(in) :: x, y
    type(float64x2) :: x2, y2
    logical :: z
    x2 = x
    y2 = y
    call renormalize(x2)
    call renormalize(y2)
    z = (x2%limbs(1) .lt. y2%limbs(1)) .or.     &
        & ((x2%limbs(1) .eq. y2%limbs(1)) .and. &
        &  (x2%limbs(2) .lt. y2%limbs(2)))
  end function

  pure function lt_fd(x, d) result(z)
    type(float64x2), intent(in) :: x
    double precision, intent(in) :: d
    logical :: z
    type(float64x2) :: x2
    x2 = x
    call renormalize(x2)
    z = (x2%limbs(1) .lt. d) .or. &
        & ((x2%limbs(1) .eq. d) .and. (x2%limbs(2) .lt. 0.0d0))
  end function

  pure function lt_df(d, x) result(z)
    double precision, intent(in) :: d
    type(float64x2), intent(in) :: x
    logical :: z
    type(float64x2) :: x2
    x2 = x
    call renormalize(x2)
    z = (d .lt. x2%limbs(1)) .or. &
        & ((d .eq. x2%limbs(1)) .and. (0.0d0 .lt. x2%limbs(2)))
  end function

  pure function gt(x, y) result(z)
    type(float64x2), intent(in) :: x, y
    logical :: z
    z = lt(y, x)
  end function

  pure function gt_fd(x, d) result(z)
    type(float64x2), intent(in) :: x
    double precision, intent(in) :: d
    logical :: z
    z = lt_df(d, x)
  end function

  pure function gt_df(d, x) result(z)
    double precision, intent(in) :: d
    type(float64x2), intent(in) :: x
    logical :: z
    z = lt_fd(x, d)
  end function

  pure function le(x, y) result(z)
    type(float64x2), intent(in) :: x, y
    logical :: z
    z = .not. lt(y, x)
  end function

  pure function le_fd(x, d) result(z)
    type(float64x2), intent(in) :: x
    double precision, intent(in) :: d
    logical :: z
    z = .not. lt_df(d, x)
  end function

  pure function le_df(d, x) result(z)
    double precision, intent(in) :: d
    type(float64x2), intent(in) :: x
    logical :: z
    z = .not. lt_fd(x, d)
  end function

  pure function ge(x, y) result(z)
    type(float64x2), intent(in) :: x, y
    logical :: z
    z = .not. lt(x, y)
  end function

  pure function ge_fd(x, d) result(z)
    type(float64x2), intent(in) :: x
    double precision, intent(in) :: d
    logical :: z
    z = .not. lt_fd(x, d)
  end function

  pure function ge_df(d, x) result(z)
    double precision, intent(in) :: d
    type(float64x2), intent(in) :: x
    logical :: z
    z = .not. lt_df(d, x)
  end function

end module
