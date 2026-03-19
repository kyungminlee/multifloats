module multifloat
  use iso_c_binding
  use, intrinsic :: ieee_arithmetic

  public :: float64x2
  public :: abs, sqrt, sign, min, max, floor, ceiling, aint, anint, exp, log, log10

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
    procedure :: precision => precision_f
    procedure :: minexponent => minexponent_f
    procedure :: maxexponent => maxexponent_f
    procedure :: tiny => tiny_f
    procedure :: huge => huge_f
    procedure :: exponent => exponent_f
    procedure :: scale => scale_f
    procedure :: floor => floor_f
    procedure :: ceiling => ceiling_f
    procedure :: aint => aint_f
    procedure :: anint => anint_f
    procedure :: abs => abs_f
    procedure :: sqrt => sqrt_f
    procedure :: fraction => fraction_f
    procedure :: set_exponent => set_exponent_f
    procedure :: spacing => spacing_f
    procedure :: rrspacing => rrspacing_f
    procedure :: exp => exp_f
    procedure :: log => log_f
    procedure :: log10 => log10_f
  end type

  interface operator (+)
    module procedure add, add_fd, add_df
  end interface

  interface operator (-)
    module procedure sub, sub_fd, sub_df, neg
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

  interface abs
    module procedure abs_f
  end interface

  interface sqrt
    module procedure sqrt_f
  end interface

  interface sign
    module procedure sign_f, sign_fd, sign_df
  end interface

  interface min
    module procedure min_ff, min_fd, min_df
  end interface

  interface max
    module procedure max_ff, max_fd, max_df
  end interface

  interface aint
    module procedure aint_f
  end interface

  interface anint
    module procedure anint_f
  end interface

  interface floor
    module procedure floor_f
  end interface

  interface ceiling
    module procedure ceiling_f
  end interface

  interface exp
    module procedure exp_f
  end interface

  interface log
    module procedure log_f
  end interface

  interface log10
    module procedure log10_f
  end interface

contains

  elemental function exp_f(x) result(z)
    class(float64x2), intent(in) :: x
    type(float64x2) :: z
    double precision :: y
    if (.not. ieee_is_finite(x%limbs(1))) then
       z%limbs(1) = exp(x%limbs(1))
       z%limbs(2) = 0.0d0
       return
    end if
    y = exp(x%limbs(1))
    ! exp(x_h + x_l) = exp(x_h) * exp(x_l) approx exp(x_h) * (1 + x_l)
    z = y + to_f64x2_d(y) * x%limbs(2)
  end function

  elemental function log_f(x) result(z)
    class(float64x2), intent(in) :: x
    type(float64x2) :: z
    double precision :: y
    if (.not. ieee_is_finite(x%limbs(1)) .or. x%limbs(1) <= 0.0d0) then
       z%limbs(1) = log(x%limbs(1))
       z%limbs(2) = 0.0d0
       return
    end if
    ! Initial guess
    y = log(x%limbs(1))
    ! One NR step: z = y + (x - exp(y)) / exp(y)
    ! which is z = y + x * exp(-y) - 1
    ! z = y + x * exp_f(-to_f64x2_d(y)) - 1.0d0
    z = to_f64x2_d(y) + to_f64x2_d(x%limbs(2)) / to_f64x2_d(x%limbs(1))
  end function

  elemental function log10_f(x) result(z)
    class(float64x2), intent(in) :: x
    type(float64x2) :: z
    z = log_f(x) / log_f(to_f64x2_d(10.0d0))
  end function

  elemental function to_f64x2_d(d) result(z)
    double precision, intent(in) :: d
    type(float64x2) :: z
    z%limbs(1) = d
    z%limbs(2) = 0.0d0
  end function

  elemental function abs_f(x) result(z)
    class(float64x2), intent(in) :: x
    type(float64x2) :: z
    if (x%limbs(1) < 0.0d0) then
      z = -x
    else
      z = x
    end if
  end function

  elemental function sqrt_f(x) result(z)
    class(float64x2), intent(in) :: x
    type(float64x2) :: z
    double precision :: s, p, e, h
    if (x%limbs(1) <= 0.0d0) then
       if (x%limbs(1) < 0.0d0) then
          z%limbs(1) = ieee_value(0.0d0, ieee_quiet_nan)
          z%limbs(2) = 0.0d0
       else
          z%limbs = 0.0d0
       end if
       return
    end if
    s = sqrt(x%limbs(1))
    ! Refinement: z = s + (x - s*s) / (2*s)
    call two_prod(s, s, p, e)
    h = (x%limbs(1) - p) + (x%limbs(2) - e)
    ! Must use double-double addition
    z%limbs(1) = s
    z%limbs(2) = h / (2.0d0 * s)
    call renormalize(z)
  end function

  elemental function sign_f(a, b) result(z)
    class(float64x2), intent(in) :: a, b
    type(float64x2) :: z
    z = abs(a)
    if (b%limbs(1) < 0.0d0) z = -z
  end function

  elemental function sign_fd(a, b) result(z)
    class(float64x2), intent(in) :: a
    double precision, intent(in) :: b
    type(float64x2) :: z
    z = abs(a)
    if (b < 0.0d0) z = -z
  end function

  elemental function sign_df(a, b) result(z)
    double precision, intent(in) :: a
    class(float64x2), intent(in) :: b
    type(float64x2) :: z
    z%limbs(1) = abs(a)
    z%limbs(2) = 0.0d0
    if (b%limbs(1) < 0.0d0) z = -z
  end function

  elemental function min_ff(a, b) result(z)
    class(float64x2), intent(in) :: a, b
    type(float64x2) :: z
    if (a < b) then
      z = a
    else
      z = b
    end if
  end function

  elemental function min_fd(a, b) result(z)
    class(float64x2), intent(in) :: a
    double precision, intent(in) :: b
    type(float64x2) :: z
    if (a < b) then
      z = a
    else
      z = b
    end if
  end function

  elemental function min_df(a, b) result(z)
    double precision, intent(in) :: a
    class(float64x2), intent(in) :: b
    type(float64x2) :: z
    if (a < b) then
      z = a
    else
      z = b
    end if
  end function

  elemental function max_ff(a, b) result(z)
    class(float64x2), intent(in) :: a, b
    type(float64x2) :: z
    if (a > b) then
      z = a
    else
      z = b
    end if
  end function

  elemental function max_fd(a, b) result(z)
    class(float64x2), intent(in) :: a
    double precision, intent(in) :: b
    type(float64x2) :: z
    if (a > b) then
      z = a
    else
      z = b
    end if
  end function

  elemental function max_df(a, b) result(z)
    double precision, intent(in) :: a
    class(float64x2), intent(in) :: b
    type(float64x2) :: z
    if (a > b) then
      z = a
    else
      z = b
    end if
  end function

  elemental function aint_f(x) result(z)
    class(float64x2), intent(in) :: x
    type(float64x2) :: z
    z%limbs(1) = aint(x%limbs(1))
    if (z%limbs(1) == x%limbs(1)) then
      z%limbs(2) = aint(x%limbs(2))
    else
      z%limbs(2) = 0.0d0
    end if
  end function

  elemental function anint_f(x) result(z)
    class(float64x2), intent(in) :: x
    type(float64x2) :: z
    z%limbs(1) = anint(x%limbs(1))
    if (z%limbs(1) == x%limbs(1)) then
      z%limbs(2) = anint(x%limbs(2))
    else
      z%limbs(2) = 0.0d0
    end if
  end function

  elemental function floor_f(x) result(z)
    class(float64x2), intent(in) :: x
    type(float64x2) :: z
    z%limbs(1) = floor(x%limbs(1))
    if (z%limbs(1) == x%limbs(1)) then
      z%limbs(2) = floor(x%limbs(2))
    else
      z%limbs(2) = 0.0d0
    end if
  end function

  elemental function ceiling_f(x) result(z)
    class(float64x2), intent(in) :: x
    type(float64x2) :: z
    z%limbs(1) = ceiling(x%limbs(1))
    if (z%limbs(1) == x%limbs(1)) then
      z%limbs(2) = ceiling(x%limbs(2))
    else
      z%limbs(2) = 0.0d0
    end if
  end function

  elemental function precision_f(x) result(z)
    class(float64x2), intent(in) :: x
    integer :: z
    z = 31
  end function

  elemental function minexponent_f(x) result(z)
    class(float64x2), intent(in) :: x
    integer :: z
    z = minexponent(1.0d0)
  end function

  elemental function maxexponent_f(x) result(z)
    class(float64x2), intent(in) :: x
    integer :: z
    z = maxexponent(1.0d0)
  end function

  elemental function tiny_f(x) result(z)
    class(float64x2), intent(in) :: x
    type(float64x2) :: z
    z%limbs(1) = tiny(1.0d0)
    z%limbs(2) = 0.0d0
  end function

  elemental function huge_f(x) result(z)
    class(float64x2), intent(in) :: x
    type(float64x2) :: z
    z%limbs(1) = huge(1.0d0)
    z%limbs(2) = 0.0d0
  end function

  elemental function exponent_f(x) result(z)
    class(float64x2), intent(in) :: x
    integer :: z
    z = exponent(x%limbs(1))
  end function

  elemental function scale_f(x, i) result(z)
    class(float64x2), intent(in) :: x
    integer, intent(in) :: i
    type(float64x2) :: z
    z%limbs = scale(x%limbs, i)
  end function

  elemental function fraction_f(x) result(z)
    class(float64x2), intent(in) :: x
    type(float64x2) :: z
    z = scale_f(x, -exponent(x%limbs(1)))
  end function

  elemental function set_exponent_f(x, i) result(z)
    class(float64x2), intent(in) :: x
    integer, intent(in) :: i
    type(float64x2) :: z
    z = scale_f(fraction_f(x), i)
  end function

  elemental function spacing_f(x) result(z)
    class(float64x2), intent(in) :: x
    type(float64x2) :: z
    z%limbs(1) = scale(spacing(x%limbs(1)), -53)
    z%limbs(2) = 0.0d0
  end function

  elemental function rrspacing_f(x) result(z)
    class(float64x2), intent(in) :: x
    type(float64x2) :: z
    z = abs(x) / spacing_f(x)
  end function

  elemental function neg(x) result(z)
    type(float64x2), intent(in) :: x
    type(float64x2) :: z
    z%limbs = -x%limbs
  end function

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
    if (.not. ieee_is_finite(x%limbs(1)) .or. .not. ieee_is_finite(y%limbs(1))) then
      z = (x%limbs(1) .eq. y%limbs(1))
      return
    end if
    z = (x%limbs(1) .eq. y%limbs(1)) .and. (x%limbs(2) .eq. y%limbs(2))
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
    if (.not. ieee_is_finite(x%limbs(1)) .or. .not. ieee_is_finite(d)) then
      z = (x%limbs(1) .eq. d)
      return
    end if
    z = (x%limbs(1) .eq. d) .and. (x%limbs(2) .eq. 0.0d0)
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
    if (.not. ieee_is_finite(x%limbs(1))) return
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
    if (.not. ieee_is_finite(x%limbs(1)) .or. .not. ieee_is_finite(y%limbs(1))) then
      z = (x%limbs(1) .lt. y%limbs(1))
      return
    end if
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
    if (.not. ieee_is_finite(x%limbs(1)) .or. .not. ieee_is_finite(d)) then
      z = (x%limbs(1) .lt. d)
      return
    end if
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
    if (.not. ieee_is_finite(d) .or. .not. ieee_is_finite(x%limbs(1))) then
      z = (d .lt. x%limbs(1))
      return
    end if
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
    if (.not. ieee_is_finite(x%limbs(1)) .or. .not. ieee_is_finite(y%limbs(1))) then
      z = (x%limbs(1) .le. y%limbs(1))
      return
    end if
    z = .not. lt(y, x)
  end function

  pure function le_fd(x, d) result(z)
    type(float64x2), intent(in) :: x
    double precision, intent(in) :: d
    logical :: z
    if (.not. ieee_is_finite(x%limbs(1)) .or. .not. ieee_is_finite(d)) then
      z = (x%limbs(1) .le. d)
      return
    end if
    z = .not. lt_df(d, x)
  end function

  pure function le_df(d, x) result(z)
    double precision, intent(in) :: d
    type(float64x2), intent(in) :: x
    logical :: z
    if (.not. ieee_is_finite(d) .or. .not. ieee_is_finite(x%limbs(1))) then
      z = (d .le. x%limbs(1))
      return
    end if
    z = .not. lt_fd(x, d)
  end function

  pure function ge(x, y) result(z)
    type(float64x2), intent(in) :: x, y
    logical :: z
    if (.not. ieee_is_finite(x%limbs(1)) .or. .not. ieee_is_finite(y%limbs(1))) then
      z = (x%limbs(1) .ge. y%limbs(1))
      return
    end if
    z = .not. lt(x, y)
  end function

  pure function ge_fd(x, d) result(z)
    type(float64x2), intent(in) :: x
    double precision, intent(in) :: d
    logical :: z
    if (.not. ieee_is_finite(x%limbs(1)) .or. .not. ieee_is_finite(d)) then
      z = (x%limbs(1) .ge. d)
      return
    end if
    z = .not. lt_fd(x, d)
  end function

  pure function ge_df(d, x) result(z)
    double precision, intent(in) :: d
    type(float64x2), intent(in) :: x
    logical :: z
    if (.not. ieee_is_finite(d) .or. .not. ieee_is_finite(x%limbs(1))) then
      z = (d .ge. x%limbs(1))
      return
    end if
    z = .not. lt_df(d, x)
  end function

end module
