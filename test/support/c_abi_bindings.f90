! Shared bind(c) mirror of the scalar-DD C ABI exported by
! src/multifloats_math.cc (declared in include/multifloats/float64x2.h):
! the dd_c layout twin of real64x2 plus interfaces for the kernels the
! ABI tests drive directly. Used by test/unit/bench_abi.f90 (overhead
! baseline) and test/integration/abi_equivalence.f90 (native vs C-ABI
! fuzz equivalence); test/integration/crosscheck_bindings.f90 reuses
! the dd_c type.
module c_abi_bindings
  use, intrinsic :: iso_c_binding, only: c_double, c_int
  implicit none

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
    pure type(dd_c) function c_powidd(a, n) bind(c, name='powidd')
      import :: dd_c, c_int; type(dd_c), value :: a; integer(c_int), value :: n
    end function
    pure type(dd_c) function c_dd_modulo(a, b) bind(c, name='modulodd')
      import :: dd_c; type(dd_c), value :: a, b
    end function
  end interface

end module c_abi_bindings
