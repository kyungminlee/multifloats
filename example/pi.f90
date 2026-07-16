! Compute π to full double-double precision as 4·atan(1) and print it.
!
! Build against an installed multifloatsf package (see README.md).
program pi_example
  use multifloats                              ! generic atan + defined I/O
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none

  type(real64x2) :: x, pi

  x  = real64x2(1.0_dp)
  pi = atan(x) * 4                             ! π to ~106 bits
  print '(a,2es24.17)', "pi limbs = ", pi%limbs(1), pi%limbs(2)
  print *, pi                                  ! defined I/O: ~32 digits
end program pi_example
