program consumer
  ! Smoke test for the source-distributed multifloatsf package. Exercises:
  !   - the type(real64x2) derived type (layout, %limbs accessor),
  !   - an inline elemental arithmetic operator (`+` on DD), and
  !   - a transcendental that routes through the C ABI (sqrt → sqrtdd).
  ! The printed values are stable to the final ULP so the driver script
  ! can diff them verbatim.
  use multifloats
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none

  type(real64x2) :: x, y, z, r

  ! x = 1.0 exactly (zero lo limb).
  x = 1.0_dp
  ! y = 2.0 exactly.
  y = 2.0_dp
  ! Inline elemental DD add: z = 3.0 exact.
  z = x + y
  ! Transcendental via bind(c): r = sqrt(2.0).
  r = sqrt(y)

  print '(a,2es24.17)', "z      = ", z%limbs(1), z%limbs(2)
  print '(a,2es24.17)', "sqrt2  = ", r%limbs(1), r%limbs(2)
  print '(a)', "ok"
end program consumer
