program test_bessel_improved
    use multifloats
    use bessel_improved_mod
    implicit none
    integer, parameter :: dp = 8
    integer, parameter :: qp = 16
    type(float64x2) :: x, y
    real(qp) :: qx, qy, qgot
    real(dp) :: rel_err

    print *, "Testing Improved Bessel J0 (libquadmath sub-intervals)"
    print *, "------------------------------------------------------"

    ! Test x=1.0
    qx = 1.0_qp
    x = float64x2(real(qx, dp))
    y = bessel_j0_improved(x)
    qy = bessel_j0(qx)
    qgot = real(y%limbs(1), qp) + real(y%limbs(2), qp)
    rel_err = real(abs(qgot - qy) / abs(qy), dp)
    print *, "x = 1.0"
    print *, "  Expected J0(1) = ", qy
    print *, "  Got J0(1)      = ", qgot
    print *, "  Relative error = ", rel_err

    ! Test x=10.0
    qx = 10.0_qp
    x = float64x2(real(qx, dp))
    y = bessel_j0_improved(x)
    qy = bessel_j0(qx)
    qgot = real(y%limbs(1), qp) + real(y%limbs(2), qp)
    rel_err = real(abs(qgot - qy) / abs(qy), dp)
    print *, "x = 10.0"
    print *, "  Expected J0(10) = ", qy
    print *, "  Got J0(10)      = ", qgot
    print *, "  Relative error  = ", rel_err

    ! Test x=30.0
    qx = 30.0_qp
    x = float64x2(real(qx, dp))
    y = bessel_j0_improved(x)
    qy = bessel_j0(qx)
    qgot = real(y%limbs(1), qp) + real(y%limbs(2), qp)
    rel_err = real(abs(qgot - qy) / abs(qy), dp)
    print *, "x = 30.0"
    print *, "  Expected J0(30) = ", qy
    print *, "  Got J0(30)      = ", qgot
    print *, "  Relative error  = ", rel_err

end program test_bessel_improved
