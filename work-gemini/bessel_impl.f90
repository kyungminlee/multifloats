module bessel_improved_mod
    use multifloats
    use, intrinsic :: ieee_arithmetic
    implicit none
    private

    public :: bessel_j0_improved

    contains

    include 'bessel_constants.f90.inc'

    pure elemental function bessel_j0_improved(x) result(res)
        type(float64x2), intent(in) :: x
        type(float64x2) :: res, z, p, q, s, c, ss, cc, z_trig, tmp, x_arg, isqrt2
        real(dp) :: xx, xinv
        if (.not. ieee_is_finite(x%limbs(1))) then
            res%limbs = [bessel_j0(x%limbs(1)), 0.0_dp]
            return
        end if
        xx = abs(x%limbs(1))
        if (xx == 0.0_dp) then
            res = MF_ONE
            return
        end if

        if (xx <= 2.0_dp) then
            z = x * x
            res = z * z * dd_neval(z, j0_J0_2N_hi, j0_J0_2N_lo) / &
                          dd_deval(z, j0_J0_2D_hi, j0_J0_2D_lo)
            res = res - z * 0.25_dp + MF_ONE
            return
        end if

        x_arg%limbs = [xx, 0.0_dp]
        s = dd_sin_full(x_arg)
        c = dd_cos_full(x_arg)
        isqrt2%limbs = [INV_SQRT2_hi, INV_SQRT2_lo]
        ss = (s - c) * isqrt2
        cc = (s + c) * isqrt2
        x_arg%limbs = [xx + xx, 0.0_dp]
        z_trig = -dd_cos_full(x_arg)
        if (s%limbs(1) * c%limbs(1) < 0.0_dp) then
            cc = (z_trig / (s - c)) * isqrt2
        else
            ss = (z_trig / (s + c)) * isqrt2
        end if

        xinv = 1.0_dp / xx
        z%limbs = [xinv * xinv, 0.0_dp]
        if (xinv <= 0.25_dp) then
            if (xinv <= 0.125_dp) then
                if (xinv <= 0.0625_dp) then
                    p = dd_neval(z, j0_P16_IN_hi, j0_P16_IN_lo) / &
                        dd_deval(z, j0_P16_ID_hi, j0_P16_ID_lo)
                    q = dd_neval(z, j0_Q16_IN_hi, j0_Q16_IN_lo) / &
                        dd_deval(z, j0_Q16_ID_hi, j0_Q16_ID_lo)
                else
                    p = dd_neval(z, j0_P8_16N_hi, j0_P8_16N_lo) / &
                        dd_deval(z, j0_P8_16D_hi, j0_P8_16D_lo)
                    q = dd_neval(z, j0_Q8_16N_hi, j0_Q8_16N_lo) / &
                        dd_deval(z, j0_Q8_16D_hi, j0_Q8_16D_lo)
                end if
            else if (xinv <= 0.1875_dp) then
                p = dd_neval(z, j0_P5_8N_hi, j0_P5_8N_lo) / &
                    dd_deval(z, j0_P5_8D_hi, j0_P5_8D_lo)
                q = dd_neval(z, j0_Q5_8N_hi, j0_Q5_8N_lo) / &
                    dd_deval(z, j0_Q5_8D_hi, j0_Q5_8D_lo)
            else
                p = dd_neval(z, j0_P4_5N_hi, j0_P4_5N_lo) / &
                    dd_deval(z, j0_P4_5D_hi, j0_P4_5D_lo)
                q = dd_neval(z, j0_Q4_5N_hi, j0_Q4_5N_lo) / &
                    dd_deval(z, j0_Q4_5D_hi, j0_Q4_5D_lo)
            end if
        else
            if (xinv <= 0.375_dp) then
                if (xinv <= 0.3125_dp) then
                    p = dd_neval(z, j0_P3r2_4N_hi, j0_P3r2_4N_lo) / &
                        dd_deval(z, j0_P3r2_4D_hi, j0_P3r2_4D_lo)
                    q = dd_neval(z, j0_Q3r2_4N_hi, j0_Q3r2_4N_lo) / &
                        dd_deval(z, j0_Q3r2_4D_hi, j0_Q3r2_4D_lo)
                else
                    p = dd_neval(z, j0_P2r7_3r2N_hi, j0_P2r7_3r2N_lo) / &
                        dd_deval(z, j0_P2r7_3r2D_hi, j0_P2r7_3r2D_lo)
                    q = dd_neval(z, j0_Q2r7_3r2N_hi, j0_Q2r7_3r2N_lo) / &
                        dd_deval(z, j0_Q2r7_3r2D_hi, j0_Q2r7_3r2D_lo)
                end if
            else if (xinv <= 0.4375_dp) then
                p = dd_neval(z, j0_P2r3_2r7N_hi, j0_P2r3_2r7N_lo) / &
                    dd_deval(z, j0_P2r3_2r7D_hi, j0_P2r3_2r7D_lo)
                q = dd_neval(z, j0_Q2r3_2r7N_hi, j0_Q2r3_2r7N_lo) / &
                    dd_deval(z, j0_Q2r3_2r7D_hi, j0_Q2r3_2r7D_lo)
            else
                p = dd_neval(z, j0_P2_2r3N_hi, j0_P2_2r3N_lo) / &
                    dd_deval(z, j0_P2_2r3D_hi, j0_P2_2r3D_lo)
                q = dd_neval(z, j0_Q2_2r3N_hi, j0_Q2_2r3N_lo) / &
                    dd_deval(z, j0_Q2_2r3D_hi, j0_Q2_2r3D_lo)
            end if
        end if
        p = MF_ONE + z * p
        q = (z * q - 0.125_dp) * xinv
        tmp%limbs = [TWO_OVER_PI_hi, TWO_OVER_PI_lo]
        x_arg%limbs = [xx, 0.0_dp]
        res = sqrt(tmp / x_arg) * (p * cc - q * ss)
    end function

end module bessel_improved_mod
