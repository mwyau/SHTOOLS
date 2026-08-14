subroutine SHExpandCC(gridcc, cilm, lmax, norm, csphase, lmax_calc, &
                      exitstatus)
!------------------------------------------------------------------------------
!
!   Expand a real, endpoint-including Clenshaw-Curtis grid into real
!   spherical-harmonic coefficients. The input grid has L+2 colatitudes
!   theta=j*pi/(L+1), including both poles, and 2*L+2 longitudes.
!
!   The direct regular-grid projection is evaluated as follows:
!       1. Transform every latitude ring with a real-to-complex FFT.
!       2. Extract each longitudinal order m through LMAX_CALC.
!       3. Continue its meridional samples with the parity
!          f_m(2*pi-theta)=(-1)^m f_m(theta).
!       4. Zero-pad the meridional Fourier spectrum by two and inverse
!          transform it onto a twice-finer colatitude grid.
!       5. Form classical Clenshaw-Curtis weights on the 2*L+3 fine nodes.
!       6. Project the interpolated orders onto associated Legendre functions.
!
!   The continuation is valid because P_lm(cos(theta)) is sin(theta)^m times
!   a polynomial in cos(theta). The fine quadrature integrates the degree-2L
!   products in the projection. The algorithm does not construct a dense
!   latitude-by-degree matrix and does not use a least-squares solve.
!
!   Calling Parameters
!
!       IN
!           gridcc      Real, non-extended grid of dimensions (L+2, 2*L+2).
!
!       OUT
!           cilm        Real cosine and sine coefficients, dimensioned at
!                       least (2, LMAX_CALC+1, LMAX_CALC+1).
!           lmax        Maximum recoverable degree inferred as L=nlat-2.
!
!       OPTIONAL (IN)
!           norm        Legendre normalization: 1 = 4-pi normalized
!                       (default), 2 = Schmidt, 3 = unnormalized, and
!                       4 = orthonormal.
!           csphase     1 excludes the Condon-Shortley phase (default);
!                       -1 includes it.
!           lmax_calc   Maximum degree recovered; default is LMAX.
!
!       OPTIONAL (OUT)
!           exitstatus  0 = no error, 1 = invalid dimensions, 2 = invalid
!                       input, 3 = allocation or FFTW-plan failure.
!
!   Notes:
!       CC has one layout and no Driscoll and Healy sampling parameter.
!       Unnormalized Legendre functions are accurate only to approximately
!       degree 15.
!
!   References:
!       Clenshaw, C. W., and A. R. Curtis (1960), Numerische Mathematik, 2,
!       197-205, https://doi.org/10.1007/BF01386223.
!       Waldvogel, J. (2006), BIT Numerical Mathematics, 46, 195-202,
!       https://doi.org/10.1007/s10543-006-0045-4.
!       Holmes, S. A., and W. E. Featherstone (2002), Journal of Geodesy, 76,
!       279-299, https://doi.org/10.1007/s00190-002-0216-2.
!
!   Copyright (c) 2026, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use FFTW3
    use SHTOOLS, only: PlmBar, PlmSchmidt, PLegendreA, PlmON
    use ftypes
    use, intrinsic :: iso_c_binding, only: C_PTR, C_NULL_PTR, C_ASSOCIATED

    implicit none

    real(dp), intent(in) :: gridcc(:,:)
    real(dp), intent(out) :: cilm(:,:,:)
    integer(int32), intent(out) :: lmax
    integer(int32), intent(in), optional :: norm, csphase, lmax_calc
    integer(int32), intent(out), optional :: exitstatus
    integer(int32) :: nlat, nlong, lmax_comp, lnorm, phase, np, &
                      ntheta_coarse, ntheta_quadrature, ntheta_fine, l, m, i, &
                      idx, istat
    real(dp) :: pi, z, integral
    real(dp), allocatable :: p(:), weights(:), row(:)
    complex(dp), allocatable :: spectrum(:), longitude(:,:), meridional(:,:), &
                              theta_input(:), theta_spectrum(:), &
                              theta_padded(:), theta_fine(:)
    type(C_PTR) :: plan_long, plan_theta_forward, plan_theta_backward
    real(dp) :: CCPlmIntegral
    external :: CCInterpolate, CCWeights, CCPlmIntegral

    if (present(exitstatus)) exitstatus = 0
    plan_long = C_NULL_PTR
    plan_theta_forward = C_NULL_PTR
    plan_theta_backward = C_NULL_PTR

    nlat = size(gridcc, 1)
    nlong = size(gridcc, 2)
    if (nlat < 2 .or. nlong /= 2*(nlat-1)) then
        print*, 'Error --- SHExpandCC'
        print*, 'GRIDCC must have shape (L+2, 2*L+2).'
        print*, 'Input dimensions are ', nlat, nlong
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if
    end if

    lmax = nlat - 2
    if (present(lmax_calc)) then
        lmax_comp = lmax_calc
    else
        lmax_comp = lmax
    end if
    if (lmax_comp < 0 .or. lmax_comp > lmax) then
        print*, 'Error --- SHExpandCC'
        print*, 'LMAX_CALC must be between 0 and LMAX.'
        print*, 'LMAX = ', lmax
        print*, 'LMAX_CALC = ', lmax_comp
        if (present(exitstatus)) then
            exitstatus = 2
            return
        else
            stop
        end if
    end if

    if (size(cilm, 1) < 2 .or. size(cilm, 2) < lmax_comp+1 .or. &
        size(cilm, 3) < lmax_comp+1) then
        print*, 'Error --- SHExpandCC'
        print*, 'CILM must be dimensioned as (2, LMAX_CALC+1, LMAX_CALC+1).'
        print*, 'Input dimensions are ', size(cilm, 1), size(cilm, 2), &
                size(cilm, 3)
        if (present(exitstatus)) then
            exitstatus = 1
            return
        else
            stop
        end if
    end if

    if (present(norm)) then
        lnorm = norm
    else
        lnorm = 1
    end if
    if (lnorm < 1 .or. lnorm > 4) then
        print*, 'Error --- SHExpandCC'
        print*, 'NORM must be between 1 and 4.'
        print*, 'Input value is ', lnorm
        if (present(exitstatus)) then
            exitstatus = 2
            return
        else
            stop
        end if
    end if

    if (present(csphase)) then
        phase = csphase
    else
        phase = 1
    end if
    if (phase /= -1 .and. phase /= 1) then
        print*, 'Error --- SHExpandCC'
        print*, 'CSPHASE must be -1 or 1.'
        print*, 'Input value is ', phase
        if (present(exitstatus)) then
            exitstatus = 2
            return
        else
            stop
        end if
    end if

    np = (lmax_comp+1)*(lmax_comp+2)/2
    ntheta_coarse = 2*(lmax+1)
    ntheta_quadrature = ntheta_coarse
    ntheta_fine = 2*ntheta_coarse
    allocate(p(np), weights(ntheta_quadrature+1), row(nlong), &
             spectrum(nlong/2+1), longitude(nlat, lmax_comp+1), &
             meridional(ntheta_quadrature+1, lmax_comp+1), &
             theta_input(ntheta_coarse), theta_spectrum(ntheta_coarse), &
             theta_padded(ntheta_fine), theta_fine(ntheta_fine), stat=istat)
    if (istat /= 0) then
        print*, 'Error --- SHExpandCC'
        print*, 'Problem allocating CC transform work arrays.'
        if (present(exitstatus)) then
            exitstatus = 3
            return
        else
            stop
        end if
    end if

    plan_long = fftw_plan_dft_r2c_1d(nlong, row, spectrum, FFTW_MEASURE)
    plan_theta_forward = fftw_plan_dft_1d(ntheta_coarse, theta_input, &
                                           theta_spectrum, FFTW_FORWARD, &
                                           FFTW_MEASURE)
    plan_theta_backward = fftw_plan_dft_1d(ntheta_fine, theta_padded, &
                                            theta_fine, FFTW_BACKWARD, &
                                            FFTW_MEASURE)
    if (.not. C_ASSOCIATED(plan_long) .or. &
        .not. C_ASSOCIATED(plan_theta_forward) .or. &
        .not. C_ASSOCIATED(plan_theta_backward)) then
        if (C_ASSOCIATED(plan_long)) call fftw_destroy_plan(plan_long)
        if (C_ASSOCIATED(plan_theta_forward)) then
            call fftw_destroy_plan(plan_theta_forward)
        end if
        if (C_ASSOCIATED(plan_theta_backward)) then
            call fftw_destroy_plan(plan_theta_backward)
        end if
        deallocate(p, weights, row, spectrum, longitude, meridional, &
                   theta_input, theta_spectrum, theta_padded, theta_fine)
        print*, 'Error --- SHExpandCC'
        print*, 'Problem creating FFTW plans.'
        if (present(exitstatus)) then
            exitstatus = 3
            return
        else
            stop
        end if
    end if

    do i = 1, nlat
        row = gridcc(i, 1:nlong)
        call fftw_execute_dft_r2c(plan_long, row, spectrum)
        longitude(i, 1:lmax_comp+1) = spectrum(1:lmax_comp+1) / dble(nlong)
    end do
    call fftw_destroy_plan(plan_long)
    plan_long = C_NULL_PTR

    do m = 0, lmax_comp
        call CCInterpolate(longitude(:, m+1), ntheta_coarse/2, m, &
                           plan_theta_forward, plan_theta_backward, &
                           theta_input, theta_spectrum, theta_padded, &
                           theta_fine, meridional(:, m+1))
    end do
    call fftw_destroy_plan(plan_theta_forward)
    call fftw_destroy_plan(plan_theta_backward)
    plan_theta_forward = C_NULL_PTR
    plan_theta_backward = C_NULL_PTR

    call CCWeights(weights, ntheta_quadrature)
    cilm = 0.0_dp
    pi = acos(-1.0_dp)
    do i = 1, ntheta_quadrature+1
        z = cos(pi*dble(i-1)/dble(ntheta_quadrature))
        select case (lnorm)
            case (1)
                call PlmBar(p, lmax_comp, z, csphase=phase, cnorm=0, &
                            exitstatus=istat)
            case (2)
                call PlmSchmidt(p, lmax_comp, z, csphase=phase, cnorm=0, &
                                exitstatus=istat)
            case (3)
                call PLegendreA(p, lmax_comp, z, csphase=phase, &
                                exitstatus=istat)
            case (4)
                call PlmON(p, lmax_comp, z, csphase=phase, cnorm=0, &
                           exitstatus=istat)
        end select
        if (istat /= 0) then
            deallocate(p, weights, row, spectrum, longitude, meridional, &
                       theta_input, theta_spectrum, theta_padded, theta_fine)
            print*, 'Error --- SHExpandCC'
            print*, 'Problem evaluating associated Legendre functions.'
            if (present(exitstatus)) then
                exitstatus = istat
                return
            else
                stop
            end if
        end if

        do m = 0, lmax_comp
            do l = m, lmax_comp
                idx = l*(l+1)/2 + m + 1
                if (m == 0) then
                    cilm(1, l+1, 1) = cilm(1, l+1, 1) + weights(i) * &
                                       real(meridional(i, 1), dp) * p(idx)
                else
                    cilm(1, l+1, m+1) = cilm(1, l+1, m+1) + &
                                          2.0_dp*weights(i) * &
                                          real(meridional(i, m+1), dp) * p(idx)
                    cilm(2, l+1, m+1) = cilm(2, l+1, m+1) - &
                                          2.0_dp*weights(i) * &
                                          aimag(meridional(i, m+1)) * p(idx)
                end if
            end do
        end do
    end do

    do m = 0, lmax_comp
        do l = m, lmax_comp
            integral = CCPlmIntegral(l, m, lnorm, 0)
            cilm(1, l+1, m+1) = cilm(1, l+1, m+1) / integral
            if (m > 0) then
                cilm(2, l+1, m+1) = cilm(2, l+1, m+1) / integral
            end if
        end do
    end do

    deallocate(p, weights, row, spectrum, longitude, meridional, theta_input, &
               theta_spectrum, theta_padded, theta_fine)
end subroutine SHExpandCC
