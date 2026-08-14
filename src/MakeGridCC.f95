subroutine MakeGridCC(gridcc, cilm, lmax, norm, csphase, lmax_calc, extend, &
                      exitstatus)
!------------------------------------------------------------------------------
!
!   Evaluate real spherical-harmonic coefficients on the endpoint-including
!   Clenshaw-Curtis grid. The grid has LMAX+2 colatitudes
!   theta=j*pi/(LMAX+1), including both poles, and 2*LMAX+2 longitudes.
!   Longitude synthesis uses an FFT for each latitude ring. Associated
!   Legendre functions are calculated with the SHTOOLS recurrences described
!   by Holmes and Featherstone (2002).
!
!   Calling Parameters
!
!       IN
!           cilm        Real cosine and sine coefficients, dimensioned at
!                       least (2, LMAX+1, LMAX+1).
!           lmax        Maximum degree defining the output grid.
!
!       OUT
!           gridcc      Real grid, dimensioned at least (LMAX+2,
!                       2*LMAX+2) when EXTEND=0 and (LMAX+2,
!                       2*LMAX+3) when EXTEND=1.
!
!       OPTIONAL (IN)
!           norm        Legendre normalization: 1 = 4-pi normalized
!                       (default), 2 = Schmidt, 3 = unnormalized, and
!                       4 = orthonormal.
!           csphase     1 excludes the Condon-Shortley phase (default);
!                       -1 includes it.
!           lmax_calc   Maximum degree synthesized; default is LMAX.
!           extend      1 appends the 360-degree longitude column; default
!                       is 0.
!
!       OPTIONAL (OUT)
!           exitstatus  0 = no error, 1 = invalid dimensions, 2 = invalid
!                       input, 3 = allocation or FFTW-plan failure.
!
!   Notes:
!       The output has one CC layout and no Driscoll and Healy sampling
!       parameter. Unnormalized Legendre functions are accurate only to
!       approximately degree 15.
!
!   References:
!       Holmes, S. A., and W. E. Featherstone (2002), Journal of Geodesy,
!       76, 279-299, https://doi.org/10.1007/s00190-002-0216-2.
!
!   Copyright (c) 2026, SHTOOLS
!   All rights reserved.
!
!------------------------------------------------------------------------------
    use FFTW3
    use SHTOOLS, only: PlmBar, PlmSchmidt, PLegendreA, PlmON
    use ftypes
    use, intrinsic :: iso_c_binding, only: C_PTR, C_ASSOCIATED

    implicit none

    real(dp), intent(out) :: gridcc(:,:)
    real(dp), intent(in) :: cilm(:,:,:)
    integer(int32), intent(in) :: lmax
    integer(int32), intent(in), optional :: norm, csphase, lmax_calc, extend
    integer(int32), intent(out), optional :: exitstatus
    integer(int32) :: l, m, i, lmax_comp, lcoeff, lnorm, phase, extend_grid, &
                      nlat, nlong, nlong_out, np, idx, istat
    real(dp) :: pi, theta, z
    real(dp), allocatable :: p(:), values(:)
    complex(dp), allocatable :: coef(:)
    type(C_PTR) :: plan

    if (present(exitstatus)) exitstatus = 0

    if (lmax < 0 .or. size(cilm(:,1,1)) < 2) then
        if (present(exitstatus)) then
            exitstatus = 1
            return
        end if
        print*, 'Error --- MakeGridCC'
        print*, 'LMAX must be nonnegative and CILM must have two components.'
        stop
    end if

    lnorm = 1
    if (present(norm)) lnorm = norm
    if (lnorm < 1 .or. lnorm > 4) then
        if (present(exitstatus)) then
            exitstatus = 2
            return
        end if
        print*, 'Error --- MakeGridCC'
        print*, 'NORM must be between 1 and 4.'
        stop
    end if

    phase = 1
    if (present(csphase)) phase = csphase
    if (phase /= -1 .and. phase /= 1) then
        if (present(exitstatus)) then
            exitstatus = 2
            return
        end if
        print*, 'Error --- MakeGridCC'
        print*, 'CSPHASE must be -1 or 1.'
        stop
    end if

    lmax_comp = lmax
    if (present(lmax_calc)) lmax_comp = lmax_calc
    if (lmax_comp < 0 .or. lmax_comp > lmax) then
        if (present(exitstatus)) then
            exitstatus = 2
            return
        end if
        print*, 'Error --- MakeGridCC'
        print*, 'LMAX_CALC must be between 0 and LMAX.'
        stop
    end if

    extend_grid = 0
    if (present(extend)) extend_grid = extend
    if (extend_grid /= 0 .and. extend_grid /= 1) then
        if (present(exitstatus)) then
            exitstatus = 2
            return
        end if
        print*, 'Error --- MakeGridCC'
        print*, 'EXTEND must be 0 or 1.'
        stop
    end if

    nlat = lmax + 2
    nlong = 2*lmax + 2
    nlong_out = nlong + extend_grid
    if (size(gridcc(:,1)) < nlat .or. size(gridcc(1,:)) < nlong_out) then
        if (present(exitstatus)) then
            exitstatus = 1
            return
        end if
        print*, 'Error --- MakeGridCC'
        print*, 'GRIDCC has incompatible dimensions.'
        stop
    end if

    np = (lmax_comp + 1)*(lmax_comp + 2)/2
    lcoeff = min(lmax_comp, size(cilm(1,:,1))-1, size(cilm(1,1,:))-1)
    allocate(p(np), values(nlong), coef(nlong/2+1), stat=istat)
    if (istat /= 0) then
        if (present(exitstatus)) then
            exitstatus = 3
            return
        end if
        print*, 'Error --- MakeGridCC'
        print*, 'Problem allocating CC synthesis work arrays.'
        stop
    end if

    plan = fftw_plan_dft_c2r_1d(nlong, coef, values, FFTW_MEASURE)
    if (.not. C_ASSOCIATED(plan)) then
        deallocate(p, values, coef)
        print*, 'Error --- MakeGridCC'
        print*, 'Problem creating the FFTW plan.'
        if (present(exitstatus)) then
            exitstatus = 3
            return
        else
            stop
        end if
    end if
    pi = acos(-1.0_dp)
    do i = 1, nlat
        theta = pi*dble(i-1)/dble(nlat-1)
        z = cos(theta)
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
            call fftw_destroy_plan(plan)
            if (present(exitstatus)) then
                exitstatus = istat
                return
            end if
            print*, 'Error --- MakeGridCC'
            print*, 'Problem evaluating associated Legendre functions.'
            stop
        end if

        coef = cmplx(0.0_dp, 0.0_dp, dp)
        do l = 0, lcoeff
            idx = l*(l+1)/2 + 1
            coef(1) = coef(1) + cilm(1,l+1,1)*p(idx)
        end do
        do m = 1, lcoeff
            do l = m, lcoeff
                idx = l*(l+1)/2 + m + 1
                coef(m+1) = coef(m+1) + 0.5_dp*p(idx)* &
                    cmplx(cilm(1,l+1,m+1), -cilm(2,l+1,m+1), dp)
            end do
        end do
        call fftw_execute_dft_c2r(plan, coef, values)
        gridcc(i,1:nlong) = values
    end do
    if (extend_grid == 1) gridcc(1:nlat,nlong_out) = gridcc(1:nlat,1)
    call fftw_destroy_plan(plan)
    deallocate(p, values, coef)
end subroutine MakeGridCC
