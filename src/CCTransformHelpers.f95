subroutine CCInterpolate(samples, ninterval, order, plan_forward, plan_backward, &
                         theta_input, theta_spectrum, theta_padded, &
                         theta_fine, interpolated)
!------------------------------------------------------------------------------
! Continue an order-m meridional function onto [0, 2*pi), zero-pad its
! Fourier spectrum by two, and return nodes 0 through pi on the fine grid.
! The signed reflection follows from the parity of P_lm(cos(theta)).
! This routine is a private implementation helper for SHExpandCC and
! SHExpandCCC; it is deliberately not exposed through SHTOOLS interfaces.
!------------------------------------------------------------------------------
    use FFTW3
    use ftypes
    use, intrinsic :: iso_c_binding, only: C_PTR

    implicit none

    integer(int32), intent(in) :: ninterval, order
    complex(dp), intent(in) :: samples(ninterval+1)
    type(C_PTR), intent(in) :: plan_forward, plan_backward
    complex(dp), intent(inout) :: theta_input(2*ninterval), &
                                  theta_spectrum(2*ninterval), &
                                  theta_padded(4*ninterval), &
                                  theta_fine(4*ninterval)
    complex(dp), intent(out) :: interpolated(2*ninterval+1)
    integer(int32) :: nsource, ntheta_fine, signm

    nsource = 2*ninterval
    ntheta_fine = 2*nsource
    if (mod(order, 2) == 0) then
        signm = 1
    else
        signm = -1
    end if

    theta_input(1:ninterval+1) = samples
    if (ninterval > 1) then
        theta_input(ninterval+2:nsource) = dble(signm) * &
                                             samples(ninterval:2:-1)
    end if
    call fftw_execute_dft(plan_forward, theta_input, theta_spectrum)

    theta_padded = cmplx(0.0_dp, 0.0_dp, dp)
    theta_padded(1:ninterval) = theta_spectrum(1:ninterval)
    theta_padded(ninterval+1) = 0.5_dp*theta_spectrum(ninterval+1)
    theta_padded(ntheta_fine-ninterval+1) = 0.5_dp*theta_spectrum(ninterval+1)
    if (ninterval > 1) then
        theta_padded(ntheta_fine-ninterval+2:ntheta_fine) = &
            theta_spectrum(ninterval+2:nsource)
    end if
    call fftw_execute_dft(plan_backward, theta_padded, theta_fine)
    interpolated = theta_fine(1:nsource+1) / dble(nsource)
end subroutine CCInterpolate


subroutine CCWeights(weights, n)
!------------------------------------------------------------------------------
! Construct classical Clenshaw-Curtis weights on x_j=cos(j*pi/n), where n is
! even and equals size(weights)-1. This direct even-n form is equivalent to
! Waldvogel (2006) and integrates constants to two. This private helper is
! used only by the CC analysis routines.
!------------------------------------------------------------------------------
    use ftypes

    implicit none

    integer(int32), intent(in) :: n
    real(dp), intent(out) :: weights(n+1)
    integer(int32) :: j, k
    real(dp) :: pi, theta, value

    pi = acos(-1.0_dp)
    weights(1) = 1.0_dp / dble(n*n-1)
    weights(n+1) = weights(1)
    do j = 1, n-1
        theta = pi*dble(j)/dble(n)
        value = 1.0_dp
        do k = 1, n/2-1
            value = value - 2.0_dp*cos(2*k*theta) / dble(4*k*k-1)
        end do
        value = value - cos(n*theta) / dble(n*n-1)
        weights(j+1) = 2.0_dp*value / dble(n)
    end do
end subroutine CCWeights


real(dp) function CCPlmIntegral(l, m, norm, cnorm) result(integral)
!------------------------------------------------------------------------------
! Return integral_{-1}^{1} P_lm(x)^2 dx for the SHTOOLS normalization used by
! the caller. CNORM distinguishes the complex and real associated-Legendre
! conventions of PlmBar, PlmSchmidt, and PlmON. This private helper is used
! only by the CC analysis routines.
!------------------------------------------------------------------------------
    use ftypes

    implicit none

    integer(int32), intent(in) :: l, m, norm, cnorm
    integer(int32) :: k
    real(dp) :: pi, factor

    integral = 0.0_dp
    select case (norm)
        case (1)
            if (cnorm == 1 .or. m == 0) then
                integral = 2.0_dp
            else
                integral = 4.0_dp
            end if
        case (2)
            if (cnorm == 1 .or. m == 0) then
                integral = 2.0_dp / dble(2*l+1)
            else
                integral = 4.0_dp / dble(2*l+1)
            end if
        case (3)
            factor = 1.0_dp
            do k = l-m+1, l+m
                factor = factor*dble(k)
            end do
            integral = 2.0_dp*factor / dble(2*l+1)
        case (4)
            pi = acos(-1.0_dp)
            if (cnorm == 1 .or. m == 0) then
                integral = 1.0_dp / (2.0_dp*pi)
            else
                integral = 1.0_dp / pi
            end if
    end select
end function CCPlmIntegral
