!> Solve the layer continuity equation using the WENO method for layer fluxes.
module MOM_continuity_WENO

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_cpu_clock, only : cpu_clock_id, cpu_clock_begin, cpu_clock_end, CLOCK_ROUTINE
use MOM_diag_mediator, only : time_type, diag_ctrl
use MOM_error_handler, only : MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser, only : get_param, log_version, param_file_type
use MOM_grid, only : ocean_grid_type
use MOM_open_boundary, only : ocean_OBC_type, OBC_segment_type, OBC_NONE
use MOM_open_boundary, only : OBC_DIRECTION_E, OBC_DIRECTION_W, OBC_DIRECTION_N, OBC_DIRECTION_S
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : BT_cont_type, porous_barrier_type
use MOM_verticalGrid, only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public weno3_reconstruction_interface
public weno5_reconstruction_interface
public weno7_reconstruction_interface
public WENO_limiter

contains

!> 3th-order weno z-type reconstruction flux
subroutine weno3_reconstruction_interface(wmR, wpL, q, h_min)
  real, intent(in) :: q(3)       !< tracer concentration for 3-stencil wide
  real, intent(in)  :: h_min     !< The minimum thickness
  real, intent(out) :: wmR, wpL  !< weno reconstruction at the cell interfaces i-1/2 and i+1/2

  real :: P1, P2         ! reconstructed polynomials
  real :: b1, b2         ! smoothness indicator
  real :: w1, w2         ! nonlinear weights
  real :: d1, d2         ! linear weights
  real :: eps,  wnorm, tau
  integer, parameter :: r = 2

  ! linear weights
  d1 = 1.0/3.0 ; d2 = 2.0/3.0
  eps = 1.0e-20

  ! Compute flux at the right side of i+1/2
  ! reconstructed polynomials
  P1 = 0.5*(-q(1) + 3.0*q(2))
  P2 = 0.5*(q(2) + q(3))

  ! smoothness indicator
  b1 = (q(2)-q(1))*(q(2)-q(1))
  b2 = (q(3)-q(2))*(q(3)-q(2))

  ! Alpha values
  tau = abs(b1-b2)
  w1 = d1*fac_fn(tau, b1)
  w2 = d2*fac_fn(tau, b2)

  ! Normalization
  wnorm = w1+w2
  wpL = (w1*P1 + w2*P2)/wnorm

  ! Compute flux at the right side of i-1/2
  ! reconstructed polynomials
  P1 = 0.5*(-q(3) + 3.0*q(2))
  P2 = 0.5*(q(2) + q(1))

  ! smoothness indicator
  b1 = (q(2)-q(3))*(q(2)-q(3))
  b2 = (q(1)-q(2))*(q(1)-q(2))

  ! Alpha values
  tau = abs(b1-b2)
  w1 = d1*fac_fn(tau, b1)
  w2 = d2*fac_fn(tau, b2)

  ! Normalization
  wnorm = w1+w2
  wmR = (w1*P1 + w2*P2)/wnorm

  !call PP_limiter(q0, wmR, wpL, h_min)

end subroutine weno3_reconstruction_interface

!> 5th-order weno z-type reconstruction flux
subroutine weno5_reconstruction_interface(wmR, wpL, q, h_min)
  real, intent(in) :: q(5)       !< tracer concentration for 5-stencil wide
  real, intent(in)  :: h_min     !< The minimum thickness
  real, intent(out) :: wmR, wpL  !< weno reconstruction at the cell interfaces i-1/2 and i+1/2

  real :: P0, P1, P2         ! reconstructed polynomials
  real :: b0, b1, b2         ! smoothness indicator
  real :: w0, w1, w2         ! nonlinear weights
  real :: d0, d1, d2         ! linear weights
  real :: a0, a1, a2
  real :: eps,  wnorm, tau
  integer, parameter :: r = 2

  ! linear weights
  d0 = 1.0/10.0 ; d1 = 6.0/10.0 ; d2 = 3.0/10.0
  eps = 1.0e-20

  ! Compute flux at left side of i+1/2
  ! First stencil
  P0 = (2.0*q(1) - 7.0*q(2) + 11.0*q(3))/6.0
  b0 = (13.0/12.0)*(q(1) - 2.0*q(2) + q(3))**2 + 0.25*(q(1) - 4.0*q(2) + 3.0*q(3))**2

  ! Second stencil
  P1 = (-q(2) + 5.0*q(3) + 2.0*q(4))/6.0
  b1 = (13.0/12.0)*(q(2) - 2.0*q(3) + q(4))**2 + 0.25*(q(2) - q(4))**2

  ! Third stencil
  P2 = (2.0*q(3) + 5.0*q(4) - q(5))/6.0
  b2 = (13.0/12.0)*(q(3) - 2.0*q(4) + q(5))**2 + 0.25*(3.0*q(3) - 4.0*q(4) + q(5))**2

  ! Alpha values
  tau = abs(b2-b0)
  w0 = d0*fac_fn(tau, b0)
  w1 = d1*fac_fn(tau, b1)
  w2 = d2*fac_fn(tau, b2)

  wnorm = w0+w1+w2
  wpL = (w0*P0 + w1*P1 + w2*P2)/wnorm

  ! Compute flux at the right side of i-1/2
  !d0 = 3.0/10.0 ; d1 = 6.0/10.0 ; d2 = 1.0/10.0
  ! First stencil
  P0 = (2.0*q(5) - 7.0*q(4) + 11.0*q(3))/6.0
  b0 = (13.0/12.0)*(q(5) - 2.0*q(4) + q(3))**2 + 0.25*(q(5) - 4.0*q(4) + 3.0*q(3))**2

  ! Second stencil
  P1 = (-q(4) + 5.0*q(3) + 2.0*q(2))/6.0
  b1 = (13.0/12.0)*(q(4) - 2.0*q(3) + q(2))**2 + 0.25*(q(4) - q(3))**2

  ! Third stencil
  P2 = (2.0*q(3) + 5.0*q(2) - q(1))/6.0
  b2 = (13.0/12.0)*(q(3) - 2.0*q(2) + q(1))**2 + 0.25*(3.0*q(3) - 4.0*q(2) + q(1))**2

  ! Alpha values
  tau = abs(b2-b0)
  w0 = d0*fac_fn(tau, b0)
  w1 = d1*fac_fn(tau, b1)
  w2 = d2*fac_fn(tau, b2)

  wnorm = w0+w1+w2
  wmR = (w0*P0 + w1*P1 + w2*P2)/wnorm

  !call PP_limiter(q0, wmR, wpL, h_min)

end subroutine weno5_reconstruction_interface

!> 7th-order weno z-type reconstruction flux
subroutine weno7_reconstruction_interface(wmR, wpL, q, h_min)
  real, intent(in) :: q(7) !< tracer concentration for 7-stencil wide
  real, intent(in)  :: h_min     !< The minimum thickness
  real, intent(out) :: wmR, wpL  !< weno reconstruction at the cell interfaces i-1/2 and i+1/2

  real :: P0, P1, P2, P3     ! reconstructed polynomials
  real :: b0, b1, b2, b3     ! smoothness indicator
  real :: w0, w1, w2, w3     ! nonlinear weights
  real :: d0, d1, d2, d3     ! nonlinear weights
  real :: a0, a1, a2, a3
  real :: eps, tau, wnorm
  integer, parameter :: r = 2

  ! linear weights
  d0 = 1.0/35.0 ;  d1 = 12.0/35.0 ; d2 = 18.0/35.0 ; d3 = 4.0/35.0
  eps = 1.0e-20

  ! Compute flux at the right side of i+1/2
  ! 1st stencil
  b0 = q(1)*(547.0*q(1) - 3882.0*q(2) + 4642.0*q(3) - 1854.0*q(4)) + &
        q(2)*(7043.0*q(2) - 17246.0*q(3) + 7042.0*q(4)) + &
        q(3)*(11003.0*q(3) - 9402.0*q(4)) + 2107.0*q(4)**2
  P0 = (-3.0*q(1) + 13.0*q(2) - 23.0*q(3) + 25.0*q(4))/12.0

  ! 2nd stencil
  b1 = q(2)*(267.0*q(2) - 1642.0*q(3) + 1602.0*q(4) - 494.0*q(5)) + &
         q(3)*(2843.0*q(3) - 5966.0*q(4) + 1922.0*q(5)) &
         + q(4)*(3443.0*q(4) - 2522.0*q(5)) + 547.0*q(5)**2
  P1 = (q(2) - 5.0*q(3) + 13.0*q(4) + 3.0*q(5))/12.0

  ! 3rd stencil
  b2 = q(3)*(547.0*q(3) - 2522.0*q(4) + 1922.0*q(5) - 494.0*q(6)) + &
         q(4)*(3443.0*q(4) - 5966.0*q(5) + 1602.0*q(6)) &
         + q(5)*(2843.0*q(5) - 1642.0*q(6)) + 267.0*q(6)**2
  P2 = (-q(3) + 7.0*q(4) + 7.0*q(5) - q(6))/12.0

  ! 4rd stencil
  b3 = q(4)*(2107.0*q(4) - 9402.0*q(5) + 7042.0*q(6) - 1854.0*q(7)) + &
         q(5)*(11003.0*q(5) - 17246.0*q(6) + 4642.0*q(7)) &
         + q(6)*(7043.0*q(6) - 3882.0*q(7)) + 547.0*q(7)**2
  P3 = (3.0*q(4) + 13.0*q(5) - 5.0*q(6) + q(7))/12.0

  ! Alpha values
  !tau = abs(b3 - 3.0 * b2 + 3.0 * b1 - b0)
  tau = abs(b3 - b0)
  w0 = d0*fac_fn(tau, b0)
  w1 = d1*fac_fn(tau, b1)
  w2 = d2*fac_fn(tau, b2)
  w3 = d3*fac_fn(tau, b3)

  ! Normalization
  wnorm = w0+w1+w2+w3
  wpL = (w0*P0 + w1*P1 + w2*P2 + w3*P3)/wnorm

  ! Compute flux at the right side of i-1/2
  !d0 = 4.0/35.0 ;  d1 = 18.0/35.0 ; d2 = 12.0/35.0 ; d3 = 1.0/35.0
  ! 1st stencil
  P0 = (-3.0*q(7) + 13.0*q(6) - 23.0*q(5) + 25.0*q(4))/12.0
  b0 = q(7)*(547.0*q(7) - 3882.0*q(6) + 4642.0*q(5) - 1854.0*q(4)) + &
        q(6)*(7043.0*q(6) - 17246.0*q(5) + 7042.0*q(4)) + &
        q(5)*(11003.0*q(5) - 9402.0*q(4)) + 2107.0*q(4)**2

  ! 2nd stencil
  P1 = ( q(6) - 5.0*q(5) + 13.0*q(4) + 3.0*q(3) )/12.0
  b1 = q(6)*(267.0*q(6) - 1642.0*q(5) + 1602.0*q(4) - 494.0*q(3)) + &
         q(5)*(2843.0*q(5) - 5966.0*q(4) + 1922.0*q(3)) + &
         q(4)*(3443.0*q(4) - 2522.0*q(3)) + 547.0*q(3)**2

  ! 3rd stencil
  P2 = ( -q(5) + 7.0*q(4) + 7.0*q(3) - q(2) )/12.0
  b2 = q(5)*(547.0*q(5) - 2522.0*q(4) + 1922.0*q(3) - 494.0*q(2)) + &
         q(4)*(3443.0*q(4) - 5966.0*q(3) + 1602.0*q(2)) + &
         q(3)*(2843.0*q(3) - 1642.0*q(2)) + 267.0*q(2)**2

  ! 4th stencil
  P3 = ( 3.0*q(4) + 13.0*q(3) - 5.0*q(2) + q(1) )/12.0
  b3 = q(4)*(2107.0*q(4) - 9402.0*q(3) + 7042.0*q(2) - 1854.0*q(1)) + &
        q(3)*(11003.0*q(3) - 17246.0*q(2) + 4642.0*q(1)) + &
        q(2)*(7043.0*q(2) - 3882.0*q(1)) + 547.0*q(1)**2

  ! Alpha values
  !tau = abs(b3 - 3.0 * b2 + 3.0 * b1 - b0)
  tau = abs(b3 - b0)
  w0 = d0*fac_fn(tau, b0)
  w1 = d1*fac_fn(tau, b1)
  w2 = d2*fac_fn(tau, b2)
  w3 = d3*fac_fn(tau, b3)

  ! Normalization
  wnorm = w0+w1+w2+w3
  wmR = (w0*P0 + w1*P1 + w2*P2 + w3*P3)/wnorm

  !call PP_limiter(q0, wmR, wpL, h_min)

end subroutine weno7_reconstruction_interface

!> This is the subroutine for the positivity-preserving limiter
!! It limits the WENO reconstruction to give a reconstruction
!! that is positive-definite.
subroutine PP_limiter(q0, wmr, wpl, h_min)
  real, intent(in) :: q0 !< tracer concentration in cell i
  real, intent(inout) :: wmr, wpl   !< weno reconstruction on the cell interface i-1/2 and i+1/2
  real, intent(in)  :: h_min     !< The minimum thickness

  real :: qmin, qmax, theta, eps
  real :: w0, P0

  w0 = 5.0/18.0
  P0 = (q0 - w0*(wmr + wpl))/(1.0 - 2.0*w0)
  qmin = min(wmr, P0, wpl)

  eps = min(h_min, q0)
  theta = min(((q0-eps)/(q0-qmin)), 1.0)
  wpl = theta*(wpl - q0) + q0
  wmr = theta*(wmr - q0) + q0

  end subroutine PP_limiter

!> This is the subroutine for the positivity-preserving limiter
!! It limits the WENO reconstruction to give a reconstruction
!! that is positive-definite.
subroutine WENO_limiter(h_in, wmr, wpl, h_min, G, iis, iie, jis, jie)
  type(ocean_grid_type),             intent(in)  :: G    !< Ocean's grid structure.
  real, dimension(SZI_(G),SZJ_(G)),  intent(in)  :: h_in !< Layer thickness [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G)),  intent(inout) :: wmr !< Left thickness in the reconstruction [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G)),  intent(inout) :: wpl !< Right thickness in the reconstruction [H ~> m or kg m-2].
  real,                              intent(in)  :: h_min !< The minimum thickness
                    !! that can be obtained by a concave parabolic fit [H ~> m or kg m-2]
  integer,                           intent(in)  :: iis      !< Start of i index range.
  integer,                           intent(in)  :: iie      !< End of i index range.
  integer,                           intent(in)  :: jis      !< Start of j index range.
  integer,                           intent(in)  :: jie      !< End of j index range.

  real :: qmin, theta, eps, a(3)
  !real, parameter :: Fmin(3) = (/ 1.0,  -0.5,  0.0 /)
  !real, parameter :: Fmax(3) = (/ 1.0,   0.5,  0.25 /)
  integer :: i,j,k
  real :: w0, P0

  w0 = 5.0/18.0

  do j=jis,jie ; do i=iis,iie
    ! This limiter prevents undershooting minima within the domain with
    ! values less than h_min.
    eps = min(h_min, h_in(i,j))

    !a(1) = (6.0*h_in(i,j) - (wmr(i,j) + wpl(i,j)))/4.0
    !a(2) = (wpl(i,j) - wmr(i,j))
    !a(3) = -6.0*h_in(i,j) + 3.0*(wmr(i,j) + wpl(i,j))

    !qmin = 0.0
    !do k = 1,3
    !  qmin = qmin + a(k)*(0.5*(1.0-sign(1.0,a(k)))*Fmax(k) + 0.5*(1.0+sign(1.0,a(k)))*Fmin(k))
    !enddo

    P0 = (h_in(i,j) - w0*(wmr(i,j) + wpl(i,j)))/(1.0 - 2.0*w0)
    qmin = min(wmr(i,j), P0, wpl(i,j))

    theta = min(((h_in(i,j)-eps)/(h_in(i,j)-qmin)), 1.0)
    wpl(i,j) = theta*(wpl(i,j) - h_in(i,j)) + h_in(i,j)
    wmr(i,j) = theta*(wmr(i,j) - h_in(i,j)) + h_in(i,j)

  enddo ; enddo

end subroutine WENO_limiter

!> Compute the factor for the WENO weights
function fac_fn(tau, b) result(fac)
  real, intent(in)  :: tau  !< Difference of the smoothness indicator [A ~> a]
  real, intent(in)  :: b    !< The smoothness indicator [A ~> a]
  real :: fac               !< The factor for the weight [nondim]

  fac = 1.0e40; if (abs(b) > 1.0e-20*tau) fac = (1 + tau / b)**2

end function fac_fn

!> \namespace mom_continuity_weno
!!
!! This module contains the subroutines that advect layer
!! thickness.  The scheme here uses Weighted Essentially Non-Oscillatory Schemes with
!! a positive definite limiter.

end module MOM_continuity_WENO
