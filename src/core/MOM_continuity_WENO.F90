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
subroutine weno3_reconstruction_interface(wmR, wpL, qm, q0, qp, h_min)

   real, intent(in) :: qm, q0, qp !< tracer concentration for 3-stencil wide
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
   P1 = 0.5*(-qm + 3.0*q0)
   P2 = 0.5*(q0 + qp)

   ! smoothness indicator
   b1 = (q0-qm)*(q0-qm)
   b2 = (qp-q0)*(qp-q0)

   ! Alpha values
   tau = abs(b1-b2)
   w1 = d1*(1.0 + (tau/(b1+eps))**r)
   w2 = d2*(1.0 + (tau/(b2+eps))**r)

   ! Normalization
   wnorm = w1+w2
   wpL = (w1*P1 + w2*P2)/wnorm

   ! Compute flux at the right side of i-1/2
   ! reconstructed polynomials
   P1 = 0.5*(-qp + 3.0*q0)
   P2 = 0.5*(q0 + qm)

   ! smoothness indicator
   b1 = (q0-qp)*(q0-qp)
   b2 = (qm-q0)*(qm-q0)

   ! Alpha values
   tau = abs(b1-b2)
   w1 = d1*(1.0 + (tau/(b1+eps))**r)
   w2 = d2*(1.0 + (tau/(b2+eps))**r)

   ! Normalization
   wnorm = w1+w2
   wmR = (w1*P1 + w2*P2)/wnorm

   !call PP_limiter(q0, wmR, wpL, h_min)

end subroutine weno3_reconstruction_interface

!> 5th-order weno z-type reconstruction flux
subroutine weno5_reconstruction_interface(wmR, wpL, qmm, qm, q0, qp, qpp, h_min)

   real, intent(in) :: qmm, qm, q0, qp, qpp !< tracer concentration for 5-stencil wide
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
   P0 = (2.0*qmm - 7.0*qm + 11.0*q0)/6.0
   b0 = (13.0/12.0)*(qmm - 2.0*qm + q0)**2 + 0.25*(qmm - 4.0*qm + 3.0*q0)**2

   ! Second stencil
   P1 = (-qm + 5.0*q0 + 2.0*qp)/6.0
   b1 = (13.0/12.0)*(qm - 2.0*q0 + qp)**2 + 0.25*(qm - qp)**2

   ! Third stencil
   P2 = (2.0*q0 + 5.0*qp - qpp)/6.0
   b2 = (13.0/12.0)*(q0 - 2.0*qp + qpp)**2 + 0.25*(3.0*q0 - 4.0*qp + qpp)**2

   ! Alpha values
   tau = abs(b2-b0)
   w0 = d0*(1.0 + (tau/(b0+eps))**r)
   w1 = d1*(1.0 + (tau/(b1+eps))**r)
   w2 = d2*(1.0 + (tau/(b2+eps))**r)

   wnorm = w0+w1+w2
   wpL = (w0*P0 + w1*P1 + w2*P2)/wnorm

   ! Compute flux at the right side of i-1/2
   ! First stencil
   P0 = (2.0*qpp - 7.0*qp + 11.0*q0)/6.0
   b0 = (13.0/12.0)*(qpp - 2.0*qp + q0)**2 + 0.25*(qpp - 4.0*qp + 3.0*q0)**2

   ! Second stencil
   P1 = (-qp + 5.0*q0 + 2.0*qm)/6.0
   b1 = (13.0/12.0)*(qp - 2.0*q0 + qm)**2 + 0.25*(qp - qm)**2

   ! Third stencil
   P2 = (2.0*q0 + 5.0*qm - qmm)/6.0
   b2 = (13.0/12.0)*(q0 - 2.0*qm + qmm)**2 + 0.25*(3.0*q0 - 4.0*qm + qmm)**2

   ! Alpha values
   tau = abs(b2-b0)
   w0 = d0*(1.0 + (tau/(b0+eps))**r)
   w1 = d1*(1.0 + (tau/(b1+eps))**r)
   w2 = d2*(1.0 + (tau/(b2+eps))**r)

   wnorm = w0+w1+w2
   wmR = (w0*P0 + w1*P1 + w2*P2)/wnorm

   !call PP_limiter(q0, wmR, wpL, h_min)

end subroutine weno5_reconstruction_interface

!> 7th-order weno z-type reconstruction flux
subroutine weno7_reconstruction_interface(wmR, wpL, qm3, qm2, qm1, q0, qp1, qp2, qp3, h_min)

   real, intent(in) :: qm3, qm2, qm1, q0, qp1, qp2, qp3 !< tracer concentration for 7-stencil wide
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
   P0 = (-3.0*qm3 + 13.0*qm2 - 23.0*qm1 + 25.0*q0)/12.0
   b0 = qm3*(547.0*qm3 - 3882.0*qm2 + 4642.0*qm1 - 1854.0*q0) + &
        qm2*(7043.0*qm2 - 17246.0*qm1 + 7042.0*q0) + &
        qm1*(11003.0*qm1 - 9402.0*q0) + 2107.0*q0**2

   ! 2nd stencil
   P1 = (qm2 - 5.0*qm1 + 13.0*q0 + 3.0*qp1)/12.0
   b1 = qm2*(267.0*qm2 - 1642.0*qm1 + 1602.0*q0 - 494.0*qp1) + &
           qm1*(2843.0*qm1 - 5966.0*q0 + 1922.0*qp1) &
           + q0*(3443.0*q0 - 2522.0*qp1) + 547.0*qp1**2

   ! 3rd stencil
   P2 = (-qm1 + 7.0*q0 + 7.0*qp1 - qp2)/12.0
   b2 = qm1*(547.0*qm1 - 2522.0*q0 + 1922.0*qp1 - 494.0*qp2) + &
           q0*(3443.0*q0 - 5966.0*qp1 + 1602.0*qp2) &
           + qp1*(2843.0*qp1 - 1642.0*qp2) + 267.0*qp2**2

   ! 4rd stencil
   P3 = (3.0*q0 + 13.0*qp1 - 5.0*qp2 + qp3)/12.0
   b3 = q0*(2107.0*q0 - 9402.0*qp1 + 7042.0*qp2 - 1854.0*qp3) + &
           qp1*(11003.0*qp1 - 17246.0*qp2 + 4642.0*qp3) &
           + qp2*(7043.0*qp2 - 3882.0*qp3) + 547.0*qp3**2

   ! Alpha values
   tau = abs(b3 + 3.0 * b2 - 3.0 * b1 - b0)
   !tau = abs(b3 - b0)
   w0 = d0*(1.0 + (tau/(b0+eps))**r)
   w1 = d1*(1.0 + (tau/(b1+eps))**r)
   w2 = d2*(1.0 + (tau/(b2+eps))**r)
   w3 = d3*(1.0 + (tau/(b3+eps))**r)

   ! Normalization
   wnorm = w0+w1+w2+w3
   wpL = (w0*P0 + w1*P1 + w2*P2 + w3*P3)/wnorm

   ! Compute flux at the right side of i-1/2
   ! 1st stencil
   P0 = (-3.0*qp3 + 13.0*qp2 - 23.0*qp1 + 25.0*q0)/12.0
   b0 = qp3*(547.0*qp3 - 3882.0*qp2 + 4642.0*qp1 - 1854.0*q0) + &
        qp2*(7043.0*qp2 - 17246.0*qp1 + 7042.0*q0) + &
        qp1*(11003.0*qp1 - 9402.0*q0) + 2107.0*q0**2

   ! 2nd stencil
   P1 = (qp2 - 5.0*qp1 + 13.0*q0 + 3.0*qm1)/12.0
   b1 = qp2*(267.0*qp2 - 1642.0*qp1 + 1602.0*q0 - 494.0*qm1) + &
           qp1*(2843.0*qp1 - 5966.0*q0 + 1922.0*qm1) &
           + q0*(3443.0*q0 - 2522.0*qm1) + 547.0*qm1**2

   ! 3rd stencil
   P2 = (-qp1 + 7.0*q0 + 7.0*qm1 - qm2)/12.0
   b2 = qp1*(547.0*qp1 - 2522.0*q0 + 1922.0*qm1 - 494.0*qm2) + &
           q0*(3443.0*q0 - 5966.0*qm1 + 1602.0*qm2) &
           + qm1*(2843.0*qm1 - 1642.0*qm2) + 267.0*qm2**2

   ! 4rd stencil
   P3 = (3.0*q0 + 13.0*qm1 - 5.0*qm2 + qm3)/12.0
   b3 = q0*(2107.0*q0 - 9402.0*qm1 + 7042.0*qm2 - 1854.0*qm3) + &
           qm1*(11003.0*qm1 - 17246.0*qm2 + 4642.0*qm3) &
           + qm2*(7043.0*qm2 - 3882.0*qm3) + 547.0*qm3**2

   ! Alpha values
   tau = abs(b3 + 3.0 * b2 - 3.0 * b1 - b0)
   !tau = abs(b3 - b0)
   w0 = d0*(1.0 + (tau/(b0+eps))**r)
   w1 = d1*(1.0 + (tau/(b1+eps))**r)
   w2 = d2*(1.0 + (tau/(b2+eps))**r)
   w3 = d3*(1.0 + (tau/(b3+eps))**r)

   ! Normalization
   wnorm = w0+w1+w2+w3
   wmR = (w0*P0 + w1*P1 + w2*P2 + w3*P3)/wnorm

   !call PP_limiter(q0, wmR, wpL, h_min)

end subroutine weno7_reconstruction_interface

!> This is the subroutine for the positivity-preserving limiter
!! It limits the WENO reconstruction to give a reconstruction
!! that is positive-definite.
subroutine PP_limiter0(q0, wmr, wpl, h_min)

   real, intent(in) :: q0 !< tracer concentration in cell i
   real, intent(inout) :: wmr, wpl   !< weno reconstruction on the cell interface i-1/2 and i+1/2
   real, intent(in)  :: h_min     !< The minimum thickness

   real :: qmin, theta, eps, a(3)
   real :: theta0, theta1
   integer :: i
   real, parameter :: Fmin(3) = (/ 1.0,  -0.5,  0.0 /)
   real, parameter :: Fmax(3) = (/ 1.0,   0.5,  0.25 /)

   eps = min(h_min, q0)

   a(1) = (6.0*q0 - (wmr + wpl))/4.0
   a(2) = (wpl - wmr)
   a(3) = -6.0*q0 + 3.0*(wmr + wpl)

   qmin = 0.0
   do i = 1,3
     qmin = qmin + a(i)*(0.5*(1.0-sign(1.0,a(i)))*Fmax(i) + 0.5*(1.0+sign(1.0,a(i)))*Fmin(i))
   enddo

   theta = min(((q0-eps)/(q0-qmin)), 1.0)
   wpl = theta*(wpl - q0) + q0
   wmr = theta*(wmr - q0) + q0

end subroutine PP_limiter0

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
  real, parameter :: Fmin(3) = (/ 1.0,  -0.5,  0.0 /)
  real, parameter :: Fmax(3) = (/ 1.0,   0.5,  0.25 /)
  integer :: i,j,k
  real :: w0, P0

  w0 = 5.0/18.0

  do j=jis,jie ; do i=iis,iie
    ! This limiter prevents undershooting minima within the domain with
    ! values less than h_min.
    eps = min(h_min, h_in(i,j))

    a(1) = (6.0*h_in(i,j) - (wmr(i,j) + wpl(i,j)))/4.0
    a(2) = (wpl(i,j) - wmr(i,j))
    a(3) = -6.0*h_in(i,j) + 3.0*(wmr(i,j) + wpl(i,j))

    qmin = 0.0
    do k = 1,3
      qmin = qmin + a(k)*(0.5*(1.0-sign(1.0,a(k)))*Fmax(k) + 0.5*(1.0+sign(1.0,a(k)))*Fmin(k))
    enddo

    !P0 = (h_in(i,j) - w0*(wmr(i,j) + wpl(i,j)))/(1.0 - 2.0*w0)
    !qmin = min(wmr(i,j), P0, wpl(i,j))

    theta = min(((h_in(i,j)-eps)/(h_in(i,j)-qmin)), 1.0)
    wpl(i,j) = theta*(wpl(i,j) - h_in(i,j)) + h_in(i,j)
    wmr(i,j) = theta*(wmr(i,j) - h_in(i,j)) + h_in(i,j)

  enddo ; enddo

end subroutine WENO_limiter

!> \namespace mom_continuity_weno
!!
!! This module contains the subroutines that advect layer
!! thickness.  The scheme here uses Weighted Essentially Non-Oscillatory Schemes with
!! a positive definite limiter.

end module MOM_continuity_WENO
