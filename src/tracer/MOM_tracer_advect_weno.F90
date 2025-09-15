!>  This module contains the subroutines of the WENO schemes that advect tracers along coordinate surfaces.
module MOM_tracer_advect_weno

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_grid,            only : ocean_grid_type
use MOM_tracer_registry, only : tracer_registry_type, tracer_type
use MOM_verticalGrid,    only : verticalGrid_type
use MOM_spatial_means, only : array_global_min_max
implicit none ; private

#include <MOM_memory.h>

public weno3_reconstruction
public weno5_reconstruction
public weno7_reconstruction
public weno9_reconstruction
public tracer_min_max_init
public PPM_reconstruction
public weno5NM_reconstruction

contains

!> 3rd weno reconstruction subroutine and limiter
pure subroutine weno3_reconstruction(wq, q, u, mu)
  real, intent(in) :: q(3) !< tracer concentration from cell i-1 to i+1
  real, intent(in) :: u           !< advection velocity
  real, intent(in) :: mu         !< cfl
  real, intent(out) :: wq         !< weno reconstruction at the cell
                                  !! interface i+1/2

  real :: wmr ! wmr : weno reconstruction on the cell interface i-1/2
  real :: wpl ! wpl : weno reconstruction on the cell interface i+1/2

  wq = 0.0

  call weno3_reconstruction_interface(wpl, wmr, q(1), q(2), q(3))
  ! Monotonicity limiter
  call WENO_monotonicity_limiter(q(1), q(2), q(3), wmr, wpl, mu)

  if (u >= 0.0) then
    wq = wpl
  else
    wq = wmr
  endif

end subroutine weno3_reconstruction

!> 3rd-order weno z-type reconstruction flux
pure subroutine weno3_reconstruction_interface(wpl, wmr, qm, q0, qp)
  real, intent(in) :: qm, q0, qp !< tracer concentration for 5-stencil wide
  real, intent(out) :: wpl, wmr  !< reconstructed weno flux

  real :: b1, b2    ! smoothness indicator
  real :: d1, d2    ! linear weights
  real :: P1, P2    ! reconstructed polynomials
  real :: w1, w2    ! nonlinear weights
  real :: a1, a2
  real :: eps, wnorm, tau
  integer, parameter :: r = 1

  ! linear weights
  d1 = 1.0/3.0 ; d2 = 2.0/3.0
  eps = 1.0e-20

  ! Compute flux at left side of i+1/2
  ! reconstructed polynomials
  P1 = 0.5*(-qm + 3.0*q0)
  P2 = 0.5*(q0 + qp)

  ! smoothness indicator
  b1 = (q0-qm)*(q0-qm)
  b2 = (qp-q0)*(qp-q0)

  ! Alpha values
  tau = abs(b2-b1)
  w1 = d1*(1.0 + (tau/(b1+eps))**r)
  w2 = d2*(1.0 + (tau/(b2+eps))**r)

  ! Normalization
  wnorm = w1+w2
  wpl = (w1*P1 + w2*P2)/wnorm
  !if ((qp-q0)*(q0-qm) <= 0.) wpl = q0
  wpl = max(min(q0,qp), wpl) ; wpl = min(max(q0,qp), wpl)

  ! Compute flux at the right side of i-1/2
  ! reconstructed polynomials
  P1 = 0.5*(-qp + 3.0*q0)
  P2 = 0.5*(q0 + qm)

  ! smoothness indicator
  b1 = (q0-qp)*(q0-qp)
  b2 = (qm-q0)*(qm-q0)

  ! Alpha values
  tau = abs(b2-b1)
  w1 = d1*(1.0 + (tau/(b1+eps))**r)
  w2 = d2*(1.0 + (tau/(b2+eps))**r)

  ! Normalization
  wnorm = w1+w2
  wmr = (w1*P1 + w2*P2)/wnorm
  !if ((qp-q0)*(q0-qm) <= 0.) wmr = q0
  wmr = max(min(q0,qm), wmr) ; wmr = min(max(q0,qm), wmr)

end subroutine weno3_reconstruction_interface

!> 5th-order weno reconstruction subroutine and limiter
pure subroutine weno5_reconstruction(wq, q, u, mu)
  real, intent(in) :: q(5)       !< tracer concentration from i-2 to  i+3 respectively
  real, intent(in) :: u          !< advection velocity
  real, intent(in) :: mu         !< cfl
  real, intent(out) :: wq        !< weno reconstruction at the interface i+1/2

  real :: wmr ! wmr : weno reconstruction on the cell interface i-1/2
  real :: wpl ! wpl : weno reconstruction on the cell interface i+1/2

  wq = 0.0

  call weno5z_reconstruction_interface(wpl, wmr, q(1), q(2), q(3), q(4), q(5))
  ! Monotonicity limiter
  call WENO_monotonicity_limiter(q(2), q(3), q(4), wmr, wpl, mu)

  if (u >= 0.0) then
    wq = wpl
  else
    wq = wmr
  endif

end subroutine weno5_reconstruction

!> 5th-order weno z-type reconstruction flux
pure subroutine weno5z_reconstruction_interface(wpl, wmr, qmm, qm, q0, qp, qpp)
  real, intent(in) :: qmm, qm, q0, qp, qpp !< tracer concentration for 5-stencil wide
  real, intent(out) :: wpl, wmr            !< reconstructed weno flux

  real :: P0, P1, P2         ! reconstructed polynomials
  real :: b0, b1, b2         ! smoothness indicator
  real :: w0, w1, w2         ! nonlinear weights
  real :: d0, d1, d2         ! linear weights
  real :: a0, a1, a2
  real :: eps,  wnorm, tau
  integer, parameter :: r = 1

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
  w0 = w0/wnorm ; w1 = w1/wnorm ; w2 = w2/wnorm

  wpl = w0*P0 + w1*P1 + w2*P2
  wpl = max(min(q0,qp), wpl) ; wpl = min(max(q0,qp), wpl)

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
  w0 = w0/wnorm ; w1 = w1/wnorm ; w2 = w2/wnorm

  wmr = w0*P0 + w1*P1 + w2*P2
  wmr = max(min(q0,qm), wmr) ; wmr = min(max(q0,qm), wmr)

end subroutine weno5z_reconstruction_interface

!> 5th-order weno reconstruction for non-uniform grid and limiter
pure subroutine weno5NM_reconstruction(wq, q, u, ds, mu)
  real, intent(in) :: q(5)    !< tracer concentration from i-2 to  i+3 respectively
  real, intent(in) :: u       !< advection velocity
  real, intent(in) :: ds(5)   !< grid sizes
  real, intent(in) :: mu         !< cfl
  real, intent(out) :: wq     !< weno reconstruction at the interface i+1/2

  real :: wmr ! wmr : weno reconstruction on the cell interface i-1/2
  real :: wpl ! wpl : weno reconstruction on the cell interface i+1/2

  call weno5NM_reconstruction_interface(wpl, wmr, q(1), q(2), q(3), q(4), q(5), ds)
  ! Monotonicity limiter
  call WENO_monotonicity_limiter(q(2), q(3), q(4), wmr, wpl, mu)

  if (u >= 0.0) then
    wq = wpl
  else
    wq = wmr
  endif

end subroutine weno5NM_reconstruction

!> 5th-order weno z-type reconstruction flux for non-uniform grid
pure subroutine weno5NM_reconstruction_interface(wpl, wmr, qmm, qm, q0, qp, qpp, dx)
  real, intent(in) :: qmm, qm, q0, qp, qpp !< tracer concentration for 5-stencil wide
  real, intent(in) :: dx(5)                !< grid sizes
  real, intent(out) :: wpl, wmr                  !< reconstructed weno flux

  real :: P0, P1, P2         ! reconstructed polynomials
  real :: b0, b1, b2         ! smoothness indicator
  real :: w0, w1, w2         ! nonlinear weights
  real :: d0, d1, d2         ! linear weights
  real :: a0, a1, a2, a3
  real :: eps,  wnorm, tau
  integer, parameter :: r = 2
  real :: qh, qhh, qhp, qhpp

  ! Gamma values in Weno reconstruction
  d0 = 1.0/10.0 ; d1 = 6.0/10.0 ; d2 = 3.0/10.0
  eps = 1.0e-6

  ! Compute flux at the right side of i+1/2
  ! First stencil
  a0 = dx(1)/(dx(1)+dx(2)) ; a1 = dx(2)/(dx(2)+dx(3))

  qh = (1.0-a1)*qm + a1*q0
  qhh = (2.0-a0)*qm - (1.0-a0)*qmm

  b0 = (13.0/12.0)*(2.0*qh - 2.0*qhh)**2 + 0.25*(4.0*q0 - 2.0*qh - 2.0*qhh)**2
  P0 = (6.0*q0 - qh - 2.0*qhh)/3.0

  ! Second stencil
  a2 = dx(3)/(dx(3)+dx(4))
  qhp = (1.0-a2)*q0 + a2*qp

  b1 = (13.0/12.0)*(2.0*qh - 4.0*q0 + 2.0*qhp)**2 + 0.25*(-2.0*qh + 2.0*qhp)**2
  P1 = (-qh + 2*q0 + 2.0*qhp)/3.0

  ! Third stencil
  a3 = dx(4)/(dx(4)+dx(5))
  qhpp = (1.0-a3)*qp + a3*qpp

  b2 = (13.0/12.0)*(2.0*qhp - 4.0*qp + 2.0*qhpp)**2 + 0.25*(-6.0*qhp + 8.0*qp - 2.0*qhpp)**2
  P2 = (2.0*qhp + 2.0*qp - qhpp)/3.0

  ! Alpha values
  tau = abs(b2-b0)
  w0 = d0*(1.0 + (tau/(b0+eps))**r)
  w1 = d1*(1.0 + (tau/(b1+eps))**r)
  w2 = d2*(1.0 + (tau/(b2+eps))**r)

  wnorm = w0+w1+w2
  wpl = (w0*P0 + w1*P1 + w2*P2)/wnorm
  wpl = max(min(q0,qp), wpl) ; wpl = min(max(q0,qp), wpl)

  ! Compute flux at the right side of i-1/2
  ! First stencil
  a0 = dx(5)/(dx(5)+dx(4)) ; a1 = dx(4)/(dx(4)+dx(3))

  qh = (1.0-a1)*qp + a1*q0
  qhh = (2.0-a0)*qp - (1.0-a0)*qpp

  b0 = (13.0/12.0)*(2.0*qh - 2.0*qhh)**2 + 0.25*(4.0*q0 - 2.0*qh - 2.0*qhh)**2
  P0 = (6.0*q0 - qh - 2.0*qhh)/3.0

  ! Second stencil
  a2 = dx(3)/(dx(3)+dx(2))
  qhp = (1.0-a2)*q0 + a2*qm

  b1 = (13.0/12.0)*(2.0*qh - 4.0*q0 + 2.0*qhp)**2 + 0.25*(-2.0*qh + 2.0*qhp)**2
  P1 = (-qh + 2*q0 + 2.0*qhp)/3.0

  ! Third stencil
  a3 = dx(2)/(dx(2)+dx(1))
  qhpp = (1.0-a3)*qm + a3*qmm

  b2 = (13.0/12.0)*(2.0*qhp - 4.0*qm + 2.0*qhpp)**2 + 0.25*(-6.0*qhp + 8.0*qm - 2.0*qhpp)**2
  P2 = (2.0*qhp + 2.0*qm - qhpp)/3.0

  ! Alpha values
  tau = abs(b2-b0)
  w0 = d0*(1.0 + (tau/(b0+eps))**r)
  w1 = d1*(1.0 + (tau/(b1+eps))**r)
  w2 = d2*(1.0 + (tau/(b2+eps))**r)

  wnorm = w0+w1+w2
  wmr = (w0*P0 + w1*P1 + w2*P2)/wnorm
  wmr = max(min(q0,qm), wmr) ; wmr = min(max(q0,qm), wmr)

end subroutine weno5NM_reconstruction_interface

!> 7th-order weno reconstruction subroutine and limiter
pure subroutine weno7_reconstruction(wq, q, u, mu)
  real, intent(in) :: q(7)       !< tracer concentration
                                 !! from i-3 to i+4 respectively
  real, intent(in) :: u          !< advection velocity
  real, intent(in) :: mu         !< cfl
  real, intent(out) :: wq        !< weno reconstruction at the interface i+1/2

  real :: wmr ! wmr : weno reconstruction on the cell interface i-1/2
  real :: wpl ! wpl : weno reconstruction on the cell interface i+1/2

  wq = 0.0

  call weno7z_reconstruction_interface(wpl, wmr, q(1), q(2), q(3), q(4), q(5), q(6), q(7))
  ! Monotonicity limiter
  call WENO_monotonicity_limiter(q(3), q(4), q(5), wmr, wpl, mu)

  if (u >= 0.0) then
    wq = wpl
  else
    wq = wmr
  endif

end subroutine weno7_reconstruction

!> 7th-order weno z-type reconstruction flux
pure subroutine weno7z_reconstruction_interface(wpl, wmr, qm3, qm2, qm1, q0, qp1, qp2, qp3)
  real, intent(in) :: qm3, qm2, qm1, q0, qp1, qp2, qp3 !< tracer concentration for 7-stencil wide
  real, intent(out) :: wpl, wmr                        !< reconstructed weno flux

  real :: P0, P1, P2, P3     ! reconstructed polynomials
  real :: b0, b1, b2, b3     ! smoothness indicator
  real :: w0, w1, w2, w3     ! nonlinear weights
  real :: d0, d1, d2, d3     ! nonlinear weights
  real :: a0, a1, a2, a3
  real :: eps, tau, wnorm
  integer, parameter :: r = 1

  d0 = 1.0/35.0 ;  d1 = 12.0/35.0 ; d2 = 18.0/35.0 ; d3 = 4.0/35.0
  eps = 1.0e-6

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
  w0 = d0*(1.0 + (tau/(b0+eps))**r)
  w1 = d1*(1.0 + (tau/(b1+eps))**r)
  w2 = d2*(1.0 + (tau/(b2+eps))**r)
  w3 = d3*(1.0 + (tau/(b3+eps))**r)

  ! Normalization
  wnorm = w0+w1+w2+w3
  w0 = w0/wnorm ; w1 = w1/wnorm ; w2 = w2/wnorm ; w3 = w3/wnorm

  wpl = w0*P0 + w1*P1 + w2*P2 + w3*P3
  wpl = max(min(q0,qp1), wpl) ; wpl = min(max(q0,qp1), wpl)

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
  w0 = d0*(1.0 + (tau/(b0+eps))**r)
  w1 = d1*(1.0 + (tau/(b1+eps))**r)
  w2 = d2*(1.0 + (tau/(b2+eps))**r)
  w3 = d3*(1.0 + (tau/(b3+eps))**r)

  ! Normalization
  wnorm = w0+w1+w2+w3
  w0 = w0/wnorm ; w1 = w1/wnorm ; w2 = w2/wnorm ; w3 = w3/wnorm

  wmr = w0*P0 + w1*P1 + w2*P2 + w3*P3
  wmr = max(min(q0,qm1), wmr) ; wmr = min(max(q0,qm1), wmr)

end subroutine weno7z_reconstruction_interface

!> 9th-order weno reconstruction subroutine and limiter
pure subroutine weno9_reconstruction(wq, q, u, mu)
  real, intent(in) :: q(9)      !< tracer concentration
                                 !! from i-4 to i+4
  real, intent(in) :: u          !< advection velocity
  real, intent(in) :: mu         !< cfl
  real, intent(out) :: wq        !< weno reconstruction at the interface i+1/2

  real :: wmr ! wmr : weno reconstruction on the cell interface i-1/2
  real :: wpl ! wpl : weno reconstruction on the cell interface i+1/2

  wq = 0.0

  call weno9_reconstruction_interface(wpl, wmr, q(1), q(2), q(3), q(4), &
          q(5), q(6), q(7), q(8), q(9))
  ! Monotonicity limiter
  call WENO_monotonicity_limiter(q(4), q(5), q(6), wmr, wpl, mu)

  if (u >= 0.0) then
    wq = wpl
  else
    wq = wmr
  endif

end subroutine weno9_reconstruction

!> 9th-order weno z-type reconstruction flux
pure subroutine weno9_reconstruction_interface(wpl, wmr, qm4, qm3, qm2, qm1, &
                q0, qp1, qp2, qp3, qp4)
  real, intent(in) :: qm4, qm3, qm2, qm1, q0, qp1, qp2, qp3, qp4 !< tracer concentration
                                                                 !! for 7-stencil wide
  real, intent(out) :: wpl, wmr                                  !< reconstructed weno flux

  real :: b0, b1, b2, b3, b4               ! smoothness indicator
  real :: d0, d1, d2, d3, d4               ! linear weights
  real :: w0, w1, w2, w3, w4               ! nonlinear weights
  real :: P0, P1, P2, P3, P4               ! reconstructed polynomials
  real :: a0, a1, a2, a3, a4
  real :: eps,wnorm, tau
  integer, parameter :: r = 2

  d0 = 1.0/126.0 ; d1 = 10.0/63.0 ; d2 = 10.0/21.0 ; d3 = 20.0/63.0 ; d4 = 5.0/126.0
  eps = 1.0e-6

  ! Compute flux at the right side of i+1/2
  call weno9_poly(P0, P1, P2, P3, P4, b0, b1, b2, b3, b4, &
          qm4, qm3, qm2, qm1, q0, qp1, qp2, qp3, qp4)

  ! Alpha values
  tau = abs(b0 - b4)
  w0 = d0*(1.0 + (tau/(b0+eps))**r)
  w1 = d1*(1.0 + (tau/(b1+eps))**r)
  w2 = d2*(1.0 + (tau/(b2+eps))**r)
  w3 = d3*(1.0 + (tau/(b3+eps))**r)
  w4 = d4*(1.0 + (tau/(b4+eps))**r)

  ! Normalization
  wnorm = w0+w1+w2+w3+w4
  wpl = (w0*P0 + w1*P1 + w2*P2 + w3*P3 + w4*P4)/wnorm
  wpl = max(min(q0,qp1), wpl) ; wpl = min(max(q0,qp1), wpl) ! Bound

  ! Compute flux at the right side of i+1/2
  call weno9_poly(P0, P1, P2, P3, P4, b0, b1, b2, b3, b4, &
          qp4, qp3, qp2, qp1, q0, qm1, qm2, qm3, qm4)

  ! Alpha values
  tau = abs(b0 - b4)
  w0 = d0*(1.0 + (tau/(b0+eps))**r)
  w1 = d1*(1.0 + (tau/(b1+eps))**r)
  w2 = d2*(1.0 + (tau/(b2+eps))**r)
  w3 = d3*(1.0 + (tau/(b3+eps))**r)
  w4 = d4*(1.0 + (tau/(b4+eps))**r)

  ! Normalization
  wnorm = w0+w1+w2+w3+w4
  wmr = (w0*P0 + w1*P1 + w2*P2 + w3*P3 + w4*P4)/wnorm
  wmr = max(min(q0,qm1), wmr) ; wmr = min(max(q0,qm1), wmr) ! Bound

end subroutine weno9_reconstruction_interface

pure subroutine weno9_poly(P0, P1, P2, P3, P4, b0, b1, b2, b3, b4, qm4, qm3, qm2, qm1, q0, &
                qp1, qp2, qp3, qp4)
  real, intent(in) :: qm4, qm3, qm2, qm1, q0, qp1, qp2, qp3, qp4 !< tracer concentration
                                                                 !! for 9-stencil wide
  real, intent(out) :: P0, P1, P2, P3, P4    !< recontructed polynomials
  real, intent(out) :: b0, b1, b2, b3, b4    !< smoothness indicator

  real :: qx, qx2, qx3, qx4, L1, L2, L3, L4

  L1 = 0.5 ; L2 = 1.0/6.0 ; L3 = 1.0/20.0 ; L4 = 1.0/70.0

  ! 1st stencil
  qx = (27.0*qm4 - 146.0*qm3 + 336.0*qm2 - 462.0*qm1 + 245.0*q0)/120.0
  qx2 = (25.0*qm4 - 128.0*qm3 + 262.0*qm2 - 240.0*qm1 + 81.0*q0)/56.0
  qx3 = (3.0*qm4 - 14.0*qm3 + 24.0*qm2 - 18.0*qm1 + 5.0*q0)/12.0
  qx4 = (qm4 - 4.0*qm3 + 6.0*qm2 - 4.0*qm1 + q0)/24.0

  P0 = q0 + qx*L1 + qx2*L2 + qx3*L3 + qx4*L4
  b0 = (qx+ 0.1*qx3)**2 + (13.0/3.0)*(qx2 + (123.0/455.0)*qx4)**2 + &
        (781.0/20.0)*(qx3)**2 + (1421461.0/2275.0)*((qx4**2)*(qx4**2))

  ! 2nd stencil
  qx = (-11.0*qm3 + 66.0*qm2 - 192.0*qm1 + 110.0*q0 + 27.0*qp1)/120.0
  qx2 = (-3.0*qm3 + 12.0*qm2 + 10.0*qm1 - 44.0*q0 + 25.0*qp1)/56.0
  qx3 = (qm3 - 6.0*qm2 + 12.0*qm1 - 10.0*q0 + 3.0*qp1)/12.0
  qx4 = (qm3 - 4.0*qm2 + 6.0*qm1 - 4.0*q0 + qp1)/24.0

  P1 = q0 + qx*L1 + qx2*L2 + qx3*L3 + qx4*L4
  b1 = (qx+ 0.1*qx3)**2 + (13.0/3.0)*(qx2 + (123.0/455.0)*qx4)**2 + &
         (781.0/20.0)*(qx3)**2 + (1421461.0/2275.0)*((qx4**2)*(qx4**2))

  ! 3rd stencil
  qx = (11.0*qm2 - 82.0*qm1 + 82.0*qp1 - 11.0*qp2)/120.0
  qx2 = (-3.0*qm2 + 40.0*qm1 - 74.0*q0 + 40.0*qp1 - 3.0*qp2)/56.0
  qx3 = (-qm2 + 2.0*qm1 - 2.0*qp1 + qp2)/12.0
  qx4 = (qm2 - 4.0*qm1 + 6.0*q0 - 4.0*qp1 + qp2)/24.0

  P2 = q0 + qx*L1 + qx2*L2 + qx3*L3 + qx4*L4
  b2 = (qx+ 0.1*qx3)**2 + (13.0/3.0)*(qx2 + (123.0/455.0)*qx4)**2 + &
         (781.0/20.0)*(qx3)**2 + (1421461.0/2275.0)*((qx4**2)*(qx4**2))

  ! 4th stencil
  qx = (-27.0*qm1 - 110.0*q0 + 192.0*qp1 - 66.0*qp2 + 11.0*qp3)/120.0
  qx2 = (25.0*qm1 - 44.0*q0 + 10.0*qp1 + 12.0*qp2 - 3.0*qp3)/56.0
  qx3 = (-3.0*qm1 + 10.0*q0 - 12.0*qp1 + 6.0*qp2 - qp3)/12.0
  qx4 = (qm1 - 4.0*q0 + 6.0*qp1 - 4.0*qp2 + qp3)/24.0

  P3 = q0 + qx*L1 + qx2*L2 + qx3*L3 + qx4*L4
  b3 = (qx+ 0.1*qx3)**2 + (13.0/3.0)*(qx2 + (123.0/455.0)*qx4)**2 + &
         (781.0/20.0)*(qx3)**2 + (1421461.0/2275.0)*((qx4**2)*(qx4**2))

  ! 5th stencil
  qx = (-245.0*q0 + 462.0*qp1 - 336.0*qp2 + 146.0*qp3 - 27.0*qp4)/120.0
  qx2 = (81.0*q0 - 240.0*qp1 + 262.0*qp2 - 128.0*qp3 + 25.0*qp4)/56.0
  qx3 = (-5.0*q0 + 18.0*qp1 - 24.0*qp2 + 14.0*qp3 - 3.0*qp4)/24.0
  qx4 = (q0 - 4.0*qp1 + 6.0*qp2 - 4.0*qp3 + qp4)/24.0

  P4 = q0 + qx*L1 + qx2*L2 + qx3*L3 + qx4*L4
  b4 = (qx+ 0.1*qx3)**2 + (13.0/3.0)*(qx2 + (123.0/455.0)*qx4)**2 + &
         (781.0/20.0)*(qx3)**2 + (1421461.0/2275.0)*((qx4**2)*(qx4**2))

end subroutine weno9_poly

!> ppm reconstruction flux
subroutine PPM_reconstruction(wq_ppm, q, u, mu, qext)
  real, intent(in) :: q(3) !< tracer concentration for 3-stencil wide
  real, intent(in) :: u               !< advection velocity
  real, intent(in) :: mu              !< cfl
  real, intent(in) :: qext            !< check local extrema
  real, intent(out) :: wq_ppm         !< reconstructed flux

  real :: aL, aR, dA, mA, a6
  real :: qm, q0, qp

  qm = q(1) ; q0 = q(2) ; qp = q(3)

  aL = ( 5.*q0 + ( 2.*qm - qp ) )/6. ! H3 estimate
  aL = max( min(q0,qm), aL) ; aL = min( max(q0,qm), aL) ! Bound
  aR = ( 5.*q0 + ( 2.*qp - qm ) )/6. ! H3 estimate
  aR = max( min(q0,qp), aR) ; aR = min( max(q0,qp), aR) ! Bound

  dA = aR - aL ; mA = 0.5*( aR + aL )
  if (qext*(qp-q0)*(q0-qm) <= 0.) then
    aL = q0 ; aR = q0 ! PCM for local extrema and boundary cells
  elseif ( dA*(q0-mA) > (dA*dA)/6. ) then
    aL = (3.*q0) - 2.*aR
  elseif ( dA*(q0-mA) < - (dA*dA)/6. ) then
    aR = (3.*q0) - 2.*aL
  endif

  a6 = 6.0*q0 - 3.0 * (aR + aL) ! Curvature

  if (u >= 0.0) then
    wq_ppm = aR - 0.5*mu*((aR-aL) - a6*(1.0 - 2.0/3.0 * mu))
  else
    wq_ppm = aL + 0.5*mu*((aR-aL) + a6*(1.0 - 2.0/3.0 * mu))
  endif

end subroutine PPM_reconstruction

!> This is the subroutine for the WENO limiter
pure subroutine WENO_monotonicity_limiter(qm, q0, qp, wmr, wpl, mu)
  real, intent(in) :: qm, q0, qp   !< tracer concentration
  real, intent(in) :: mu           !< cfl
  real, intent(inout) :: wmr, wpl  !< reconstructed weno fluxes at right and left
                                   !! of cell i

  real :: dA, mA, a6, wtmp

  dA = wpl - wmr ; mA = 0.5*(wpl+wmr)
  if ((qp-q0)*(q0-qm) <= 0.0) then
    wpl = q0 ; wmr = q0
  elseif ( dA*(q0 - mA) > (dA*dA)/6.0 ) then
    wmr = (3.0*q0) - 2.0*wpl
  elseif ( dA*(q0 - mA) < - (dA*dA)/6.0 ) then
    wpl = (3.0*q0) - 2.0*wmr
  endif

  a6 = 6.0*q0 - 3.0 * (wpl + wmr)
  wtmp = wpl - 0.5*mu*((wpl-wmr) - a6*(1.0 - 2.0/3.0 * mu))
  wmr = wmr + 0.5*mu*((wpl-wmr) + a6*(1.0 - 2.0/3.0 * mu))
  wpl = wtmp

end subroutine WENO_monotonicity_limiter

!> Find the min and max of the tracer to use in max-principle limiter with WENO reconstructions
subroutine tracer_min_max_init(Reg, G, GV)
  type(tracer_registry_type), pointer    :: Reg  !< pointer to the tracer registry
  type(ocean_grid_type),      intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),    intent(in) :: GV   !< ocean vertical grid structure

  type(tracer_type), pointer :: Tr=>NULL()
  integer :: k, nz, m, ntr
  real :: Tmin, Tmax ! Global min and max of tracer

  nz = GV%ke ; ntr = Reg%ntr

  do m=1,ntr
    if (Reg%Tr(m)%advect_scheme > 2) then
      Tr => Reg%Tr(m)
      do k=1,nz
        call array_global_min_max(Tr%t(:,:,k), G, 1, Tmin, Tmax)
        Tr%Tmingg = Tmin
        Tr%Tmaxgg = Tmax
      enddo
    endif
  enddo

end subroutine tracer_min_max_init

!> \namespace mom_tracer_advect
!!
!!  This program contains the subroutines that advect tracers
!!  horizontally (i.e. along layers) using high-order WENO schemes (Balsara et al., 2016)
!!  using the Z-type smoothness indicators (Borges et al., 2008).
!!  We followed Suresh & Huynh (1997) and Balsara & Shu (2000 for the monotonicity preserving
!!  used along with the WENO schemes.
!!
!!  This scheme conserves the total amount of tracer while avoiding
!!  spurious maxima and minima of the tracer concentration

end module MOM_tracer_advect_weno
