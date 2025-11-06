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
  real, intent(in) :: q(3)   !< tracer concentration from cell i-1 to i+1
  real, intent(in) :: u      !< advection velocity
  real, intent(in) :: mu     !< cfl
  real, intent(out) :: wq    !< weno flux

  real :: wpl ! weno reconstruction on the cell interface

  wq = 0.0
  if (u > 0.0) then
    call weno3_reconstruction_interface(wpl, q(1), q(2), q(3))
  elseif (u < 0.0) then
    call weno3_reconstruction_interface(wpl, q(3), q(2), q(1))
  endif

  wq = u*wpl

end subroutine weno3_reconstruction

!> 3rd-order weno z-type reconstruction flux
pure subroutine weno3_reconstruction_interface(wpl, qm, q0, qp)
  real, intent(in) :: qm, q0, qp !< tracer concentration for 5-stencil wide
  real, intent(out) :: wpl       !< reconstruction value

  real :: b1, b2    ! smoothness indicator
  real :: d1, d2    ! linear weights
  real :: P1, P2    ! reconstructed polynomials
  real :: w1, w2    ! nonlinear weights
  real :: a1, a2
  real :: wnorm, tau

  ! linear weights
  d1 = 1.0/3.0 ; d2 = 2.0/3.0

  ! reconstructed polynomials
  P1 = 0.5*(-qm + 3.0*q0)
  P2 = 0.5*(q0 + qp)

  ! smoothness indicator
  b1 = (q0-qm)*(q0-qm)
  b2 = (qp-q0)*(qp-q0)

  ! Alpha values
  tau = abs(b2-b1)
  w1 = d1*weight_fac(tau, b1)
  w2 = d2*weight_fac(tau, b2)

  ! Normalization
  wnorm = 1.0/(w1 + w2)
  w1 = w1*wnorm ; w2 = w2*wnorm
  wpl = w1*P1 + w2*P2
  !if ((qp-q0)*(q0-qm) <= 0.) wpl = q0

end subroutine weno3_reconstruction_interface

!> 5th-order weno reconstruction subroutine and limiter
pure subroutine weno5_reconstruction(wq, q, u, mu)
  real, intent(in) :: q(5)       !< tracer concentration from i-2 to  i+3 respectively
  real, intent(in) :: u          !< advection velocity
  real, intent(in) :: mu         !< cfl
  real, intent(out) :: wq        !< weno flux

  real :: wpl ! weno reconstruction on the cell interface

  wq = 0.0
  if ( u > 0.0) then
    call weno5z_reconstruction_interface(wpl, q(1), q(2), q(3), q(4), q(5), mu)
  elseif ( u < 0.0) then
    call weno5z_reconstruction_interface(wpl, q(5), q(4), q(3), q(2), q(1), mu)
  endif

  wq = u*wpl

end subroutine weno5_reconstruction

!> 5th-order weno z-type reconstruction flux
pure subroutine weno5z_reconstruction_interface(wpl, qmm, qm, q0, qp, qpp, mu)
  real, intent(in) :: qmm, qm, q0, qp, qpp !< tracer concentration for 5-stencil wide
  real, intent(in) :: mu                   !< cfl
  real, intent(out) :: wpl                 !< reconstruction value

  real :: P0, P1, P2         ! reconstructed polynomials
  real :: b0, b1, b2         ! smoothness indicator
  real :: w0, w1, w2         ! nonlinear weights
  real :: d0, d1, d2         ! linear weights
  real :: a0, a1, a2
  real :: eps,  wnorm, tau
  real, parameter :: C1_6 = 1.0/6.0  ! [nondim]
  real :: dm1, dd0, dd1, dm4p, dm4m, mm1, mm2
  real :: qul, qmd, qlc, qmin, qmax, alpha

  ! linear weights
  d0 = 1.0/10.0 ; d1 = 6.0/10.0 ; d2 = 3.0/10.0

  ! Compute flux at left side of i+1/2
  ! First stencil
  P0 = ((2.0*qmm - 7.0*qm) + 11.0*q0)*C1_6
  b0 = (qmm*(4.0*qmm - 19.0*qm) + 11.0*q0) + (qm*(25.0*qm - 31.0*q0) + 10.0*(q0*q0))

  ! Second stencil
  P1 = ((-qm + 5.0*q0) + 2.0*qp)*C1_6
  b1 = (qm*(4.0*qm - 13.0*q0) - qp) + (q0*(13.0*q0 - 13.0*qp) + 4.0*(qp*qp))

  ! Third stencil
  P2 = ((2.0*q0 + 5.0*qp) - qpp)*C1_6
  b2 = (q0*(10.0*q0 - 31.0*qp) + 11.0*qpp) + (qp*(25.0*qp - 19.0*qpp) + 4.0*(qpp*qpp))

  ! Alpha values
  tau = abs(b2-b0)
  w0 = d0*weight_fac(tau, b0)
  w1 = d1*weight_fac(tau, b1)
  w2 = d2*weight_fac(tau, b2)

  wnorm = 1.0/((w0 + w1) + w2)
  w0 = w0*wnorm ; w1 = w1*wnorm ; w2 = w2*wnorm
  wpl = w0*P0 + w1*P1 + w2*P2

  ! Apply monotonicity preserving limiter based on Suresh & Huynh (1997)
  alpha = (1.0-mu)/mu
  !alpha = 2.0

  qul = q0 + alpha*(q0-qm)

  dm1 = qmm - 2.0*qm + q0
  dd0 = qp  - 2.0*q0 + qm
  dd1 = qpp - 2.0*qp + q0

  mm1 = 4.0*dd0-dd1 ; mm2 = 4.0*dd1-dd0
  dm4m = minmod2(mm1, mm2)

  mm1 = 4.0*dm1-dd0 ; mm2 = 4.0*dd0-dm1
  dm4m = minmod2(mm1, mm2)

  qmd = 0.5*(q0 + qp) - 0.5*dm4p
  qlc = 0.5*(q0+qul) + 0.5*alpha*dm4m

  qmin = max(min(q0,qp,qmd),min(q0,qul,qlc))
  qmax = min(max(q0,qp,qmd),max(q0,qul,qlc))
  wpl = max( min(qmin,qmax), wpl) ; wpl = min( max(qmin,qmax), wpl)

end subroutine weno5z_reconstruction_interface

!> 5th-order weno reconstruction for non-uniform grid and limiter
pure subroutine weno5NM_reconstruction(wq, q, u, ds, mu)
  real, intent(in) :: q(5)    !< tracer concentration for 5-stencil wide
  real, intent(in) :: u       !< advection velocity
  real, intent(in) :: ds(5)   !< grid sizes
  real, intent(in) :: mu      !< cfl
  real, intent(out) :: wq     !< reconstruction value

  real :: wpl ! weno reconstruction on the cell interface

  wq = 0.0
  if (u > 0.0) then
    call weno5NM_reconstruction_interface(wpl, q(1), q(2), q(3), q(4), q(5), ds, mu)
  elseif (u < 0.0) then
    call weno5NM_reconstruction_interface(wpl, q(5), q(4), q(3), q(2), q(1), ds, mu)
  endif

  wq = u*wpl

end subroutine weno5NM_reconstruction

!> 5th-order weno z-type reconstruction flux for non-uniform grid
pure subroutine weno5NM_reconstruction_interface(wpl, qmm, qm, q0, qp, qpp, dx, mu)
  real, intent(in) :: qmm, qm, q0, qp, qpp !< tracer concentration for 5-stencil wide
  real, intent(in) :: dx(5)                !< grid sizes
  real, intent(in) :: mu                   !< cfl
  real, intent(out) :: wpl                 !< reconstruction value

  real :: P0, P1, P2         ! reconstructed polynomials
  real :: b0, b1, b2         ! smoothness indicator
  real :: w0, w1, w2         ! nonlinear weights
  real :: d0, d1, d2         ! linear weights
  real :: a0, a1, a2, a3
  real :: eps,  wnorm, tau
  real :: qh, qhh, qhp, qhpp
  real :: dm1, dd0, dd1, dm4p, dm4m, mm1, mm2
  real :: qul, qmd, qlc, qmin, qmax, alpha

  ! Gamma values in Weno reconstruction
  d0 = 1.0/10.0 ; d1 = 6.0/10.0 ; d2 = 3.0/10.0

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
  w0 = d0*weight_fac(tau, b0)
  w1 = d1*weight_fac(tau, b1)
  w2 = d2*weight_fac(tau, b2)

  wnorm = 1.0/((w0 + w1) + w2)
  w0 = w0*wnorm ; w1 = w1*wnorm ; w2 = w2*wnorm
  wpl = w0*P0 + w1*P1 + w2*P2

  ! Apply monotonicity preserving limiter based on Suresh & Huynh (1997)
  alpha = (1.0-mu)/mu
  !alpha = 2.0

  qul = q0 + alpha*(q0-qm)

  dm1 = qmm - 2.0*qm + q0
  dd0 = qp  - 2.0*q0 + qm
  dd1 = qpp - 2.0*qp + q0

  mm1 = 4.0*dd0-dd1 ; mm2 = 4.0*dd1-dd0
  dm4p = minmod2(mm1, mm2)

  mm1 = 4.0*dm1-dd0 ; mm2 = 4.0*dd0-dm1
  dm4m = minmod2(mm1, mm2)

  qmd = 0.5*(q0 + qp) - 0.5*dm4p
  qlc = 0.5*(q0+qul) + 0.5*alpha*dm4m

  qmin = max(min(q0,qp,qmd),min(q0,qul,qlc))
  qmax = min(max(q0,qp,qmd),max(q0,qul,qlc))
  wpl = max( min(qmin,qmax), wpl) ; wpl = min( max(qmin,qmax), wpl)

end subroutine weno5NM_reconstruction_interface

!> 7th-order weno reconstruction subroutine and limiter
pure subroutine weno7_reconstruction(wq, q, u, mu)
  real, intent(in) :: q(7)       !< tracer concentration for 7-stencil wide
  real, intent(in) :: u          !< advection velocity
  real, intent(in) :: mu         !< cfl
  real, intent(out) :: wq        !< weno flux

  real :: wpl ! weno reconstruction on the cell interface

  wq = 0.0
  if (u > 0.0) then
    call weno7z_reconstruction_interface(wpl, q(1), q(2), q(3), q(4), q(5), q(6), q(7), mu)
  elseif (u < 0.0) then
    call weno7z_reconstruction_interface(wpl, q(7), q(6), q(5), q(4), q(3), q(2), q(1), mu)
  endif

  wq = u*wpl

end subroutine weno7_reconstruction

!> 7th-order weno z-type reconstruction flux
pure subroutine weno7z_reconstruction_interface(wpl, qm3, qm2, qm1, q0, qp1, qp2, qp3, mu)
  real, intent(in) :: qm3, qm2, qm1, q0, qp1, qp2, qp3 !< tracer concentration for 7-stencil wide
  real, intent(in) :: mu                               !< cfl
  real, intent(out) :: wpl                             !< reconstruction value

  real :: P0, P1, P2, P3     ! reconstructed polynomials
  real :: b0, b1, b2, b3     ! smoothness indicator
  real :: w0, w1, w2, w3     ! nonlinear weights
  real :: d0, d1, d2, d3     ! nonlinear weights
  real :: a0, a1, a2, a3
  real :: eps, tau, wnorm
  real, parameter :: C1_12 = 1.0/12.0  ! [nondim]
  real :: dm1, dd0, dd1, dm4p, dm4m, mm1, mm2
  real :: qul, qmd, qlc, qmin, qmax, alpha

  d0 = 1.0/35.0 ;  d1 = 12.0/35.0 ; d2 = 18.0/35.0 ; d3 = 4.0/35.0

  ! 1st stencil
  P0 = (((-3.0*qm3 + 13.0*qm2) - 23.0*qm1) + 25.0*q0)*C1_12
  b0 = ((qm3*((0.547*qm3 - 3.882*qm2) + (4.642*qm1 - 1.854*q0))) + &
      (qm2*((7.043*qm2 - 17.246*qm1) + 7.042*q0))) + &
      ((qm1*(11.003*qm1 - 9.402*q0)) + 2.107*(q0*q0))

  ! 2nd stencil
  P1 = (((qm2 - 5.0*qm1) + 13.0*q0) + 3.0*qp1)*C1_12
  b1 = ((qm2*((0.267*qm2 - 1.642*qm1) + (1.602*q0 - 0.494*qp1))) + &
         (qm1*((2.843*qm1 - 5.966*q0) + 1.922*qp1))) + &
         ((q0*(3.443*q0 - 2.522*qp1)) + 0.547*(qp1*qp1))

  ! 3rd stencil
  P2 = (((-qm1 + 7.0*q0) + 7.0*qp1) - qp2)*C1_12
  b2 = ((qm1*((0.547*qm1 - 2.522*q0) + (1.922*qp1 - 0.494*qp2))) + &
         (q0*((3.443*q0 - 5.966*qp1) + 1.602*qp2))) + &
         ((qp1*(2.843*qp1 - 1.642*qp2)) + 0.267*(qp2*qp2))

  ! 4rd stencil
  P3 = (((3.0*q0 + 13.0*qp1) - 5.0*qp2) + qp3)*C1_12
  b3 = ((q0*((2.107*q0 - 9.402*qp1) + (7.042*qp2 - 1.854*qp3))) + &
         (qp1*((11.003*qp1 - 17.246*qp2) + 4.642*qp3))) + &
         ((qp2*(7.043*qp2 - 3.882*qp3)) + 0.547*(qp3*qp3))

  ! Alpha values
  !tau = abs(b3 - 3.0*b2 + 3.0*b1 - b0)
  tau = abs(b3-b0)
  w0 = d0*weight_fac(tau, b0)
  w1 = d1*weight_fac(tau, b1)
  w2 = d2*weight_fac(tau, b2)
  w3 = d3*weight_fac(tau, b3)

  ! Normalization
  wnorm = 1.0/((w0 + w1) + (w2 + w3))
  w0 = w0*wnorm ; w1 = w1*wnorm ; w2 = w2*wnorm ; w3 = w3*wnorm
  wpl = w0*P0 + w1*P1 + w2*P2 + w3*P3

  ! Apply monotonicity preserving limiter based on Suresh & Huynh (1997)
  alpha = (1.0-mu)/mu
  !alpha = 2.0

  qul = q0 + alpha*(q0-qm1)

  dm1 = qm2 - 2.0*qm1 + q0
  dd0 = qp1  - 2.0*q0 + qm1
  dd1 = qp2 - 2.0*qp1 + q0

  mm1 = 4.0*dd0-dd1 ; mm2 = 4.0*dd1-dd0
  dm4m = minmod2(mm1, mm2)

  mm1 = 4.0*dm1-dd0 ; mm2 = 4.0*dd0-dm1
  dm4m = minmod2(mm1, mm2)

  qmd = 0.5*(q0 + qp1) - 0.5*dm4p
  qlc = 0.5*(q0+qul) + 0.5*alpha*dm4m

  qmin = max(min(q0,qp1,qmd),min(q0,qul,qlc))
  qmax = min(max(q0,qp1,qmd),max(q0,qul,qlc))
  wpl = max( min(qmin,qmax), wpl) ; wpl = min( max(qmin,qmax), wpl)

end subroutine weno7z_reconstruction_interface

!> 9th-order weno reconstruction subroutine and limiter
pure subroutine weno9_reconstruction(wq, q, u, mu)
  real, intent(in) :: q(9)       !< tracer concentrationi for 9-stencil wide
  real, intent(in) :: u          !< advection velocity
  real, intent(in) :: mu         !< cfl
  real, intent(out) :: wq        !< weno flux

  real :: wpl ! weno reconstruction on the cell interface

  wq = 0.0
  if (u > 0.0) then
    call weno9_reconstruction_interface(wpl, q(1), q(2), q(3), q(4), &
          q(5), q(6), q(7), q(8), q(9), mu)
  else
    call weno9_reconstruction_interface(wpl, q(9), q(8), q(7), q(6), &
          q(5), q(4), q(3), q(2), q(1), mu)
  endif

  wq = u*wpl

end subroutine weno9_reconstruction

!> 9th-order weno z-type reconstruction flux
pure subroutine weno9_reconstruction_interface(wpl, qm4, qm3, qm2, qm1, &
                q0, qp1, qp2, qp3, qp4, mu)
  real, intent(in) :: qm4, qm3, qm2, qm1, q0, qp1, qp2, qp3, qp4 !< tracer concentration
                                                                 !! for 9-stencil wide
  real, intent(in) :: mu                                    !< cfl
  real, intent(out) :: wpl                                  !< reconstruction value

  real :: b0, b1, b2, b3, b4               ! smoothness indicator
  real :: d0, d1, d2, d3, d4               ! linear weights
  real :: w0, w1, w2, w3, w4               ! nonlinear weights
  real :: P0, P1, P2, P3, P4               ! reconstructed polynomials
  real :: a0, a1, a2, a3, a4
  real :: eps,wnorm, tau
  real :: dm1, dd0, dd1, dm4p, dm4m, mm1, mm2
  real :: qul, qmd, qlc, qmin, qmax, alpha

  d0 = 1.0/126.0 ; d1 = 10.0/63.0 ; d2 = 10.0/21.0 ; d3 = 20.0/63.0 ; d4 = 5.0/126.0

  ! Compute flux at the right side of i+1/2
  call weno9_poly(P0, P1, P2, P3, P4, b0, b1, b2, b3, b4, &
          qm4, qm3, qm2, qm1, q0, qp1, qp2, qp3, qp4)

  ! Alpha values
  tau = abs(b0 - b4)
  w0 = d0*weight_fac(tau, b0)
  w1 = d1*weight_fac(tau, b1)
  w2 = d2*weight_fac(tau, b2)
  w3 = d3*weight_fac(tau, b3)
  w4 = d4*weight_fac(tau, b4)

  ! Normalization
  wnorm = 1.0/(((w0 + w1) + (w2 + w3)) + w4)
  w0 = w0*wnorm ; w1 = w1*wnorm ; w2 = w2*wnorm
  w3 = w3*wnorm ; w4 = w4*wnorm
  wpl = w0*P0 + w1*P1 + w2*P2 + w3*P3 + w4*P4

  ! Apply monotonicity preserving limiter based on Suresh & Huynh (1997)
  alpha = (1.0-mu)/mu
  !alpha = 2.0

  qul = q0 + alpha*(q0-qm1)

  dm1 = qm2 - 2.0*qm1 + q0
  dd0 = qp1  - 2.0*q0 + qm1
  dd1 = qp2 - 2.0*qp1 + q0

  mm1 = 4.0*dd0-dd1 ; mm2 = 4.0*dd1-dd0
  dm4p = minmod2(mm1, mm2)

  mm1 = 4.0*dm1-dd0 ; mm2 = 4.0*dd0-dm1
  dm4m = minmod2(mm1, mm2)

  qmd = 0.5*(q0 + qp1) - 0.5*dm4p
  qlc = 0.5*(q0+qul) + 0.5*alpha*dm4m

  qmin = max(min(q0,qp1,qmd),min(q0,qul,qlc))
  qmax = min(max(q0,qp1,qmd),max(q0,qul,qlc))
  wpl = max( min(qmin,qmax), wpl) ; wpl = min( max(qmin,qmax), wpl)

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
    wq_ppm = u*(aR - 0.5*mu*((aR-aL) - a6*(1.0 - 2.0/3.0 * mu)))
  else
    wq_ppm = u*(aL + 0.5*mu*((aR-aL) + a6*(1.0 - 2.0/3.0 * mu)))
  endif

end subroutine PPM_reconstruction

!> Compute the factor for the WENO weights
pure function weight_fac(tau, b) result(fac)
  real, intent(in)  :: tau  !< Difference of the smoothness indicator [A ~> a]
  real, intent(in)  :: b    !< The smoothness indicator [A ~> a]
  real :: fac               !< The factor for the weight [nondim]

  fac = 1.0e40; if (abs(b) > 1.0e-20*tau) fac = (1.0 + tau / b)**2

end function weight_fac

pure elemental function minmod4(a,b,c,d) result(r)
  real, intent(in) :: a, b, c, d
  real :: r

  if ( (a*b <= 0.0) .or. (a*c <= 0.0) .or. (a*d <= 0.0) ) then
    r = 0.0
  else
    r = sign( min( min(abs(a),abs(b)), min(abs(c),abs(d)) ), a )
  end if
end function minmod4

pure elemental function minmod2(a,b) result(r)
  real, intent(in) :: a, b
  real :: r

  ! 0 if opposite sign or either is zero; otherwise sign(a)*min(|a|,|b|)
  if (a*b <= 0.0) then
    r = 0.0
  else
    r = sign(min(abs(a),abs(b)), a)
  end if
end function minmod2

subroutine tracer_min_max_init(Reg, G, GV, local_advect_scheme)
  type(tracer_registry_type), pointer    :: Reg  !< pointer to the tracer registry
  type(ocean_grid_type),      intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),    intent(in) :: GV   !< ocean vertical grid structure
  integer :: local_advect_scheme(Reg%ntr)    !< contains the list of the advection for each tracer

  type(tracer_type), pointer :: Tr=>NULL()
  integer :: k, nz, m, ntr
  real :: Tmin, Tmax ! Global min and max of tracer
  real, dimension(SZK_(GV)) :: Tr_min, Tr_max


  nz = GV%ke ; ntr = Reg%ntr

  do m=1,ntr
    if (local_advect_scheme(m) > 2) then
      Tr => Reg%Tr(m)
      do k=1,nz
        call array_global_min_max(Tr%t(:,:,k), G, 1, Tmin, Tmax)
        Tr_min(k) = Tmin
        Tr_max(k) = Tmax
      enddo

      Tr%Tmingg = minval(Tr_min)
      Tr%Tmaxgg = maxval(Tr_max)

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
