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

contains

!> 3rd weno reconstruction subroutine and limiter
pure subroutine weno3_reconstruction(wq,  qm2, qm, q0, qp, qp2, qp3, u, qmin, qmax)

   real, intent(in) :: qm2, qm, q0, qp, qp2, qp3 !< tracer concentration from cell i-2 to i+3
                                                 !! respectively
   real, intent(in) :: u                         !< advection velocity
   real, intent(in) :: qmin, qmax                !< global min and max of tracer concentration
   real, intent(out) :: wq                       !< weno reconstruction at the cell
                                                 !! interface i+1/2

   real :: wmr ! wmr : weno reconstruction on the cell interface i-1/2
   real :: wpl ! wpl : weno reconstruction on the cell interface i+1/2
   real :: w0  ! w0 : 1st weight of N Gauss-Legendre quadrature weights over [-1/2,1/2], 
               ! N is weno reconstruction order

   w0 = 1.0
   wq = 0.0

   if (u > 0.0) then 
      call weno3_reconstruction_interface(wmr, qp2, qp, q0, qm, qm2)
      call weno3_reconstruction_interface(wpl, qm2, qm, q0, qp, qp2)
      ! maximum-principle limiter
      call PP_limiter(wq, q0, wmr, wpl, w0, qmin, qmax)
   elseif (u < 0.0) then
      call weno3_reconstruction_interface(wpl, qp3, qp2, qp, q0, qm)
      call weno3_reconstruction_interface(wmr, qm, q0, qp, qp2, qp3)
      ! maximum-principle limiter
      call PP_limiter(wq, qp, wmr, wpl, w0, qmin, qmax)
   endif

end subroutine weno3_reconstruction

!> 3rd-order weno reconstruction flux
pure subroutine weno3_reconstruction_interface(wq, qmm, qm, q0, qp, qpp)

   real, intent(in) :: qmm, qm, q0, qp, qpp !< tracer concentration for 5-stencil wide
   real, intent(out) :: wq                  !< reconstructed weno flux

   real :: P1, P2, P3                       ! reconstructed polynomials
   real :: b1, b2, b3                       ! smoothness indicator
   real :: w1, w2, w3, sw                   ! nonlinear weights
   real :: d1, d2, d3, ds                   ! linear weights
   real :: a1, a2, L1, L2                   ! Polynomial coefficients
   real :: tau, s0, s1, O01, O11, O0, O1
   real :: eps, k0, k1
   integer, parameter :: r = 2
   real :: dm1, dd0, dd1, dm4p, dm4m, mm1
   real :: mm2, qul, qmd, qlc, qmin, qmax, md

   eps = 1.0e-4

   ds = 16.0
   d1 = 10.0/ds ; d2 = 5.0/ds ; d3 = 1.0/ds

   L1 = 1.0/2.0; L2 = 1.0/6.0
   
   a1 = 0.5*(qp-qm)
   a2 = 0.5*(qm-2.0*q0+qp)

   ! Smoothness indicator

   b1 = a1**2 + (13.0/3.0)*(a2**2)
   b2 = a1**2

   k0 = (q0-qm)**2 ; k1 = (qp-q0)**2
   O01 = 10.0
   if (k0 >= k1) O01 = 1.0
   O11 = 11.0-O01

   O0 = O01/(O01+O11) ; O1 = 1.0-O0
   s0 = O0*(1.0 + (abs(k0-k1)**r/(k0+eps))) ; s1 = O1*(1.0 + (abs(k0-k1)**r/(k1+eps)))
   b3 = ((s0*(q0-qm)+s1*(qp-q0))**2)/(s0+s1)**2

   tau = ((abs(b1-b2)+abs(b1-b3))/2.0)**r

   ! Nonlinear weights
   w1 = d1*(1.0 + tau/(b1+eps))
   w2 = d2*(1.0 + tau/(b2+eps))
   w3 = d3*(1.0 + tau/(b3+eps))

   sw = w1+w2+w3
   w1 = w1/sw ; w2 = w2/sw ; w3 = w3/sw

   P1 = q0 + a1*L1 + a2*L2 ! Quadratic polynomial
   P2 = q0 + a1*L1
   P3 = q0

   ! Final reconstruction polynomial
   P1 = P1/d1 - d2*P2/d1 - d3*P3/d1
   wq = w1*P1 + w2*P2 + w3*P3
        
   ! Apply the monotonicity-preserving 

   dm1 = qmm - 2.0*qm + q0
   dd0 = qp - 2.0*q0 + qm
   dd1 = qpp - 2.0*qp + q0

   mm1 = 0.5*(sign(1.0,4.0*dd0-dd1) + sign(1.0,4.0*dd1-dd0))*min(abs(4.0*dd0-dd1),abs(4.0*dd1-dd0))
   mm2 = 0.5*(sign(1.0,dd0) + sign(1.0,dd1))*min(abs(dd0),abs(dd1))
   dm4p = 0.5*(sign(1.0,mm1) + sign(1.0,mm2))*min(abs(mm1),abs(mm2))

   mm1 = 0.5*(sign(1.0,4.0*dm1-dd0) + sign(1.0,4.0*dd0-dm1))*min(abs(4.0*dm1-dd0),abs(4.0*dd0-dm1))
   mm2 = 0.5*(sign(1.0,dm1) + sign(1.0,dd0))*min(abs(dm1),abs(dd0))
   dm4m = 0.5*(sign(1.0,mm1) + sign(1.0,mm2))*min(abs(mm1),abs(mm2))

   qul = q0 + 2.0*(q0-qm)
   qmd = 0.5*(q0 + qp) - 0.5*dm4p
   qlc = 0.5*(3.0*q0-qm) + (4.0/3.0)*dm4m
   !qlc = 0.5*(3.0*q0+qm) - (4.0/3.0)*dm4m

   qmin = max(min(q0,qp,qmd),min(q0,qul,qlc))
   qmax = min(max(q0,qp,qmd),max(q0,qul,qlc))

   md = 0.5*(sign(1.0,qmin-wq) + sign(1.0,qmax-wq))*min(abs(qmin-wq),abs(qmax-wq))
   wq = wq + md

end subroutine weno3_reconstruction_interface

!> 3rd-order weno z-type reconstruction flux
subroutine weno3z_reconstruction_interface(wq, qmm, qm, q0, qp, qpp)

   real, intent(in) :: qmm, qm, q0, qp, qpp !< tracer concentration for 5-stencil wide
   real, intent(out) :: wq                  !< reconstructed weno flux

   real :: b1, b2    ! smoothness indicator
   real :: d1, d2    ! linear weights
   real :: P1, P2    ! reconstructed polynomials
   real :: w1, w2    ! nonlinear weights
   real :: a1, a2
   real :: eps, wnorm, tau
   real :: dm1, dd0, dd1, dm4p, dm4m, mm1, mm2
   real :: qul, qmd, qlc, qmin, qmax, md

   ! reconstructed polynomials
   P1 = 0.5*(-qm + 3.0*q0)
   P2 = 0.5*(q0 + qp)

   ! smoothness indicator
   b1 = (q0-qm)*(q0-qm)
   b2 = (qp-q0)*(qp-q0)

   d1 = 1.0/3.0 ; d2 = 2.0/3.0

   ! Alpha values
   eps = 1.0e-6
   tau = abs(b2-b1)
   a1 = d1*(1.0 + (tau/(b1+eps))**2)
   a2 = d2*(1.0 + (tau/(b2+eps))**2)

   ! Normalization
   wnorm = a1+a2
   w1 = a1 / wnorm
   w2 = a2 / wnorm

   wq = w1*P1 + w2*P2

   ! Monotonicity Preserving
   dm1 = qmm - 2.0*qm + q0
   dd0 = qp - 2.0*q0 + qm
   dd1 = qpp - 2.0*qp + q0

   mm1 = 0.5*(sign(1.0,4.0*dd0-dd1) + sign(1.0,4.0*dd1-dd0))*min(abs(4.0*dd0-dd1),abs(4.0*dd1-dd0))
   mm2 = 0.5*(sign(1.0,dd0) + sign(1.0,dd1))*min(abs(dd0),abs(dd1))
   dm4p = 0.5*(sign(1.0,mm1) + sign(1.0,mm2))*min(abs(mm1),abs(mm2))

   mm1 = 0.5*(sign(1.0,4.0*dm1-dd0) + sign(1.0,4.0*dd0-dm1))*min(abs(4.0*dm1-dd0),abs(4.0*dd0-dm1))
   mm2 = 0.5*(sign(1.0,dm1) + sign(1.0,dd0))*min(abs(dm1),abs(dd0))
   dm4m = 0.5*(sign(1.0,mm1) + sign(1.0,mm2))*min(abs(mm1),abs(mm2))

   qul = q0 + 2.0*(q0-qm)
   qmd = 0.5*(q0 + qp) - 0.5*dm4p
   qlc = 0.5*(3.0*q0-qm) + (4.0/3.0)*dm4m
   !qlc = 0.5*(3.0*q0+qm) - (4.0/3.0)*dm4m

   qmin = max(min(q0,qp,qmd),min(q0,qul,qlc))
   qmax = min(max(q0,qp,qmd),max(q0,qul,qlc))

   md = 0.5*(sign(1.0,qmin-wq) + sign(1.0,qmax-wq))*min(abs(qmin-wq),abs(qmax-wq))
   wq = wq + md

end subroutine weno3z_reconstruction_interface

!> 5th-order weno reconstruction subroutine and limiter
pure subroutine weno5_reconstruction(wq, qm2, qm, q0, qp, qp2, qp3, u, qmin, qmax)

   real, intent(in) :: qm2, qm, q0, qp, qp2, qp3 !< tracer concentration from i-2 to  i+3 respectively
   real, intent(in) :: u                         !< advection velocity
   real, intent(in) :: qmin, qmax                !< global min and max of tracer concentration
   real, intent(out) :: wq                       !< weno reconstruction at the interface i+1/2

   real :: wmr ! wmr : weno reconstruction on the cell interface i-1/2
   real :: wpl ! wpl : weno reconstruction on the cell interface i+1/2
   real :: w0  ! w0 : 1st weight of N Gauss-Legendre quadrature weights over [-1/2,1/2],
               ! N is weno reconstruction order

   w0 = 5.0/18.0
   wq = 0.0

   if (u > 0.0) then
      call weno5_reconstruction_interface(wmr, qp2, qp, q0, qm, qm2)  ! i-1/2
      call weno5_reconstruction_interface(wpl, qm2, qm, q0, qp, qp2)  ! i+1/2
      ! maximum-principle limiter
      call PP_limiter(wq, q0, wmr, wpl, w0, qmin, qmax)
   elseif (u < 0.0) then
      call weno5_reconstruction_interface(wpl, qp3, qp2, qp, q0, qm)
      call weno5_reconstruction_interface(wmr, qm, q0, qp, qp2, qp3)
      ! maximum-principle limiter
      call PP_limiter(wq, qp, wmr, wpl, w0, qmin, qmax)
   endif

end subroutine weno5_reconstruction

!> 5th-order weno reconstruction flux
pure subroutine weno5_reconstruction_interface(wq, qmm, qm, q0, qp, qpp)

   real, intent(in) :: qmm, qm, q0, qp, qpp !< tracer concentration for 5-stencil wide
   real, intent(out) :: wq                  !< reconstructed weno flux

   real :: P1, P2, P3, P4, P5               ! reconstructed polynomials
   real :: b1, b2, b3, b4, b5               ! smoothness indicator
   real :: w1, w2, w3, w4, w5, sw           ! nonlinear weights
   real :: d1, d2, d3, d4, d5, ds           ! linear weights

   real :: a1, a2, a3, a4, L1, L2, L3, L4   ! Polynomial coefficients
   real :: tau, s0, s1, O01, O11, O0, O1
   real :: eps, k0, k1
   integer, parameter :: r = 2
   real :: dm1, dd0, dd1, dm4p, dm4m, mm1, mm2
   real :: qul, qmd, qlc, qmin, qmax, md

   eps = 1.0e-4

   ds = 141.0
   d1 = 100.0/ds; d2 = 25.0/ds; d3 = 10.0/ds ; d4 = 5.0/ds ; d5 = 1.0/ds

   L1 = 1.0/2.0; L2 = 1.0/6.0; L3 = 1.0/20.0 ; L4 = 1.0/70.0

   a1 = (11.0*qmm - 82.0*qm + 82.0*qp - 11.0*qpp)/120.0
   a2 = (-3.0*qmm + 40.0*qm - 74.0*q0 + 40.0*qp - 3.0*qpp)/56.0
   a3 = (-qmm + 2.0*qm - 2.0*qp + qpp)/12.0
   a4 = (qmm - 4.0*qm + 6.0*q0 - 4.0*qp + qpp)/24.0

   b1 = (a1+(a3/10.0))**2 + (13.0/3.0)*(a2+(123.0*a4/455.0))**2 &
           + (781.0/20.0)*(a3**2) + (1421461.0/2275.0)*(a4**2)

   b2 = (a1+(a3/10.0))**2 + (13.0/3.0)*(a2**2) + (781.0/20.0)*(a3**2)
   b3 = a1**2 + (13.0/3.0)*(a2**2) 
   b4 = a1**2

   k0 = (q0-qm)**2 ; k1 = (qp-q0)**2
   O01 = 10.0
   if (k0 >= k1) O01 = 1.0
   O11 = 11.0-O01

   O0 = O01/(O01+O11) ; O1 = 1.0-O0
   s0 = O0*(1.0 + (abs(k0-k1)**r/(k0+eps))) ; s1 = O1*(1.0 + (abs(k0-k1)**r/(k1+eps)))
   b5 = ((s0*(q0-qm)+s1*(qp-q0))**2)/(s0+s1)**2

   tau = ((abs(b1-b2)+abs(b1-b3)+abs(b1-b4)+abs(b1-b5))/4.0)**r

   w1 = d1*(1.0 + tau/(b1+eps))
   w2 = d2*(1.0 + tau/(b2+eps))
   w3 = d3*(1.0 + tau/(b3+eps))
   w4 = d4*(1.0 + tau/(b4+eps))
   w5 = d5*(1.0 + tau/(b5+eps))

   sw = w1+w2+w3+w4+w5
   w1 = w1/sw ; w2 = w2/sw ; w3 = w3/sw ; w4 = w4/sw ; w5 = w5/sw

   P1 = q0 + a1*L1 + a2*L2 + a3*L3 + a4*L4 
   P2 = q0 + a1*L1 + a2*L2 + a3*L3
   P3 = q0 + a1*L1 + a2*L2 
   P4 = q0 + a1*L1 
   P5 = q0 

   P1 = P1/d1 - d2*P2/d1 - d3*P3/d1 - d4*P4/d1 - d5*P5/d1

   wq = w1*P1 + w2*P2 + w3*P3 + w4*P4 + w5*P5

   ! Apply the monotonicity-preserving

   dm1 = qmm - 2.0*qm + q0
   dd0 = qp  - 2.0*q0 + qm
   dd1 = qpp - 2.0*qp + q0

   mm1 = 0.5*(sign(1.0,4.0*dd0-dd1) + sign(1.0,4.0*dd1-dd0))*min(abs(4.0*dd0-dd1),abs(4.0*dd1-dd0))
   mm2 = 0.5*(sign(1.0,dd0) + sign(1.0,dd1))*min(abs(dd0),abs(dd1))
   dm4p = 0.5*(sign(1.0,mm1) + sign(1.0,mm2))*min(abs(mm1),abs(mm2))

   mm1 = 0.5*(sign(1.0,4.0*dm1-dd0) + sign(1.0,4.0*dd0-dm1))*min(abs(4.0*dm1-dd0),abs(4.0*dd0-dm1))
   mm2 = 0.5*(sign(1.0,dm1) + sign(1.0,dd0))*min(abs(dm1),abs(dd0))
   dm4m = 0.5*(sign(1.0,mm1) + sign(1.0,mm2))*min(abs(mm1),abs(mm2))

   qul = q0 + 2.0*(q0-qm)
   qmd = 0.5*(q0 + qp) - 0.5*dm4p
   qlc = 0.5*(3.0*q0-qm) + (4.0/3.0)*dm4m
   !qlc = 0.5*(3.0*q0+qm) - (4.0/3.0)*dm4m

   qmin = max(min(q0,qp,qmd),min(q0,qul,qlc))
   qmax = min(max(q0,qp,qmd),max(q0,qul,qlc))

   md = 0.5*(sign(1.0,qmin-wq) + sign(1.0,qmax-wq))*min(abs(qmin-wq),abs(qmax-wq))
   wq = wq + md

end subroutine weno5_reconstruction_interface

!> 5th-order weno z-type reconstruction flux
subroutine weno5z_reconstruction_interface(wq, qmm, qm, q0, qp, qpp)

   real, intent(in) :: qmm, qm, q0, qp, qpp !< tracer concentration for 5-stencil wide
   real, intent(out) :: wq                  !< reconstructed weno flux

   real :: P0, P1, P2         ! reconstructed polynomials
   real :: b0, b1, b2         ! smoothness indicator
   real :: w0, w1, w2         ! nonlinear weights
   real :: d0, d1, d2         ! linear weights
   real :: a0, a1, a2
   real :: eps,  wnorm, tau
   integer, parameter :: r = 2
   real :: dm1, dd0, dd1, dm4p, dm4m, mm1, mm2
   real :: qul, qmd, qlc, qmin, qmax, md

   ! Gamma values in Weno reconstruction
   d0 = 1.0/10.0
   d1 = 6.0/10.0
   d2 = 3.0/10.0

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
   eps = 1.0e-6
   tau = abs(b2-b0)
   a0 = d0*(1.0 + (tau/(b0+eps))**r)
   a1 = d1*(1.0 + (tau/(b1+eps))**r)
   a2 = d2*(1.0 + (tau/(b2+eps))**r)

   wnorm = a0+a1+a2
   w0 = a0/wnorm
   w1 = a1/wnorm
   w2 = a2/wnorm

   wq = w0*P0 + w1*P1 + w2*P2

   ! Monotonicity Preserving
   dm1 = qmm - 2.0*qm + q0
   dd0 = qp  - 2.0*q0 + qm
   dd1 = qpp - 2.0*qp + q0

   mm1 = 0.5*(sign(1.0,4.0*dd0-dd1) + sign(1.0,4.0*dd1-dd0))*min(abs(4.0*dd0-dd1),abs(4.0*dd1-dd0))
   mm2 = 0.5*(sign(1.0,dd0) + sign(1.0,dd1))*min(abs(dd0),abs(dd1))
   dm4p = 0.5*(sign(1.0,mm1) + sign(1.0,mm2))*min(abs(mm1),abs(mm2))

   mm1 = 0.5*(sign(1.0,4.0*dm1-dd0) + sign(1.0,4.0*dd0-dm1))*min(abs(4.0*dm1-dd0),abs(4.0*dd0-dm1))
   mm2 = 0.5*(sign(1.0,dm1) + sign(1.0,dd0))*min(abs(dm1),abs(dd0))
   dm4m = 0.5*(sign(1.0,mm1) + sign(1.0,mm2))*min(abs(mm1),abs(mm2))

   qul = q0 + 2.0*(q0-qm)
   qmd = 0.5*(q0 + qp) - 0.5*dm4p
   qlc = 0.5*(3.0*q0-qm) + (4.0/3.0)*dm4m
   !qlc = 0.5*(3.0*q0+qm) - (4.0/3.0)*dm4m

   qmin = max(min(q0,qp,qmd),min(q0,qul,qlc))
   qmax = min(max(q0,qp,qmd),max(q0,qul,qlc))

   md = 0.5*(sign(1.0,qmin-wq) + sign(1.0,qmax-wq))*min(abs(qmin-wq),abs(qmax-wq))
   wq = wq + md

end subroutine weno5z_reconstruction_interface

!> 7th-order weno reconstruction subroutine and limiter
pure subroutine weno7_reconstruction(wq, qm3, qm2, qm1, q0, qp1, qp2, qp3, qp4, u, qmin, qmax) 

   real, intent(in) :: qm3, qm2, qm1, q0, qp1, qp2, qp3, qp4 !< tracer concentration 
                                                             !! from i-3 to i+4 respectively
   real, intent(inout) :: u                         !< advection velocity
   real, intent(in) :: qmin, qmax                !< global min and max of tracer concentration
   real, intent(out) :: wq                       !< weno reconstruction at the interface i+1/2

   real :: wmr ! wmr : weno reconstruction on the cell interface i-1/2
   real :: wpl ! wpl : weno reconstruction on the cell interface i+1/2
   real :: w0  ! w0 : 1st weight of N Gauss-Legendre quadrature weights over [-1/2,1/2],
               ! N is weno reconstruction order

   w0 = 5.0/18.0 !(322.0-13.0*sqrt(70.0))/1800.0
   wq = 0.0
   
   if (u > 0.0) then
      call weno7_reconstruction_interface(wmr, qp3, qp2, qp1, q0, qm1, qm2, qm3)
      call weno7_reconstruction_interface(wpl, qm3, qm2, qm1, q0, qp1, qp2, qp3)
      ! maximum-principle limiter
      call PP_limiter(wq, q0, wmr, wpl, w0, qmin, qmax)
   elseif (u < 0.0) then
      call weno7_reconstruction_interface(wpl, qp4, qp3, qp2, qp1, q0, qm1, qm2)
      call weno7_reconstruction_interface(wmr, qm2, qm1, q0, qp1, qp2, qp3, qp4)
      ! maximum-principle limiter
      call PP_limiter(wq, qp1, wmr, wpl, w0, qmin, qmax)
   endif

end subroutine weno7_reconstruction

!> 7th-order weno reconstruction flux
pure subroutine weno7_reconstruction_interface(wq, qm3, qmm, qm, q0, qp, qpp, qp3)

   real, intent(in) :: qm3, qmm, qm, q0, qp, qpp, qp3 !< tracer concentration for 7-stencil wide
   real, intent(out) :: wq                            !< reconstructed weno flux

   real :: P1, P2, P3, P4, P5, P6, P7       ! reconstructed polynomials
   real :: b1, b2, b3, b4, b5, b6, b7       ! smoothness indicator
   real :: w1, w2, w3, w4, w5, w6, w7, sw   ! nonlinear weights
   real :: d1, d2, d3, d4, d5, d6, d7, ds   ! linear weights

   real :: a1, a2, a3, a4, a5, a6           ! Polynomial coefficients
   real :: L1, L2, L3, L4, L5, L6
   real :: tau, s0, s1, O01, O11, O0, O1
   real :: eps, k0, k1
   integer, parameter :: r = 3
   real :: dm1, dd0, dd1, dm4p, dm4m, mm1, mm2
   real :: qul, qmd, qlc, qmin, qmax, md, tau_tmp

   eps = 1.0e-4

   ds = 1266.0
   d1 = 1000.0/ds; d2 = 125.0/ds ; d3 = 100.0/ds; d4 = 25.0/ds ; d5 = 10.0/ds
   d6 = 5.0/ds ; d7 = 1.0/ds

   L1 = 1.0/2.0; L2 = 1.0/6.0; L3 = 1.0/20.0 ; L4 = 1.0/70.0
   L5 = 1.0/252.0; L6 = 1.0/924.0
        
   a1 = (-7843.0*qm + 1688.0*qmm - 191.0*qm3 + 7843.0*qp - 1688.0*qpp + 191.0*qp3)/10080.0
   a2 = (8385.0*qm - 1014.0*qmm + 79.0*qm3 - 14900.0*q0 + 8385.0*qp - 1014.0*qpp + 79.0*qp3)/10080.0
   a3 = (61.0*qm - 38.0*qmm + 5.0*qm3 - 61.0*qp + 38.0*qpp - 5.0*qp3)/216.0
   a4 = (-459.0*qm + 144.0*qmm - 13.0*qm3 + 656.0*q0 - 459.0*qp + 144.0*qpp - 13.0*qp3)/1584.0
   a5 = (-5.0*qm + 4.0*qmm - qm3 + 5.0*qp - 4.0*qpp + qp3)/240.0
   a6 = (15.0*qm - 6.0*qmm + qm3 - 20.0*q0 + 15.0*qp - 6.0*qpp + qp3)/720.0

   b1 = (a1+(a3/10.0)+(a5/126.0))**2 + (13.0/3.0)*(a2+(123.0*a4/455.0)+(85.0*a6/2002.0))**2 &
           + (781.0/20.0)*(a3+(26045.0*a5/49203.0))**2 + (1421461.0/2275.0)*(a4+(81596225.0*a6/93816426.0))**2 &
           + (21520059541.0/1377684.0)*(a5**2) + (15510384942580921.0/27582029244.0)*(a6**2)

   b2 = (a1+(a3/10.0)+(a5/126.0))**2 + (13.0/3.0)*(a2+(123.0*a4/455.0))**2 &
           + (781.0/20.0)*(a3+(26045.0*a5/49203.0))**2 + (1421461.0/2275.0)*(a4**2) &
           + (21520059541.0/1377684.0)*(a5**2)

   b3 = (a1+(a3/10.0))**2 + (13.0/3.0)*(a2+(123.0*a4/455.0))**2 &
           + (781.0/20.0)*a3**2 + (1421461.0/2275.0)*(a4**2) 

   b4 = (a1+(a3/10.0))**2 + (13.0/3.0)*(a2**2) + (781.0/20.0)*(a3**2)
   b5 = a1**2 + (13.0/3.0)*(a2**2)
   b6 = a1**2

   k0 = (q0-qm)**2 ; k1 = (qp-q0)**2
   O01 = 10.0
   if (k0 >= k1) O01 = 1.0
   O11 = 11.0-O01

   O0 = O01/(O01+O11) ; O1 = 1.0-O0
   s0 = O0*(1.0 + (abs(k0-k1)*abs(k0-k1)*abs(k0-k1)/(k0+eps)))
   s1 = O1*(1.0 + (abs(k0-k1)*abs(k0-k1)*abs(k0-k1)/(k1+eps)))
   b7 = ((s0*(q0-qm)+s1*(qp-q0))**2)/(s0+s1)**2
   
   tau_tmp = ((abs(b1-b2)+abs(b1-b3)+abs(b1-b4)+abs(b1-b5)+abs(b1-b6)+abs(b1-b7))/6.0)
   tau = tau_tmp*tau_tmp*tau_tmp

   w1 = d1*(1.0 + tau/(b1+eps))
   w2 = d2*(1.0 + tau/(b2+eps))
   w3 = d3*(1.0 + tau/(b3+eps))
   w4 = d4*(1.0 + tau/(b4+eps))
   w5 = d5*(1.0 + tau/(b5+eps))
   w6 = d6*(1.0 + tau/(b6+eps))
   w7 = d7*(1.0 + tau/(b7+eps))

   sw = w1+w2+w3+w4+w5+w6+w7
   w1 = w1/sw ; w2 = w2/sw ; w3 = w3/sw ; w4 = w4/sw
   w5 = w5/sw ; w6 = w6/sw ; w7 = w7/sw

   P1 = q0 + a1*L1 + a2*L2 + a3*L3 + a4*L4 + a5*L5 + a6*L6
   P2 = q0 + a1*L1 + a2*L2 + a3*L3 + a4*L4 + a5*L5
   P3 = q0 + a1*L1 + a2*L2 + a3*L3 + a4*L4
   P4 = q0 + a1*L1 + a2*L2 + a3*L3
   P5 = q0 + a1*L1 + a2*L2
   P6 = q0 + a1*L1
   P7 = q0

   P1 = P1/d1 - d2*P2/d1 - d3*P3/d1 - d4*P4/d1 - d5*P5/d1 - d6*P6/d1 - d7*P7/d1
   !P1 = (P1 - d2*P2 - d3*P3 - d4*P4 - d5*P5 - d6*P6 - d7*P7)/d1
   wq = w1*P1 + w2*P2 + w3*P3 + w4*P4 + w5*P5 + w6*P6 + w7*P7

   ! Apply the monotonicity-preserving
   dm1 = qmm - 2.0*qm + q0
   dd0 = qm  - 2.0*q0 + qp
   dd1 = q0  - 2.0*qp + qpp

   mm1 = 0.5*(sign(1.0,4.0*dd0-dd1) + sign(1.0,4.0*dd1-dd0))*min(abs(4.0*dd0-dd1),abs(4.0*dd1-dd0))
   mm2 = 0.5*(sign(1.0,dd0) + sign(1.0,dd1))*min(abs(dd0),abs(dd1))
   dm4p = 0.5*(sign(1.0,mm1) + sign(1.0,mm2))*min(abs(mm1),abs(mm2))

   mm1 = 0.5*(sign(1.0,4.0*dm1-dd0) + sign(1.0,4.0*dd0-dm1))*min(abs(4.0*dm1-dd0),abs(4.0*dd0-dm1))
   mm2 = 0.5*(sign(1.0,dm1) + sign(1.0,dd0))*min(abs(dm1),abs(dd0))
   dm4m = 0.5*(sign(1.0,mm1) + sign(1.0,mm2))*min(abs(mm1),abs(mm2))

   qul = q0 + 2.0*(q0-qm)
   qmd = 0.5*(q0 + qp) - 0.5*dm4p
   qlc = 0.5*(3.0*q0-qm) + (4.0/3.0)*dm4m
   !qlc = 0.5*(3.0*q0+qm) - (4.0/3.0)*dm4m

   qmin = max(min(q0,qp,qmd),min(q0,qul,qlc))
   qmax = min(max(q0,qp,qmd),max(q0,qul,qlc))

   md = 0.5*(sign(1.0,qmin-wq) + sign(1.0,qmax-wq))*min(abs(qmin-wq),abs(qmax-wq))
   wq = wq + md

end subroutine weno7_reconstruction_interface

!> 7th-order weno z-type reconstruction flux
subroutine weno7z_reconstruction_interface(wq, qm3, qm2, qm1, q0, qp1, qp2, qp3)

   real, intent(in) :: qm3, qm2, qm1, q0, qp1, qp2, qp3 !< tracer concentration for 7-stencil wide
   real, intent(out) :: wq                              !< reconstructed weno flux

   real :: P0, P1, P2, P3     ! reconstructed polynomials
   real :: b0, b1, b2, b3     ! smoothness indicator
   real :: w0, w1, w2, w3     ! nonlinear weights
   real :: d0, d1, d2, d3     ! nonlinear weights
   real :: a0, a1, a2, a3
   real :: eps, tau, wnorm
   integer, parameter :: r = 2
   real :: dm1, dd0, dd1, dm4p, dm4m, mm1, mm2
   real :: qul, qmd, qlc, qmin, qmax, md

   d0 = 1.0/35.0
   d1 = 12.0/35.0
   d2 = 18.0/35.0
   d3 = 4.0/35.0

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
   eps = 1.0e-6
   tau = abs(b3 + 3.0 * b2 - 3.0 * b1 - b0)
   a0 = d0*(1.0 + (tau/(b0+eps))**r)
   a1 = d1*(1.0 + (tau/(b1+eps))**r)
   a2 = d2*(1.0 + (tau/(b2+eps))**r)
   a3 = d3*(1.0 + (tau/(b3+eps))**r)

   ! Normalization
   wnorm = a0+a1+a2+a3
   w0 = a0/wnorm
   w1 = a1/wnorm
   w2 = a2/wnorm
   w3 = a3/wnorm

   wq = w0*P0 + w1*P1 + w2*P2 + w3*P3

   ! Monotonicity Preserving
   dm1 = qm2 - 2.0*qm1 + q0
   dd0 = qm1 - 2.0*q0  + qp1
   dd1 = q0 - 2.0*qp1 + qp2

   mm1 = 0.5*(sign(1.0,4.0*dd0-dd1) + sign(1.0,4.0*dd1-dd0))*min(abs(4.0*dd0-dd1),abs(4.0*dd1-dd0))
   mm2 = 0.5*(sign(1.0,dd0) + sign(1.0,dd1))*min(abs(dd0),abs(dd1))
   dm4p = 0.5*(sign(1.0,mm1) + sign(1.0,mm2))*min(abs(mm1),abs(mm2))

   mm1 = 0.5*(sign(1.0,4.0*dm1-dd0) + sign(1.0,4.0*dd0-dm1))*min(abs(4.0*dm1-dd0),abs(4.0*dd0-dm1))
   mm2 = 0.5*(sign(1.0,dm1) + sign(1.0,dd0))*min(abs(dm1),abs(dd0))
   dm4m = 0.5*(sign(1.0,mm1) + sign(1.0,mm2))*min(abs(mm1),abs(mm2))

   qul = q0 + 2.0*(q0-qm1)
   qmd = 0.5*(q0 + qp1) - 0.5*dm4p
   qlc = 0.5*(3.0*q0-qm1) + (4.0/3.0)*dm4m
   !qlc = 0.5*(3.0*q0+qm1) - (4.0/3.0)*dm4m

   qmin = max(min(q0,qp1,qmd),min(q0,qul,qlc))
   qmax = min(max(q0,qp1,qmd),max(q0,qul,qlc))

   md = 0.5*(sign(1.0,qmin-wq) + sign(1.0,qmax-wq))*min(abs(qmin-wq),abs(qmax-wq))
   wq = wq + md

end subroutine weno7z_reconstruction_interface

pure subroutine weno9_reconstruction(wq, qm4, qm3, qm2, qm1, q0, qp1, qp2, qp3, qp4, qp5, u, qmin, qmax)

   real, intent(in) :: qm4, qm3, qm2, qm1, q0, qp1, qp2, qp3, qp4, qp5 !< tracer concentration
                                                                       !! from i-4 to i+5
   real, intent(in) :: u                         !< advection velocity
   real, intent(in) :: qmin, qmax                !< global min and max of tracer concentration
   real, intent(out) :: wq                       !< weno reconstruction at the interface i+1/2

   real :: wmr ! wmr : weno reconstruction on the cell interface i-1/2
   real :: wpl ! wpl : weno reconstruction on the cell interface i+1/2
   real :: w0  ! w0 : 1st weight of N Gauss-Legendre quadrature weights over [-1/2,1/2], 
               ! N is weno reconstruction order

   w0 = 0.5*0.1294849661688697
   wq = 0.0

   if (u > 0.0) then
      call weno9_reconstruction_interface(wmr, qp4, qp3, qp2, qp1, q0, qm1, qm2, qm3, qm4)
      call weno9_reconstruction_interface(wpl, qm4, qm3, qm2, qm1, q0, qp1, qp2, qp3, qp4)
      ! maximum-principle limiter
      call PP_limiter(wq, q0, wmr, wpl, w0, qmin, qmax)
   elseif (u < 0.0) then
      call weno9_reconstruction_interface(wpl, qp5, qp4, qp3, qp2, qp1, q0, qm1, qm2, qm3)
      call weno9_reconstruction_interface(wmr, qm3, qm2, qm1, q0, qp1, qp2, qp3, qp4, qp5)
      ! maximum-principle limiter
      call PP_limiter(wq, qp1, wmr, wpl, w0, qmin, qmax)
   endif

end subroutine weno9_reconstruction

pure subroutine weno9_reconstruction_interface(wq, qm4, qm3, qm2, qm1, q0, qp1, qp2, qp3, qp4)

   real, intent(in) :: qm4, qm3, qm2, qm1, q0, qp1, qp2, qp3, qp4 !< tracer concentration 
                                                                  !! for 7-stencil wide
   real, intent(out) :: wq                                        !< reconstructed weno flux

   real :: b0, b1, b2, b3, b4               ! smoothness indicator
   real :: d0, d1, d2, d3, d4               ! linear weights
   real :: w0, w1, w2, w3, w4               ! nonlinear weights
   real :: P0, P1, P2, P3, P4               ! reconstructed polynomials
   real :: a0, a1, a2, a3, a4                
   real :: eps,wnorm, tau 
   integer, parameter :: r = 2
   real :: dm1, dd0, dd1, dm4p, dm4m, s1, s2
   real :: qul, qmd, qlc, qmin, qmax, md

   call weno9_poly(P0, P1, P2, P3, P4, b0, b1, b2, b3, b4, qm4, qm3, qm2, qm1, q0, qp1, qp2, qp3, qp4)

   d0 = 1.0/126.0
   d1 = 10.0/63.0
   d2 = 10.0/21.0
   d3 = 20.0/63.0
   d4 = 5.0/126.0

   ! Alpha values
   eps = 1.0e-6
   tau = abs(b0 - b4)
   a0 = d0*(1.0 + (tau/(b0+eps))**r)
   a1 = d1*(1.0 + (tau/(b1+eps))**r)
   a2 = d2*(1.0 + (tau/(b2+eps))**r)
   a3 = d3*(1.0 + (tau/(b3+eps))**r)
   a4 = d4*(1.0 + (tau/(b4+eps))**r)

   ! Normalization
   wnorm = 1.0/(a0+a1+a2+a3+a4)
   w0 = a0*wnorm
   w1 = a1*wnorm
   w2 = a2*wnorm
   w3 = a3*wnorm
   w4 = a4*wnorm

   wq = w0*P0 + w1*P1 + w2*P2 + w3*P3 + w4*P4

   ! Apply the monotonicity-preserving
   dm1 = qm2 - 2.0*qm1 + q0
   dd0 = qp1 - 2.0*q0 + qm1
   dd1 = qp2 - 2.0*qp1 + q0

   s1 = 0.5*(sign(1.0,4.0*dd0-dd1) + sign(1.0,4.0*dd1-dd0))*min(abs(4.0*dd0-dd1),abs(4.0*dd1-dd0))
   s2 = 0.5*(sign(1.0,dd0) + sign(1.0,dd1))*min(abs(dd0),abs(dd1))
   dm4p = 0.5*(sign(1.0,s1) + sign(1.0,s2))*min(abs(s1),abs(s2))

   s1 = 0.5*(sign(1.0,4.0*dm1-dd0) + sign(1.0,4.0*dd0-dm1))*min(abs(4.0*dm1-dd0),abs(4.0*dd0-dm1))
   s2 = 0.5*(sign(1.0,dm1) + sign(1.0,dd0))*min(abs(dm1),abs(dd0))
   dm4m = 0.5*(sign(1.0,s1) + sign(1.0,s2))*min(abs(s1),abs(s2))

   qul = q0 + 2.0*(q0-qm1)
   qmd = 0.5*(q0 + qp1) - 0.5*dm4p
   qlc = 0.5*(3.0*q0-qm1) + (4.0/3.0)*dm4m
   !qlc = 0.5*(3.0*q0+qm1) - (4.0/3.0)*dm4m

   qmin = max(min(q0,qp1,qmd),min(q0,qul,qlc))
   qmax = min(max(q0,qp1,qmd),max(q0,qul,qlc))

   md = 0.5*(sign(1.0,qmin-wq) + sign(1.0,qmax-wq))*min(abs(qmin-wq),abs(qmax-wq))
   wq = wq + md

end subroutine weno9_reconstruction_interface

pure subroutine weno9_poly(P0, P1, P2, P3, P4, b0, b1, b2, b3, b4, qm4, qm3, qm2, qm1, q0, qp1, qp2, qp3, qp4)

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

!> This is the subroutine for the maximum-principle preserving limiter
pure subroutine PP_limiter(wq, q0, wmr, wpl, w0, qmin_g, qmax_g)

   real, intent(in) :: q0         !< tracer concentration in cell i
   real, intent(in) :: wmr, wpl   !< weno reconstruction on the cell interface i-1/2 and i+1/2
   real, intent(in) :: w0         !< w0 : 1st weight of N Gauss-Legendre quadrature weights
                                  !< over [-1/2,1/2],
                                  !< N is weno reconstruction order
   real, intent(in) :: qmin_g, qmax_g !< global min and max of tracer concentration 
                                      !! at the initial time
   real, intent(out) :: wq            !< weno reconstruction at the interface i+1/2 
                                      !! after applying the limiter

   real :: theta, eps, qmin, qmax
   real :: P0, w1, wq_ppm

   wq = 0.0

   P0 = (q0 - w0*(wmr + wpl))/(1.0 - 2.0*w0)

   qmin = min(wmr, P0, wpl)
   qmax = max(wmr, P0, wpl)

   eps = min(1.0e-2, q0)
   !eps = 1.0e-2
   theta = min(abs((qmax_g-q0)/(qmax-q0)), abs((qmin_g-q0+eps)/(qmin-q0)), 1.0)

   wq = theta*(wpl - q0) + q0

end subroutine PP_limiter

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
