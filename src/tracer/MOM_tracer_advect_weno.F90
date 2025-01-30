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
public weno5_reconstruction2

contains

!> 3rd weno reconstruction subroutine and limiter
subroutine weno3_reconstruction(wq,  qm2, qm, q0, qp, qp2, qp3, u, qmin, qmax, wppm)

   real, intent(in) :: qm2, qm, q0, qp, qp2, qp3 ! tracer concentration from i-2 to  i+3 respectively
   real, intent(in) :: u                         ! advection velocity
   real, intent(in) :: qmin, qmax                ! global min and max of tracer concentration
   real, intent(in) :: wppm
   real, intent(out) :: wq                       ! weno reconstruction at the interface i+1/2

   real :: wmr ! wmr : weno reconstruction on the cell interface i-1/2
   real :: wpl ! wpl : weno reconstruction on the cell interface i+1/2
   real :: w0  ! w0 : 1st weight of N Gauss-Legendre quadrature weights over [-1/2,1/2], 
               ! N is weno reconstruction order
   
   !w0 = 1.0
   w0 = 1.0/6.0

   if(u >= 0.0) then 
      call weno3_reconstruction_interface(wmr, qp2, qp, q0, qm, qm2)
      call weno3_reconstruction_interface(wpl, qm2, qm, q0, qp, qp2)
      ! maximum-principle limiter
      !call PP_limiter(wq, q0, wmr, wpl, w0, qmin, qmax)
      call apply_MP2(wq, qm2, qm, q0, qp, qp2, wmr, wpl, w0, qmin, qmax, wppm)
   else
      call weno3_reconstruction_interface(wpl, qp3, qp2, qp, q0, qm)
      call weno3_reconstruction_interface(wmr, qm, q0, qp, qp2, qp3)
      ! maximum-principle limiter
      !call PP_limiter(wq, qp, wmr, wpl, w0, qmin, qmax)
      call apply_MP2(wq, qp3, qp2, qp, q0, qm, wmr, wpl, w0, qmin, qmax, wppm)
   endif

end subroutine weno3_reconstruction

!> 3rd-order weno reconstruction subroutine
subroutine weno3_reconstruction_interface(wq, qmm, qm, q0, qp, qpp)

   real, intent(in) :: qmm, qm, q0, qp, qpp
   real, intent(out) :: wq

   real :: a1, a2, b1, b2, h1, h2, nu
   real :: eps, wnorm, w_1, w_2, P1, P2, tau
   
   call weno3_weights(b1, b2, qm, q0, qp)
   call weno3_poly(P1, P2, qm, q0, qp)

   wq = 0.0
   h1 = 1.0/3.0
   h2 = 2.0/3.0

   ! Alpha values
   eps = 1.0e-20
   tau = abs(b2-b1)
   a1 = h1*(1.0 + (tau/(b1+eps))**2)
   a2 = h2*(1.0 + (tau/(b2+eps))**2)

   ! Normalization
   wnorm = a1+a2
   w_1 = a1 / wnorm
   w_2 = a2 / wnorm

   wq = w_1*P1 + w_2*P2

   ! Monotonicity Preserving
   !call apply_MP(wq, qmm, qm, q0, qp, qpp)

end subroutine weno3_reconstruction_interface

subroutine weno3_poly(P1, P2, qm, q0, qp)


   real, intent(in) :: qm, q0, qp
   real, intent(out) :: P1, P2

   P1 = 0.5*(-qm + 3.0*q0)
   P2 = 0.5*(q0 + qp)

end subroutine weno3_poly

subroutine weno3_weights(b1, b2, qm, q0, qp)

   real, intent(in) :: qm, q0, qp
   real, intent(out) :: b1, b2

   b1 = (q0-qm)*(q0-qm)
   b2 = (qp-q0)*(qp-q0)

end subroutine weno3_weights

!> 5th-order weno reconstruction subroutine and limiter
subroutine weno5_reconstruction(wq, qm2, qm, q0, qp, qp2, qp3, u, qmin, qmax, wppm)

   real, intent(in) :: qm2, qm, q0, qp, qp2, qp3 ! tracer concentration from i-2 to  i+3 respectively
   real, intent(in) :: u                         ! advection velocity
   real, intent(in) :: qmin, qmax                ! global min and max of tracer concentration
   real, intent(in) :: wppm
   real, intent(out) :: wq                       ! weno reconstruction at the interface i+1/2

   real :: wmr ! wmr : weno reconstruction on the cell interface i-1/2
   real :: wpl ! wpl : weno reconstruction on the cell interface i+1/2
   real :: w0  ! w0 : 1st weight of N Gauss-Legendre quadrature weights over [-1/2,1/2], 
               ! N is weno reconstruction order

   !w0 = 5.0/18.0
   w0 = 1.0/12.0

   if(u >= 0.0) then
      call weno5_reconstruction_interface(wmr, qp2, qp, q0, qm, qm2)  ! i-1/2
      call weno5_reconstruction_interface(wpl, qm2, qm, q0, qp, qp2)  ! i+1/2
      ! maximum-principle limiter
      !call PP_limiter(wq, q0, wmr, wpl, w0, qmin, qmax)
      call apply_MP2(wq, qm2, qm, q0, qp, qp2, wmr, wpl, w0, qmin, qmax, wppm)
   else
      call weno5_reconstruction_interface(wpl, qp3, qp2, qp, q0, qm)  
      call weno5_reconstruction_interface(wmr, qm, q0, qp, qp2, qp3)
      ! maximum-principle limiter
      !call PP_limiter(wq, qp, wmr, wpl, w0, qmin, qmax)
      call apply_MP2(wq, qp3, qp2, qp, q0, qm, wmr, wpl, w0, qmin, qmax, wppm)
   endif

end subroutine weno5_reconstruction

subroutine weno5_reconstruction2(wq, qm2, qm, q0, qp, qp2, qp3, u, qmin, qmax, dx)

   real, intent(in) :: qm2, qm, q0, qp, qp2, qp3 ! tracer concentration from i-2 to  i+3 respectively
   real, intent(in) :: u                         ! advection velocity
   real, intent(in) :: qmin, qmax                ! global min and max of tracer concentration
   real, intent(in) :: dx(6)
   real, intent(out) :: wq                       ! weno reconstruction at the interface i+1/2

   real :: wmr ! wmr : weno reconstruction on the cell interface i-1/2
   real :: wpl ! wpl : weno reconstruction on the cell interface i+1/2
   real :: w0  ! w0 : 1st weight of N Gauss-Legendre quadrature weights over [-1/2,1/2],
               ! N is weno reconstruction order
   real :: ds(5), ds1(6)


   !w0 = 5.0/18.0
   w0 = 1.0/12.0

   if(u >= 0.0) then
      ds = dx(1:5)
      call weno5_reconstruction_interface2(wmr, qp2, qp, q0, qm, qm2, ds)  ! i-1/2
      call weno5_reconstruction_interface2(wpl, qm2, qm, q0, qp, qp2, ds)  ! i+1/2
      ! maximum-principle limiter
      !call PP_limiter(wq, q0, wmr, wpl, w0, qmin, qmax)
      !call apply_MP2(wq, qm2, qm, q0, qp, qp2, wmr, wpl, w0, qmin, qmax, wppm)
   else
      ds1 = dx(6:1:-1)
      ds = ds1(1:5)
      call weno5_reconstruction_interface2(wpl, qp3, qp2, qp, q0, qm, ds)
      call weno5_reconstruction_interface2(wmr, qm, q0, qp, qp2, qp3, ds)
      ! maximum-principle limiter
      !call PP_limiter(wq, qp, wmr, wpl, w0, qmin, qmax)
      !call apply_MP2(wq, qp3, qp2, qp, q0, qm, wmr, wpl, w0, qmin, qmax, wppm)
   endif

   wq = wpl

end subroutine weno5_reconstruction2

subroutine weno5_reconstruction_interface(wq, qmm, qm, q0, qp, qpp)

   real, intent(in) :: qmm, qm, q0, qp, qpp
   real, intent(out) :: wq

   real :: a0, a1, a2, b0, b1, b2, d0, d1, d2, nu
   real :: eps,  wnorm, w0, w1, w2, P0, P1, P2, tau
   real :: c0, c1, c2, bb, z0, n0, n1, n2
   integer :: r

   r = 2 

   call weno5_weights(b0, b1, b2, qmm, qm, q0, qp, qpp)
   call weno5_poly(P0, P1, P2, qmm, qm, q0, qp, qpp)

   ! Gamma values in Weno reconstruction
   d0 = 1.0/10.0
   d1 = 6.0/10.0
   d2 = 3.0/10.0

   ! Alpha values
   eps = 1.0e-40
   tau = abs(b2-b0)
   a0 = d0*(1.0 + (tau/(b0+eps))**r)
   a1 = d1*(1.0 + (tau/(b1+eps))**r)
   a2 = d2*(1.0 + (tau/(b2+eps))**r)

   wnorm = 1.0/(a0+a1+a2)
   w0 = a0*wnorm
   w1 = a1*wnorm
   w2 = a2*wnorm

   wq = w0*P0 + w1*P1 + w2*P2

   ! Monotonicity Preserving
   !call apply_MP(wq, qmm, qm, q0, qp, qpp)

end subroutine weno5_reconstruction_interface

subroutine weno5_reconstruction_interface2(wq, qmm, qm, q0, qp, qpp, dx)

   real, intent(in) :: qmm, qm, q0, qp, qpp, dx(5)
   real, intent(out) :: wq

   real :: a0, a1, a2, b0, b1, b2, d0, d1, d2, nu
   real :: eps,  wnorm, w0, w1, w2, P0, P1, P2, tau
   real :: c0, c1, c2, bb, z0, n0, n1, n2
   integer :: r

   r = 2

   call weno5NM_poly(P0, P1, P2, b0, b1, b2, qmm, qm, q0, qp, qpp, dx)

   ! Gamma values in Weno reconstruction
   d0 = 1.0/10.0
   d1 = 6.0/10.0
   d2 = 3.0/10.0

   ! Alpha values
   eps = 1.0e-20
   tau = abs(b2-b0)
   a0 = d0*(1.0 + (tau/(b0+eps))**r)
   a1 = d1*(1.0 + (tau/(b1+eps))**r)
   a2 = d2*(1.0 + (tau/(b2+eps))**r)

   wnorm = 1.0/(a0+a1+a2)
   w0 = a0*wnorm
   w1 = a1*wnorm
   w2 = a2*wnorm

   wq = w0*P0 + w1*P1 + w2*P2

   ! Monotonicity Preserving
   !call apply_MP(wq, qmm, qm, q0, qp, qpp)

end subroutine weno5_reconstruction_interface2

subroutine weno5_weights(b0, b1, b2, qmm, qm, q0, qp, qpp)

   real, intent(in) :: qmm, qm, q0, qp, qpp
   real, intent(out) :: b0, b1, b2

   ! First stencil

   b0 = (13.0/12.0)*(qmm - 2.0*qm + q0)**2 + 0.25*(qmm - 4.0*qm + 3.0*q0)**2 

   ! Second stencil

   b1 = (13.0/12.0)*(qm - 2.0*q0 + qp)**2 + 0.25*(qm - qp)**2

   ! Third stencil

   b2 = (13.0/12.0)*(q0 - 2.0*qp + qpp)**2 + 0.25*(3.0*q0 - 4.0*qp + qpp)**2

end subroutine weno5_weights

subroutine weno5_poly(P0, P1, P2, qmm, qm, q0, qp, qpp)

   real, intent(in) :: qmm, qm, q0, qp, qpp
   real, intent(out) :: P0, P1, P2

   ! First stencil

   P0 = (2.0*qmm - 7.0*qm + 11.0*q0)/6.0

   ! Second stencil

   P1 = (-qm + 5.0*q0 + 2.0*qp)/6.0

   ! Third stencil

   P2 = (2.0*q0 + 5.0*qp - qpp)/6.0

end subroutine weno5_poly

subroutine weno5NM_poly(P0, P1, P2, b0, b1, b2, qmm, qm, q0, qp, qpp, dx)

   real, intent(in) :: qmm, qm, q0, qp, qpp, dx(5)
   real, intent(out) :: P0, P1, P2, b0, b1, b2

   real :: a0, a1, a2, a3, qh, qhh, qhp, qhpp

   ! First stencil

   a0 = dx(1)/(dx(1)+dx(2))
   a1 = dx(2)/(dx(2)+dx(3))

   qh = (1.0-a1)*qm + a1*q0
   qhh = (2.0-a0)*qm - (1.0-a0)*qmm

   b0 = (13.0/12.0)*(2.0*qh - 2.0*qhh)**2 + 0.25*(4.0*q0 - 2.0*qh - 2.0*qhh)**2

   P0 = (6.0*q0 - qh - 2.0*qhh)/3.0

   ! Second stencil

   a2 = dx(3)/(dx(3)+dx(4))

   !qh = (1.0-a1)*qm + a1*q0
   qhp = (1.0-a2)*q0 + a2*qp

   b1 = (13.0/12.0)*(2.0*qh - 4.0*q0 + 2.0*qhp)**2 + 0.25*(-2.0*qh + 2.0*qhp)**2


   P1 = (-qh + 2*q0 + 2.0*qhp)/3.0

   ! Third stencil

   a3 = dx(4)/(dx(4)+dx(5))

   qhpp = (1.0-a3)*qp + a3*qpp

   b2 = (13.0/12.0)*(2.0*qhp - 4.0*qp + 2.0*qhpp)**2 + 0.25*(-6.0*qhp + 8.0*qp - 2.0*qhpp)**2

   P2 = (2.0*qhp + 5.0*qp - qhpp)/3.0

end subroutine weno5NM_poly

subroutine weno7_reconstruction(wq, qm3, qm2, qm1, q0, qp1, qp2, qp3, qp4, u, qmin, qmax, wppm)

   real, intent(in) :: qm3, qm2, qm1, q0, qp1, qp2, qp3, qp4 ! tracer concentration from i-3 to  i+4 respectively
   real, intent(in) :: u                         ! advection velocity
   real, intent(in) :: qmin, qmax                ! global min and max of tracer concentration
   real, intent(in) :: wppm
   real, intent(out) :: wq                       ! weno reconstruction at the interface i+1/2

   real :: wmr ! wmr : weno reconstruction on the cell interface i-1/2
   real :: wpl ! wpl : weno reconstruction on the cell interface i+1/2
   real :: w0  ! w0 : 1st weight of N Gauss-Legendre quadrature weights over [-1/2,1/2], 
               ! N is weno reconstruction order

   !w0 = (322.0-13.0*sqrt(70.0))/1800.0
   w0 = 1.0/20.0

   if(u >= 0.0) then
      call weno7_reconstruction_interface(wmr, qp3, qp2, qp1, q0, qm1, qm2, qm3)
      call weno7_reconstruction_interface(wpl, qm3, qm2, qm1, q0, qp1, qp2, qp3)
      ! maximum-principle limiter
      !call PP_limiter(wq, q0, wmr, wpl, w0, qmin, qmax)
      call apply_MP2(wq, qm2, qm1, q0, qp1, qp2, wmr, wpl, w0, qmin, qmax, wppm)
   else
      call weno7_reconstruction_interface(wpl, qp4, qp3, qp2, qp1, q0, qm1, qm2)
      call weno7_reconstruction_interface(wmr, qm2, qm1, q0, qp1, qp2, qp3, qp4)
      ! maximum-principle limiter
      !call PP_limiter(wq, qp1, wmr, wpl, w0, qmin, qmax)
      call apply_MP2(wq, qp3, qp2, qp1, q0, qm1, wmr, wpl, w0, qmin, qmax, wppm)
   endif

end subroutine weno7_reconstruction

subroutine weno7_reconstruction_interface(wq, qm3, qm2, qm1, q0, qp1, qp2, qp3)

   real, intent(in) :: qm3, qm2, qm1, q0, qp1, qp2, qp3
   real, intent(out) :: wq

   real :: b0, b1, b2, b3, d0, d1, d2, d3
   real :: eps, a0, a1, a2, a3, wnorm, w0, w1, w2, w3, tau
   real :: P0, P1, P2, P3, nu, psi
   integer :: r

   r = 2

   call weno7_weights(b0, b1, b2, b3, qm3, qm2, qm1, q0, qp1, qp2, qp3)
   call weno7_poly(P0, P1, P2, P3, qm3, qm2, qm1, q0, qp1, qp2, qp3)

   d0 = 1.0/35.0
   d1 = 12.0/35.0
   d2 = 18.0/35.0
   d3 = 4.0/35.0

   ! Alpha values
   eps = 1.0e-20
   tau = abs(b3 + 3.0 * b2 - 3.0 * b1 - b0)
   a0 = d0*(1.0 + (tau/(b0+eps))**r)
   a1 = d1*(1.0 + (tau/(b1+eps))**r)
   a2 = d2*(1.0 + (tau/(b2+eps))**r)
   a3 = d3*(1.0 + (tau/(b3+eps))**r)

   ! Normalization
   wnorm = 1.0/(a0+a1+a2+a3)
   w0 = a0*wnorm
   w1 = a1*wnorm
   w2 = a2*wnorm
   w3 = a3*wnorm

   wq = w0*P0 + w1*P1 + w2*P2 + w3*P3
   
   ! Monotonicity Preserving
   !call apply_MP(wq, qm2, qm1, q0, qp1, qp2)

end subroutine weno7_reconstruction_interface

subroutine weno7_poly(P0, P1, P2, P3, qm3, qm2, qm1, q0, qp1, qp2, qp3)

   real, intent(in) :: qm3, qm2, qm1, q0, qp1, qp2, qp3
   real, intent(out) :: P0, P1, P2, P3

   P0 = (-3.0*qm3 + 13.0*qm2 - 23.0*qm1 + 25.0*q0)/12.0

   P1 = (qm2 - 5.0*qm1 + 13.0*q0 + 3.0*qp1)/12.0

   P2 = (-qm1 + 7.0*q0 + 7.0*qp1 - qp2)/12.0

   P3 = (3.0*q0 + 13.0*qp1 - 5.0*qp2 + qp3)/12.0

end subroutine weno7_poly

subroutine weno7_weights(b0, b1, b2, b3, qm3, qm2, qm1, q0, qp1, qp2, qp3)


   real, intent(in) :: qm3, qm2, qm1, q0, qp1, qp2, qp3
   real, intent(out) :: b0, b1, b2, b3

   ! 1st stencil

   b0 = qm3*(547.0*qm3 - 3882.0*qm2 + 4642.0*qm1 - 1854.0*q0) + &
        qm2*(7043.0*qm2 - 17246.0*qm1 + 7042.0*q0) + &
        qm1*(11003.0*qm1 - 9402.0*q0) + 2107.0*q0**2
   
   ! 2nd stencil

   b1 = qm2*(267.0*qm2 - 1642.0*qm1 + 1602.0*q0 - 494.0*qp1) + &
           qm1*(2843.0*qm1 - 5966.0*q0 + 1922.0*qp1) &
           + q0*(3443.0*q0 - 2522.0*qp1) + 547.0*qp1**2
   
   ! 3rd stencil

   b2 = qm1*(547.0*qm1 - 2522.0*q0 + 1922.0*qp1 - 494.0*qp2) + &
           q0*(3443.0*q0 - 5966.0*qp1 + 1602.0*qp2) &
           + qp1*(2843.0*qp1 - 1642.0*qp2) + 267.0*qp2**2
   
   ! 4rd stencil
  
   b3 = q0*(2107.0*q0 - 9402.0*qp1 + 7042.0*qp2 - 1854.0*qp3) + &
           qp1*(11003.0*qp1 - 17246.0*qp2 + 4642.0*qp3) &
           + qp2*(7043.0*qp2 - 3882.0*qp3) + 547.0*qp3**2
   
end subroutine weno7_weights

subroutine weno9_reconstruction(wq, qm4, qm3, qm2, qm1, q0, qp1, qp2, qp3, qp4, qp5, u, qmin, qmax)

   real, intent(in) :: qm4, qm3, qm2, qm1, q0, qp1, qp2, qp3, qp4, qp5 ! tracer concentration from cell i-4 to i+5
   real, intent(in) :: u                         ! advection velocity
   real, intent(in) :: qmin, qmax                ! global min and max of tracer concentration
   real, intent(out) :: wq                       ! weno reconstruction at the interface i+1/2

   real :: wmr ! wmr : weno reconstruction on the cell interface i-1/2
   real :: wpl ! wpl : weno reconstruction on the cell interface i+1/2
   real :: w0  ! w0 : 1st weight of N Gauss-Legendre quadrature weights over [-1/2,1/2], 
               ! N is weno reconstruction order

   w0 = 0.5*0.1294849661688697
   !w0 = 1.0/30.0

   if(u >= 0.0) then
      call weno9_reconstruction_interface(wmr, qp4, qp3, qp2, qp1, q0, qm1, qm2, qm3, qm4)
      call weno9_reconstruction_interface(wpl, qm4, qm3, qm2, qm1, q0, qp1, qp2, qp3, qp4)
      ! maximum-principle limiter
      call PP_limiter(wq, q0, wmr, wpl, w0, qmin, qmax)
   else
      call weno9_reconstruction_interface(wpl, qp5, qp4, qp3, qp2, qp1, q0, qm1, qm2, qm3)
      call weno9_reconstruction_interface(wmr, qm3, qm2, qm1, q0, qp1, qp2, qp3, qp4, qp5)
      ! maximum-principle limiter
      call PP_limiter(wq, qp1, wmr, wpl, w0, qmin, qmax)
   endif

end subroutine weno9_reconstruction

subroutine weno9_reconstruction_interface(wq, qm4, qm3, qm2, qm1, q0, qp1, qp2, qp3, qp4)

   real, intent(in) :: qm4, qm3, qm2, qm1, q0, qp1, qp2, qp3, qp4
   real, intent(out) :: wq

   real :: b0, b1, b2, b3, b4, d0, d1, d2, d3, d4
   real :: eps, a0, a1, a2, a3, a4, wnorm, w0, w1, w2, w3, w4, tau
   real :: P0, P1, P2, P3, P4
   integer :: r

   r = 2

   call weno9_poly(P0, P1, P2, P3, P4, b0, b1, b2, b3, b4, qm4, qm3, qm2, qm1, q0, qp1, qp2, qp3, qp4)

   d0 = 1.0/126.0
   d1 = 10.0/63.0
   d2 = 10.0/21.0
   d3 = 20.0/63.0
   d4 = 5.0/126.0

   ! Alpha values
   eps = 1.0e-20
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

end subroutine weno9_reconstruction_interface

subroutine weno9_poly(P0, P1, P2, P3, P4, b0, b1, b2, b3, b4, qm4, qm3, qm2, qm1, q0, qp1, qp2, qp3, qp4)

   real, intent(in) :: qm4, qm3, qm2, qm1, q0, qp1, qp2, qp3, qp4
   real, intent(out) :: P0, P1, P2, P3, P4, b0, b1, b2, b3, b4

   real :: qx, qx2, qx3, qx4, L1, L2, L3, L4

   L1 = 0.5 ; L2 = 1.0/6.0 ; L3 = 1.0/20.0 ; L4 = 1.0/70.0

   ! stencil 1
   qx = (27.0*qm4 - 146.0*qm3 + 336.0*qm2 - 462.0*qm1 + 245.0*q0)/120.0
   qx2 = (25.0*qm4 - 128.0*qm3 + 262.0*qm2 - 240.0*qm1 + 81.0*q0)/56.0
   qx3 = (3.0*qm4 - 14.0*qm3 + 24.0*qm2 - 18.0*qm1 + 5.0*q0)/12.0
   qx4 = (qm4 - 4.0*qm3 + 6.0*qm2 - 4.0*qm1 + q0)/24.0

   P0 = q0 + qx*L1 + qx2*L2 + qx3*L3 + qx4*L4
   b0 = (qx+ 0.1*qx3)**2 + (13.0/3.0)*(qx2 + (123.0/455.0)*qx4)**2 + &
          (781.0/20.0)*(qx3)**2 + (1421461.0/2275.0)*(qx4)**4

   ! stencil 2
   qx = (-11.0*qm3 + 66.0*qm2 - 192.0*qm1 + 110.0*q0 + 27.0*qp1)/120.0
   qx2 = (-3.0*qm3 + 12.0*qm2 + 10.0*qm1 - 44.0*q0 + 25.0*qp1)/56.0
   qx3 = (qm3 - 6.0*qm2 + 12.0*qm1 - 10.0*q0 + 3.0*qp1)/12.0
   qx4 = (qm3 - 4.0*qm2 + 6.0*qm1 - 4.0*q0 + qp1)/24.0

   P1 = q0 + qx*L1 + qx2*L2 + qx3*L3 + qx4*L4
   b1 = (qx+ 0.1*qx3)**2 + (13.0/3.0)*(qx2 + (123.0/455.0)*qx4)**2 + &
           (781.0/20.0)*(qx3)**2 + (1421461.0/2275.0)*(qx4)**4

   ! stencil 3
   qx = (11.0*qm2 - 82.0*qm1 + 82.0*qp1 - 11.0*qp2)/120.0
   qx2 = (-3.0*qm2 + 40.0*qm1 - 74.0*q0 + 40.0*qp1 - 3.0*qp2)/56.0
   qx3 = (-qm2 + 2.0*qm1 - 2.0*qp1 + qp2)/12.0
   qx4 = (qm2 - 4.0*qm1 + 6.0*q0 - 4.0*qp1 + qp2)/24.0

   P2 = q0 + qx*L1 + qx2*L2 + qx3*L3 + qx4*L4
   b2 = (qx+ 0.1*qx3)**2 + (13.0/3.0)*(qx2 + (123.0/455.0)*qx4)**2 + &
           (781.0/20.0)*(qx3)**2 + (1421461.0/2275.0)*(qx4)**4

   ! stencil 4
   qx = (-27.0*qm1 - 110.0*q0 + 192.0*qp1 - 66.0*qp2 + 11.0*qp3)/120.0
   qx2 = (25.0*qm1 - 44.0*q0 + 10.0*qp1 + 12.0*qp2 - 3.0*qp3)/56.0
   qx3 = (-3.0*qm1 + 10.0*q0 - 12.0*qp1 + 6.0*qp2 - qp3)/12.0
   qx4 = (qm1 - 4.0*q0 + 6.0*qp1 - 4.0*qp2 + qp3)/24.0

   P3 = q0 + qx*L1 + qx2*L2 + qx3*L3 + qx4*L4
   b3 = (qx+ 0.1*qx3)**2 + (13.0/3.0)*(qx2 + (123.0/455.0)*qx4)**2 + &
           (781.0/20.0)*(qx3)**2 + (1421461.0/2275.0)*(qx4)**4

   ! stencil 5
   qx = (-245.0*q0 + 462.0*qp1 - 336.0*qp2 + 146.0*qp3 - 27.0*qp4)/120.0
   qx2 = (81.0*q0 - 240.0*qp1 + 262.0*qp2 - 128.0*qp3 + 25.0*qp4)/56.0
   qx3 = (-5.0*q0 + 18.0*qp1 - 24.0*qp2 + 14.0*qp3 - 3.0*qp4)/24.0
   qx4 = (q0 - 4.0*qp1 + 6.0*qp2 - 4.0*qp3 + qp4)/24.0

   P4 = q0 + qx*L1 + qx2*L2 + qx3*L3 + qx4*L4
   b4 = (qx+ 0.1*qx3)**2 + (13.0/3.0)*(qx2 + (123.0/455.0)*qx4)**2 + &
           (781.0/20.0)*(qx3)**2 + (1421461.0/2275.0)*(qx4)**4

end subroutine weno9_poly

!> This is the subroutine for the monotonic-preserving limiter
subroutine apply_MP(wq, qmm, qm, q0, qp, qpp)

   real, intent(in) :: qmm, qm, q0, qp, qpp ! tracer concentration from i-2 to i+2 respectively
   real, intent(inout) :: wq                ! weno reconstruction at the interface i+1/2 after MP limiter

   real :: d0, d1, dm1, beta, ka
   real :: dm4, qlc, qmd, qul, qmin, qmax, md, dm
   real :: cp, qmp, mmp, tmp

   ka = 2.0
   call minmod(mmp,qp-q0, ka*(q0-qm))
   qmp = q0 + mmp

   !cp = (wq - q0)*(qmp - q0) 

   cp = (wq - q0)*(wq - qmp)

   if(cp <= 1.0e-10) then

      tmp = wq

   else

      dm1 = qmm - 2.0*qm + q0
      d0 = qp - 2.0*q0 + qm
      d1 = qpp - 2.0*qp + q0

      call minmod3(dm4, 4.0*d0-d1, 4.0*d1-d0, d0, d1)
      call minmod3(dm, 4.0*dm1-d0, 4.0*d0-dm1, dm1, d0)
    
      beta = 4.0
      qul = q0 + ka*(q0-qm)
      qmd = 0.5*(q0 + qp) - 0.5*dm4
      qlc = q0 + 0.5*(q0-qm) + (beta/3.0)*dm
      !qlc = 0.5*(q0 + qul) + 0.5*ka*dm

      qmin = max(min(q0,qp,qmd),min(q0,qul,qlc))
      qmax = min(max(q0,qp,qmd),max(q0,qul,qlc))
   
      call median(md,wq,qmin,qmax)
      tmp = md

      !tmp = q0

   endif

   wq = tmp 

end subroutine apply_MP

!> This is the subroutine for the monotonic-preserving limiter
subroutine apply_MP2(wq, qmm, qm, q0, qp, qpp, wmr, wpl, w0, qmin_g, qmax_g, wppm)

   real, intent(in) :: qmm, qm, q0, qp, qpp ! tracer concentration from i-2 to i+2 respectively
   real, intent(out) :: wq                ! weno reconstruction at the interface i+1/2 after MP limiter
   real, intent(in) :: wmr, wpl   ! weno reconstruction on the cell interface i-1/2 and i+1/2
   real, intent(in) :: w0         ! w0 : 1st weight of N Gauss-Legendre quadrature weights over [-1/2,1/2],
                                  ! N is weno reconstruction order
   real, intent(in) :: qmin_g, qmax_g ! global min and max of tracer concentration at the initial time
   real, intent(in) :: wppm

   real :: d0, d1, dm1, beta, ka
   real :: dm4, qlc, qmd, qul, qmin, qmax, md, dm
   real :: cp, qmp, mmp, tmp
   real :: qmin2, qmax2, theta, eps, P0, s0

   ka = 2.0
   call minmod(mmp,qp-q0, ka*(q0-qm))
   qmp = q0 + mmp

   !cp = (wq - q0)*(qmp - q0)

   cp = (wpl - q0)*(wpl - qmp)

   !if(cp <= 0.0) then

   !   tmp = wpl

   !else

      dm1 = qmm - 2.0*qm + q0
      d0 = qp - 2.0*q0 + qm
      d1 = qpp - 2.0*qp + q0

      call minmod3(dm4, 4.0*d0-d1, 4.0*d1-d0, d0, d1)
      call minmod3(dm, 4.0*dm1-d0, 4.0*d0-dm1, dm1, d0)

      beta = 4.0
      qul = q0 + ka*(q0-qm)
      qmd = 0.5*(q0 + qp) - 0.5*dm4
      qlc = q0 + 0.5*(q0-qm) + (beta/3.0)*dm
      !qlc = 0.5*(q0 + qul) + 0.5*ka*dm

      qmin2 = max(min(q0,qp,qmd),min(q0,qul,qlc))
      qmax2 = min(max(q0,qp,qmd),max(q0,qul,qlc))

      call median(md,wpl,qmin2,qmax2)
      tmp = md

      !tmp = wppm

   !endif

   ! MPP limiter
   P0 = (q0 - w0*wmr - w0*wpl)/(1.0 - 2.0*w0)

   qmin = min(wmr, P0, wpl)
   qmax = max(wmr, P0, wpl)

   !print*, wmr, wpl, qmin, qmax

   !eps = min(1.0e-13, (q0 - qmin)**5, q0)
   !eps = min(1.0e-13, q0)
   eps = 1.0e-13
   if(abs(qmin-q0) <= eps .or. abs(qmax-q0) <= eps) then
      theta = 0.0
   else
      theta = min(abs((qmax_g-q0)/(qmax-q0)), abs((qmin_g-q0+eps)/(qmin-q0)), 1.0)
      !theta = min(abs((qmax_g-q0)/(qmax-q0)), abs((qmin_g-q0)/(qmin-q0)), 1.0)
   endif
   
   !if(tmp <= 0.0 ) then
   !   theta = 0.0
   !else
   !   theta = 1.0
   !endif

   !if(qmin < eps) then
   !   theta = 0.0
   !else 
   !   s0 = abs(q0/(q0 - qmin))
   !   !if (s0 > 1.0) s0 = 0.0
   !   theta = min(s0, 1.0)
   !endif
   
   !s0 = q0

   !if(theta == 1.0 .and. tmp < 0.0) then
   !  print*, 'theta = ', theta, q0
   !  print*, 'min1 = ', qmin, tmp
   !endif 

   !if(theta < 1.0) print*,'theta = ', theta
   !theta = 0.0
   wq = theta*(tmp - q0) + q0

   !wq = tmp

end subroutine apply_MP2

!> This is the subroutine for the maximum-principle preserving limiter
subroutine PP_limiter(wq, q0, wmr, wpl, w0, qmin_g, qmax_g)

   real, intent(in) :: q0 ! tracer concentration in cell i
   real, intent(in) :: wmr, wpl   ! weno reconstruction on the cell interface i-1/2 and i+1/2
   real, intent(in) :: w0         ! w0 : 1st weight of N Gauss-Legendre quadrature weights over [-1/2,1/2], 
                                  ! N is weno reconstruction order
   real, intent(in) :: qmin_g, qmax_g ! global min and max of tracer concentration at the initial time
   real, intent(out) :: wq            ! weno reconstruction at the interface i+1/2 after applying the limiter

   real :: qmin, qmax, theta, eps
   real :: P0, w1, wq_ppm

   wq = 0.0

   P0 = (q0 - w0*wmr - w0*wpl)/(1.0 - 2.0*w0)

   qmin = min(wmr, P0, wpl)
   qmax = max(wmr, P0, wpl)

   !eps = min(1.0e-13, (q0 - qmin)**5, q0)
   !eps = min(1.0e-13, q0)
   eps = 1.0e-13
   theta = min(abs((qmax_g-q0)/(qmax-q0)), abs((qmin_g-q0+eps)/(qmin-q0)), 1.0)
   !theta = min(abs((qmax_g-q0)/(qmax-q0)), 1.0)

   wq = theta*(wpl - q0) + q0

   !wq = wpl

end subroutine PP_limiter

!> This subroutine computes the median of three variables
subroutine median(md,a,b,c)

   real, intent(in) :: a,b,c
   real, intent(out) :: md

   real :: mm, d1, d2

   d1 = b - a
   d2 = c - a
   
   call minmod(mm,d1,d2)
   md = a + mm

end subroutine median

!> This subroutine computes the minmod of two variables
subroutine minmod(mab,a,b)

   real, intent(in) :: a,b
   real, intent(out) :: mab
   
   real :: s1,s2

   mab = 0.5*(sign(1.0,a) + sign(1.0,b))*min(abs(a),abs(b))
   
end subroutine minmod 

!> This subroutine computes the minmod of four variables
subroutine minmod3(mab,a1,a2,a3,a4)

   real, intent(in) :: a1,a2,a3,a4
   real, intent(out) :: mab
   real :: s1,s2,s3,s4

   s1 = 0.5*(sign(1.0,a1) + sign(1.0,a2))*min(abs(a1),abs(a2))

   s2 = 0.5*(sign(1.0,a3) + sign(1.0,a4))*min(abs(a3),abs(a4))

   mab = 0.5*(sign(1.0,s1) + sign(1.0,s2))*min(abs(s1),abs(s2))

end subroutine minmod3

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
    Tr => Reg%Tr(m)

    do k=1,nz

      call array_global_min_max(Tr%t(:,:,k), G, 1, Tmin, Tmax)

      Tr%Tmingg = Tmin
      Tr%Tmaxgg = Tmax

      !print*, m, k, Tr%Tmingg, Tr%Tmaxgg

    enddo
  enddo

end subroutine tracer_min_max_init


subroutine PPM_reconstruction(wq, Tm1, T0, Tp1, Tp2, mu, u, mask)


  real, intent(in) :: Tm1, T0, Tp1, Tp2, mu, u, mask 
  real, intent(out) :: wq

  real :: aL, aR, dA, mA, a6
  real :: Tm, Tc, Tp

  if(u >= 0.0) then
     Tm = Tm1; Tc = T0; Tp = Tp1
  else
     Tm = T0; Tc = Tp1; Tp = Tp2
  endif  


  aL = ( 5.*Tc + ( 2.*Tm - Tp ) )/6. ! H3 estimate
  aL = max( min(Tc,Tm), aL) ; aL = min( max(Tc,Tm), aL) ! Bound
  aR = ( 5.*Tc + ( 2.*Tp - Tm ) )/6. ! H3 estimate
  aR = max( min(Tc,Tp), aR) ; aR = min( max(Tc,Tp), aR) ! Bound

  dA = aR - aL ; mA = 0.5*( aR + aL )
  if (mask*(Tp-Tc)*(Tc-Tm) <= 0.) then
     aL = Tc ; aR = Tc ! PCM for local extrema and boundary cells
  elseif ( dA*(Tc-mA) > (dA*dA)/6. ) then
     aL = (3.*Tc) - 2.*aR
  elseif ( dA*(Tc-mA) < - (dA*dA)/6. ) then
     aR = (3.*Tc) - 2.*aL
  endif

  a6 = 6.*Tc - 3. * (aR + aL) ! Curvature

  if (u >= 0.0) then
     wq = aR - 0.5*mu*((aR - aL) - a6*(1.0 - (2.0/3.0)*mu))
  else
     wq = aL + 0.5*mu*((aR - aL) + a6*(1.0 - (2.0/3.0)*mu))
  endif

end subroutine PPM_reconstruction

!> \namespace mom_tracer_advect
!!
!!    This program contains the subroutines that advect tracers
!!  horizontally (i.e. along layers) using high-order WENO schemes (Balsara et al., 2016) 
!!  using the Z-type smoothness indicators (Borges et al., 2008). 
!!  We followed Suresh & Huynh (1997) and Balsara & Shu (2000 for the monotonicity preserving used along with the WENO schemes.
!!
!!  This scheme conserves the total amount of tracer while avoiding
!!  spurious maxima and minima of the tracer concentration

end module MOM_tracer_advect_weno
