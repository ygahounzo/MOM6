!>  This module contains the subroutines of the WENO schemes that advect tracers along coordinate surfaces.
module MOM_tracer_advect_weno

! This file is part of MOM6. See LICENSE.md for the license.

implicit none ; private

#include <MOM_memory.h>

public weno5_reconstruction
public weno7_reconstruction
public PPM_reconstruction

contains

!> 5th-order weno reconstruction subroutine
pure subroutine weno5_reconstruction(wq, q, u, cfl)
  real, intent(in) :: q(5)   !< tracer concentration from cell i-2 to i+2 [conc]
  real, intent(in) :: u      !< advective flux [H L2 ~> m3 or kg]
  real, intent(in) :: cfl(3) !< absolute value of the advective upwind-cell CFL number [nondim]
  real, intent(out) :: wq    !< weno tracer concentration at the cell interface i+1/2  [conc]

  if (u >= 0.0) then
    call weno5z_reconstruction_interface(wq, q(1), q(2), q(3), q(4), q(5), cfl, u)
  else
    call weno5z_reconstruction_interface(wq, q(5), q(4), q(3), q(2), q(1), cfl, u)
  endif

end subroutine weno5_reconstruction

!> 5th-order weno z-type reconstruction flux
pure subroutine weno5z_reconstruction_interface(wpl, qmm, qm, q0, qp, qpp, cfl, u)
  real, intent(in) :: qmm, qm, q0, qp, qpp !< tracer concentration [conc]
  real, intent(in) :: cfl(3) !< absolute value of the advective upwind-cell CFL number [nondim]
  real, intent(in) :: u      !< advective flux [H L2 ~> m3 or kg]
  real, intent(out) :: wpl   !< weno tracer concentrations at the cell interface +1/2 [conc]

  real :: P0, P1, P2         ! reconstructed polynomials
  real :: b0, b1, b2         ! smoothness indicator
  real :: w0, w1, w2         ! nonlinear weights
  real :: tau                ! Difference of smoothness indicators
  real, parameter :: C1_6 = 1.0/6.0             ! The ration of 1/6 [nondim]
  real, parameter :: d0 = 1.0/10.0              ! The ratio of 1/10 [nondim]
  real, parameter :: d1 = 6.0/10.0              ! The ratio of 3/5 [nondim]
  real, parameter :: d2 = 3.0/10.0              ! The ratio of 3/10 [nondim]
  real :: wnorm                                 ! Temporary variable
  real :: dm1, dd0, dd1, dm4p, dm4m             ! Temporary variables
  real :: qul, qmd, qlc, qmin, qmax, alpha      ! Temporary variables
  real :: qmp, Tm, Tp, wpm, sg, q0_min, q0_max  ! Temporary variables

  if (((maxval(cfl) > 0.4) .or. (abs(maxval(cfl)-minval(cfl)) > 0.1))) then
    if (u < 0.0) then
      Tm = qp ; Tp = qm
    else
      Tm = qm ; Tp = qp
    endif
    call PPM_reconstruction(wpl, Tm, q0, Tp, u, cfl(2), 1.)
    return
  endif

  ! Compute flux at left side of i+1/2
  ! First stencil
  P0 = ((2.0*qmm - 7.0*qm) + 11.0*q0)*C1_6
  b0 = (qmm*((4.0*qmm - 19.0*qm) + 11.0*q0)) + (qm*(25.0*qm - 31.0*q0) + 10.0*(q0*q0))

  ! Second stencil
  P1 = ((-qm + 5.0*q0) + 2.0*qp)*C1_6
  b1 = (qm*((4.0*qm - 13.0*q0) + 5.0*qp)) + (q0*(13.0*q0 - 13.0*qp) + 4.0*(qp*qp))

  ! Third stencil
  P2 = ((2.0*q0 + 5.0*qp) - qpp)*C1_6
  b2 = (q0*((10.0*q0 - 31.0*qp) + 11.0*qpp)) + (qp*(25.0*qp - 19.0*qpp) + 4.0*(qpp*qpp))

  ! Alpha values
  tau = abs(b2-b0)
  w0 = d0*weight_fac(tau, b0)
  w1 = d1*weight_fac(tau, b1)
  w2 = d2*weight_fac(tau, b2)

  wnorm = 1./((w0 + w1) + w2)
  w0 = w0*wnorm
  w1 = w1*wnorm
  w2 = w2*wnorm

  wpl = ((w0 * P0) + (w1 * P1)) + (w2 * P2)

  ! Apply the monotonicity preserving limiter based on the hybrid approach from He et al. (2016)
  alpha   = (1.0 - cfl(2))/cfl(2)
  qul = q0 + (alpha*(q0-qm))
  qmp = q0 + (minmod2((qp-q0),(qul-q0)))

  dm1 = (qmm - 2.*qm) + q0
  dd0 = (qp  - 2.*q0) + qm
  dd1 = (q0  - 2.*qp) + qpp

  dm4p = minmod4( (4.*dd0 - dd1), (4.*dd1 - dd0), dd0, dd1 )
  dm4m = minmod4( (4.*dm1 - dd0), (4.*dd0 - dm1), dm1, dd0 )

  qmd = 0.5*((qp  + q0) - dm4p)
  qlc = 0.5*((qul + q0) + (dm4m*alpha))

  qmin = max(min(q0,qp,qmd),min(q0,qul,qlc))
  qmax = min(max(q0,qp,qmd),max(q0,qul,qlc))
  q0_min = min(q0,qmp) ; q0_max = max(q0,qmp)
  sg = sign(1.,((wpl-qmin)*(wpl-qmax)))

  if ((sg > 0.0) .or. ((qmax-qmin) > (q0_max-q0_min))) then
    if (u < 0.0) then
      Tm = qp ; Tp = qm
    else
      Tm = qm ; Tp = qp
    endif
    call PPM_reconstruction(wpl, Tm, q0, Tp, u, cfl(2), 1.)
  endif

end subroutine weno5z_reconstruction_interface

!> 7th-order weno reconstruction subroutine
pure subroutine weno7_reconstruction(wq, q, u, cfl)
  real, intent(in) :: q(7)  !< tracer concentration from i-3 to i+3 [conc]
  real, intent(in) :: u     !< advective flux [H L2 ~> m3 or kg]
  real, intent(in) :: cfl(3)!< absolute value of the advective upwind-cell CFL number [nondim]
  real, intent(out) :: wq   !< weno tracer concentration at the cell interface i+1/2  [conc]

  if (u >= 0.0) then
    call weno7z_reconstruction_interface(wq, q(1), q(2), q(3), q(4), q(5), q(6), q(7), cfl, u)
  else
    call weno7z_reconstruction_interface(wq, q(7), q(6), q(5), q(4), q(3), q(2), q(1), cfl, u)
  endif

end subroutine weno7_reconstruction

!> 7th-order weno z-type reconstruction flux
pure subroutine weno7z_reconstruction_interface(wpl, qm3, qm2, qm1, q0, qp1, qp2, qp3, cfl, u)
  real, intent(in) :: qm3, qm2, qm1, q0, qp1, qp2, qp3 !< tracer concentration [conc]
  real, intent(in) :: cfl(3)  !< absolute value of the advective upwind-cell CFL number [nondim]
  real, intent(in) :: u       !< advective flux [H L2 ~> m3 or kg]
  real, intent(out) :: wpl    !< weno tracer concentrations at the cell interface +1/2 [conc]

  real :: P0, P1, P2, P3     ! reconstructed polynomials
  real :: b0, b1, b2, b3     ! smoothness indicator
  real :: w0, w1, w2, w3     ! nonlinear weights
  real :: tau                ! Difference of smoothness indicators
  real, parameter :: C1_12 = 1.0/12.0  ! [nondim]
  real, parameter :: d0 = 1.0/35.0     ! The ratio of 1/35 [nondim]
  real, parameter :: d1 = 12.0/35.0    ! The ratio of 12/35 [nondim]
  real, parameter :: d2 = 18.0/35.0    ! The ratio of 18/35 [nondim]
  real, parameter :: d3 = 4.0/35.0     ! The ratio of 4/35 [nondim]
  real :: wnorm                        ! Temporary variable
  real :: dm1, dd0, dd1, dm4p, dm4m         ! Temporary variables
  real :: qul, qmd, qlc, qmin, qmax, alpha  ! Temporary variables
  real :: sg, Tm, Tp, q0_min, q0_max, qmp   ! Temporary variables

  if (((maxval(cfl) > 0.4) .or. (abs(maxval(cfl)-minval(cfl)) > 0.1))) then
    if (u < 0.0) then
      Tm = qp1 ; Tp = qm1
    else
      Tm = qm1 ; Tp = qp1
    endif
    call PPM_reconstruction(wpl, Tm, q0, Tp, u, cfl(2), 1.)
    return
  endif

  ! 1st stencil
  P0 = (((-3.0*qm3 + 13.0*qm2) - 23.0*qm1) + 25.0*q0)*C1_12
  b0 = ((qm3*((547.0*qm3 - 3882.0*qm2) + (4642.0*qm1 - 1854.0*q0))) + &
      (qm2*((7043.0*qm2 - 17246.0*qm1) + 7042.0*q0))) + &
      ((qm1*(11003.0*qm1 - 9402.0*q0)) + 2107.0*(q0*q0))

  ! 2nd stencil
  P1 = (((qm2 - 5.0*qm1) + 13.0*q0) + 3.0*qp1)*C1_12
  b1 = ((qm2*((267.0*qm2 - 1642.0*qm1) + (1602.0*q0 - 494.0*qp1))) + &
         (qm1*((2843.0*qm1 - 5966.0*q0) + 1922.0*qp1))) + &
         ((q0*(3443.0*q0 - 2522.0*qp1)) + 547.0*(qp1*qp1))

  ! 3rd stencil
  P2 = (((-qm1 + 7.0*q0) + 7.0*qp1) - qp2)*C1_12
  b2 = ((qm1*((547.0*qm1 - 2522.0*q0) + (1922.0*qp1 - 494.0*qp2))) + &
         (q0*((3443.0*q0 - 5966.0*qp1) + 1602.0*qp2))) + &
         ((qp1*(2843.0*qp1 - 1642.0*qp2)) + 267.0*(qp2*qp2))

  ! 4rd stencil
  P3 = (((3.0*q0 + 13.0*qp1) - 5.0*qp2) + qp3)*C1_12
  b3 = ((q0*((2107.0*q0 - 9402.0*qp1) + (7042.0*qp2 - 1854.0*qp3))) + &
         (qp1*((11003.0*qp1 - 17246.0*qp2) + 4642.0*qp3))) + &
         ((qp2*(7043.0*qp2 - 3882.0*qp3)) + 547.0*(qp3*qp3))

  ! Alpha values
  tau = abs((b0-b3) + (3.*(b1-b2)))
  w0 = d0*weight_fac(tau, b0)
  w1 = d1*weight_fac(tau, b1)
  w2 = d2*weight_fac(tau, b2)
  w3 = d3*weight_fac(tau, b3)

  ! Normalization
  wnorm = 1./((w0 + w1) + (w2 + w3))
  w0 = w0 * wnorm
  w1 = w1 * wnorm
  w2 = w2 * wnorm
  w3 = w3 * wnorm

  wpl = (((w0 * P0) + (w1 * P1)) + (w2 * P2)) + (w3 * P3)

  ! Apply the monotonicity preserving limiter based on the hybrid approach from He et al. (2016)
  alpha   = (1. - cfl(2))/cfl(2)
  qul = q0 + (alpha*(q0-qm1))
  qmp = q0 + (minmod2((qp1-q0),(qul-q0)))

  dm1 = (qm2 - 2.*qm1) + q0
  dd0 = (qp1 - 2.*q0)  + qm1
  dd1 = (q0  - 2.*qp1) + qp2

  dm4p = minmod4( (4.*dd0 - dd1), (4.*dd1 - dd0), dd0, dd1 )
  dm4m = minmod4( (4.*dm1 - dd0), (4.*dd0 - dm1), dm1, dd0 )

  qmd = 0.5*((qp1 + q0) - dm4p)
  qlc = 0.5*((qul + q0) + (dm4m*alpha))

  qmin = max(min(q0,qp1,qmd),min(q0,qul,qlc))
  qmax = min(max(q0,qp1,qmd),max(q0,qul,qlc))
  q0_min = min(q0,qmp) ; q0_max = max(q0,qmp)
  sg = sign(1.,((wpl-qmin)*(wpl-qmax)))

  if ((sg > 0.0) .or. ((qmax-qmin) > (q0_max-q0_min))) then
    call weno5z_reconstruction_interface(wpl, qm2, qm1, q0, qp1, qp2, cfl, u)
  endif

end subroutine weno7z_reconstruction_interface

!> ppm reconstruction flux
pure subroutine PPM_reconstruction(wq_ppm, qm, q0, qp, u, cfl, qext)
  real, intent(in) :: qm, q0, qp !< tracer concentration for 3-stencil wide [conc]
  real, intent(in) :: u        !< advective flux [H L2 ~> m3 or kg]
  real, intent(in) :: cfl      !< absolute value of the advective upwind-cell CFL number [nondim]
  real, intent(in) :: qext     !< check local extrema
  real, intent(out) :: wq_ppm  !< PPM tracer concentration at the cell interface i+1/2  [conc]

  real :: aL, aR, dA, mA, a6 ! local variables

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

  a6 = 6.*q0 - 3. * (aR + aL) ! Curvature

  if (u >= 0.0) then
    wq_ppm = (aR - 0.5 * cfl * ((aR - aL) - a6 * (1. - 2./3. * cfl)))
  else
    wq_ppm = (aL + 0.5 * cfl * ((aR - aL) + a6 * (1. - 2./3. * cfl)))
  endif

end subroutine PPM_reconstruction

pure elemental function minmod2(a,b) result(r)
  real, intent(in) :: a, b !< values to find the minmod for
  real :: r   ! minmod value

  ! 0 if opposite sign or either is zero; otherwise sign(a)*min(|a|,|b|)
  if ((a * b) <= 0.0) then
    r = 0.0
  else
    r = sign(min(abs(a),abs(b)), a)
  end if
end function minmod2

pure elemental function minmod4(a,b,c,d) result(r)
  real, intent(in) :: a, b, c, d !< values to find the minmod for
  real :: r ! minmod values

  if ( ((a*b) <= 0.0) .or. ((a*c) <= 0.0) .or. ((a*d) <= 0.0) ) then
    r = 0.
  else
    r = sign( min( min(abs(a),abs(b)), min(abs(c),abs(d)) ), a )
  end if
end function minmod4

!> Compute the factor for the WENO weights
pure function weight_fac(tau, b) result(factor)
  real, intent(in)  :: tau  !< Difference of the smoothness indicator [A ~> a]
  real, intent(in)  :: b    !< The smoothness indicator [A ~> a]
  real :: factor            !< The factor for the weight [nondim]

  ! factor = (1. + (tau / (b + 1.0e-20))**2)
  factor = 1.0e40; if (abs(b) > 1.0e-20*tau) factor = (1.0 + (tau / b)**2)

end function weight_fac

!> \namespace mom_tracer_advect
!!
!!  This program contains the subroutines that advect tracers
!!  horizontally (i.e. along layers) using high-order WENO schemes (Balsara et al., 2016)
!!  using the Z-type smoothness indicators (Borges et al., 2008).
!!  For monotonicity preservation we switch from WENO7 to WENO5 and then
!!  to PPM:H3 where necessary, i.e. this is a hybrid WENO scheme (He et al., 2016).
!!  The hybrid approach applies WENO to retain high-order accuracy, while reverting to PPM:H3
!!  only when monotonicity constraints are violated or when tracer CFL conditions exceed
!!  the monotonicity-preserving (MP) stability limit (CFL > 0.4).
!!  In practice, this fallback is not frequent when DT_TRACER_ADVECT < DT_THERM.
!!
!!  This scheme conserves the total amount of tracer while avoiding
!!  spurious maxima and minima of the tracer concentration

end module MOM_tracer_advect_weno
