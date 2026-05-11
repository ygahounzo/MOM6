!>  This module contains the subroutines of the WENO schemes that advect tracers along coordinate surfaces.
module MOM_tracer_advect_weno

! This file is part of MOM6. See LICENSE.md for the license.

	implicit none ; private

#include <MOM_memory.h>

	public weno5_reconstruction
	public weno7_reconstruction
	public PPM_reconstruction

contains

!> 5th-order weno z-type reconstruction flux
pure subroutine weno5_reconstruction(wq, q, u, cfl)
	real, intent(in) :: q(7)   !< tracer concentration from cell i-3 to i+3 [conc]
	real, intent(in) :: u      !< advective flux [H L2 ~> m3 or kg]
	real, intent(in) :: cfl(3) !< absolute value of the advective upwind-cell CFL number [nondim]
	real, intent(out) :: wq    !< weno flux  [conc]

	real :: P0, P1, P2         ! reconstructed polynomials
	real :: b0, b1, b2         ! smoothness indicator
	real :: w0, w1, w2         ! nonlinear weights
	real :: tau                ! Difference of smoothness indicators
	real, parameter :: C1_6 = 1.0/6.0             ! The ration of 1/6 [nondim]
	real, parameter :: d0 = 1.0/10.0              ! The ratio of 1/10 [nondim]
	real, parameter :: d1 = 6.0/10.0             ! The ratio of 3/5 [nondim]
	real, parameter :: d2 = 3.0/10.0             ! The ratio of 3/10 [nondim]
	real :: wnorm                                 ! Temporary variable
	real :: dm1, dd0, dd1, dm4p, dm4m             ! Temporary variables
	real :: qul, qmd, qlc, qmin, qmax, alpha      ! Temporary variables
	real :: qmp, Tm, Tp, wpm, sg, q0_min, q0_max  ! Temporary variables
	real :: dm2, dd2, delta, cfl_diff, wpl, wmr, dA, mA, a6, disct
	logical :: lim, disc

	lim = .false.

	! Left state at i+1/2
	P0 = ((2.0*q(2) - 7.0*q(3)) + 11.0*q(4))*C1_6
	b0 = q(2)*((4.0*q(2) - 19.0*q(3)) + 11.0*q(4)) + (q(3)*(25.0*q(3) - 31.0*q(4)) + 10.0*(q(4)*q(4)))

	P1 = ((-q(3) + 5.0*q(4)) + 2.0*q(5))*C1_6
	b1 = q(3)*((4.0*q(3) - 13.0*q(4)) + 5.0*q(5)) + (q(4)*(13.0*q(4) - 13.0*q(5)) + 4.0*(q(5)*q(5)))

	P2 = ((2.0*q(4) + 5.0*q(5)) - q(6))*C1_6
	b2 = q(4)*((10.0*q(4) - 31.0*q(5)) + 11.0*q(6)) + (q(5)*(25.0*q(5) - 19.0*q(6)) + 4.0*(q(6)*q(6)))

	disc = (abs(b2-b0) >= min(b0, b1, b2) .or. (abs(maxval(cfl) - minval(cfl)) > 0.0)) ! check for discontinuity

	! Nonlinear weights
	tau = abs(b2 - b0)
	w0 = d0*weight_fac(tau, b0)
	w1 = d1*weight_fac(tau, b1)
	w2 = d2*weight_fac(tau, b2)

	wnorm = 1.0 / (w0 + w1 + w2)
	wpl = (w0*P0 + w1*P1 + w2*P2) * wnorm

	! MP limiter (He et al. 2016)
	alpha = 2.0
	qul = q(4) + alpha*(q(4) - q(3))
	qmp = q(4) + minmod2((q(5)-q(4)), (qul-q(4)))

	dm2 = q(3) - 2.0*q(2) + q(1)
	dm1 = q(2) - 2.0*q(3) + q(4)
	dd0 = q(5) - 2.0*q(4) + q(3)
	dd1 = q(4) - 2.0*q(5) + q(6)
	dd2 = q(5) - 2.0*q(6) + q(7)

	dm4p = minmod6( (4.0*dd0 - dd1), (4.0*dd1 - dd0), dd0, dd1, dm1, dd2 )
	dm4m = minmod6( (4.0*dm1 - dd0), (4.0*dd0 - dm1), dm1, dd0, dm2, dd1 )
	qmd = 0.5*((q(5) + q(4)) - dm4p)
	qlc = 0.5*((qul + q(4)) + alpha*dm4m)

	qmin = max(min(q(4), q(5), qmd), min(q(4), qul, qlc))
	qmax = min(max(q(4), q(5), qmd), max(q(4), qul, qlc))
	q0_min = min(q(4), qmp) ; q0_max = max(q(4), qmp)

	if (((wpl-qmin)*(wpl-qmax) > 0.0) .or. (((qmax-qmin) > (q0_max-q0_min)) .and. disc )) then
		lim = .true.
	endif

	! Right state at i-1/2
	P0 = ((2.0*q(6) - 7.0*q(5)) + 11.0*q(4))*C1_6
	b0 = q(6)*((4.0*q(6) - 19.0*q(5)) + 11.0*q(4)) + &
			(q(5)*(25.0*q(5) - 31.0*q(4)) + 10.0*(q(4)*q(4)))

	P1 = ((-q(5) + 5.0*q(4)) + 2.0*q(3))*C1_6
	b1 = q(5)*((4.0*q(5) - 13.0*q(4)) + 5.0*q(3)) + &
			(q(4)*(13.0*q(4) - 13.0*q(3)) + 4.0*(q(3)*q(3)))

	P2 = ((2.0*q(4) + 5.0*q(3)) - q(2))*C1_6
	b2 = q(4)*((10.0*q(4) - 31.0*q(3)) + 11.0*q(2)) + &
			(q(3)*(25.0*q(3) - 19.0*q(2)) + 4.0*(q(2)*q(2)))

	disc = (abs(b2-b0) >= min(b0, b1, b2) &
				.or. (abs(maxval(cfl) - minval(cfl)) > 0.0)) ! check for discontinuity

	! Nonlinear weights
	tau = abs(b2 - b0)
	w0 = d0*weight_fac(tau, b0)
	w1 = d1*weight_fac(tau, b1)
	w2 = d2*weight_fac(tau, b2)

	wnorm = 1.0 / (w0 + w1 + w2)
	wmr = (w0*P0 + w1*P1 + w2*P2) * wnorm

	! MP limiter (He et al. 2016)
	qul = q(4) + alpha*(q(4) - q(5))
	qmp = q(4) + minmod2((q(3)-q(4)), (qul-q(4)))

	dm2 = q(5) - 2.0*q(6) + q(7)
	dm1 = q(6) - 2.0*q(5) + q(4)
	dd0 = q(3) - 2.0*q(4) + q(5)
	dd1 = q(4) - 2.0*q(3) + q(2)
	dd2 = q(3) - 2.0*q(2) + q(1)

	dm4p = minmod6( (4.0*dd0 - dd1), (4.0*dd1 - dd0), dd0, dd1, dm1, dd2 )
	dm4m = minmod6( (4.0*dm1 - dd0), (4.0*dd0 - dm1), dm1, dd0, dm2, dd1 )

	qmd = 0.5*((q(3) + q(4)) - dm4p)
	qlc = 0.5*((qul + q(4)) + alpha*dm4m)

	qmin = max(min(q(4), q(3), qmd), min(q(4), qul, qlc))
	qmax = min(max(q(4), q(3), qmd), max(q(4), qul, qlc))
	q0_min = min(q(4), qmp) ; q0_max = max(q(4), qmp)

	if ( (maxval(cfl) > 0.4) .or. ((wmr-qmin)*(wmr-qmax) > 0.0) .or. &
			(((qmax-qmin) > (q0_max-q0_min)) .and. disc )) then
		lim = .true.
	endif

	if (lim) then
		wpl = ((-q(3) + 5.0*q(4)) + 2.0*q(5))*C1_6
		wpl = min(max(q(4), q(5)), wpl) ; wpl = max(min(q(4), q(5)), wpl)
		wmr = ((-q(5) + 5.0*q(4)) + 2.0*q(3))*C1_6
		wmr = min(max(q(4), q(3)), wmr) ; wmr = max(min(q(4), q(3)), wmr)
		dA = wpl - wmr ; mA = 0.5*( wpl + wmr )
		if ((q(5)-q(4))*(q(4)-q(3)) <= 0.) then
			wmr = q(4) ; wpl = q(4)
		elseif ( dA*(q(4)-mA) > (dA*dA)/6. ) then
			wmr = (3.*q(4)) - 2.*wpl
		elseif ( dA*(q(4)-mA) < - (dA*dA)/6. ) then
			wpl = (3.*q(4)) - 2.*wmr
		endif
	endif

	a6 = 6.*q(4) - 3. * (wpl + wmr) ! Curvature
	if (u >= 0.0) then
			wq = (wpl - 0.5 * cfl(2) * ((wpl - wmr) - a6 * (1. - 2./3. * cfl(2))))
	else
			wq = (wmr + 0.5 * cfl(2) * ((wpl - wmr) + a6 * (1. - 2./3. * cfl(2))))
	endif

end subroutine weno5_reconstruction

!> 7th-order weno z-type reconstruction flux
pure subroutine weno7_reconstruction(wq, q, u, cfl)
	real, intent(in) :: q(7)  !< tracer concentration from i-3 to i+3 [conc]
	real, intent(in) :: u     !< advective flux [H L2 ~> m3 or kg]
	real, intent(in) :: cfl(3)!< absolute value of the advective upwind-cell CFL number [nondim]
	real, intent(out) :: wq   !< weno flux  [conc]

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
	real :: dm2, dd2, delta, Dq6, cfl_diff, wpl, wmr, dA, mA, a6
	logical :: lim, disc

	lim = .false.

	! Left state at i+1/2
	P0 = (((-3.0*q(1) + 13.0*q(2)) - 23.0*q(3)) + 25.0*q(4))*C1_12
	b0 = q(1)*((547.0*q(1) - 3882.0*q(2)) + (4642.0*q(3) - 1854.0*q(4))) + &
			q(2)*((7043.0*q(2) - 17246.0*q(3)) + 7042.0*q(4)) + &
			q(3)*(11003.0*q(3) - 9402.0*q(4)) + 2107.0*(q(4)*q(4))

	P1 = (((q(2) - 5.0*q(3)) + 13.0*q(4)) + 3.0*q(5))*C1_12
	b1 = q(2)*((267.0*q(2) - 1642.0*q(3)) + (1602.0*q(4) - 494.0*q(5))) + &
			q(3)*((2843.0*q(3) - 5966.0*q(4)) + 1922.0*q(5)) + &
			q(4)*(3443.0*q(4) - 2522.0*q(5)) + 547.0*(q(5)*q(5))

	P2 = (((-q(3) + 7.0*q(4)) + 7.0*q(5)) - q(6))*C1_12
	b2 = q(3)*((547.0*q(3) - 2522.0*q(4)) + (1922.0*q(5) - 494.0*q(6))) + &
			q(4)*((3443.0*q(4) - 5966.0*q(5)) + 1602.0*q(6)) + &
			q(5)*(2843.0*q(5) - 1642.0*q(6)) + 267.0*(q(6)*q(6))

	P3 = (((3.0*q(4) + 13.0*q(5)) - 5.0*q(6)) + q(7))*C1_12
	b3 = q(4)*((2107.0*q(4) - 9402.0*q(5)) + (7042.0*q(6) - 1854.0*q(7))) + &
			q(5)*((11003.0*q(5) - 17246.0*q(6)) + 4642.0*q(7)) + &
			q(6)*(7043.0*q(6) - 3882.0*q(7)) + 547.0*(q(7)*q(7))

	! Nonlinear weights
	tau = abs(b0 - b3)
	w0 = d0*weight_fac(tau, b0)
	w1 = d1*weight_fac(tau, b1)
	w2 = d2*weight_fac(tau, b2)
	w3 = d3*weight_fac(tau, b3)

	wnorm = 1.0 / (w0 + w1 + w2 + w3)
	w0 = w0 * wnorm
	w1 = w1 * wnorm
	w2 = w2 * wnorm
	w3 = w3 * wnorm

	disc = (((abs(w0-d0) + abs(w1-d1) + abs(w2-d2) + abs(w3-d3)) > 1.1) &
						.or. (abs(maxval(cfl) - minval(cfl)) > 0.0)) ! check for discontinuity

	wpl = (w0*P0 + w1*P1 + w2*P2 + w3*P3) !* wnorm

	! MP limiter (He et al. 2016)
	alpha = 2.0
	qul = q(4) + alpha*(q(4) - q(3))
	qmp = q(4) + minmod2((q(5)-q(4)), (qul-q(4)))

	dm2 = q(1) - 2.0*q(2) + q(3)
	dm1 = q(2) - 2.0*q(3) + q(4)
	dd0 = q(5) - 2.0*q(4) + q(3)
	dd1 = q(4) - 2.0*q(5) + q(6)
	dd2 = q(5) - 2.0*q(6) + q(7)

	dm4p = minmod6( (4.0*dd0 - dd1), (4.0*dd1 - dd0), dd0, dd1, dm1, dd2 )
	dm4m = minmod6( (4.0*dm1 - dd0), (4.0*dd0 - dm1), dm1, dd0, dm2, dd1 )

	qmd = 0.5*((q(5) + q(4)) - dm4p)
	qlc = 0.5*((qul + q(4)) + alpha*dm4m)

	qmin = max(min(q(4), q(5), qmd), min(q(4), qul, qlc))
	qmax = min(max(q(4), q(5), qmd), max(q(4), qul, qlc))
	q0_min = min(q(4), qmp) ; q0_max = max(q(4), qmp)

	if ( ((wpl-qmin)*(wpl-qmax) > 0.0) .or. (((qmax-qmin) > (q0_max-q0_min)) .and. disc ) ) then
			lim = .true.
	endif

	! Right state at i-1/2
	P0 = (((-3.0*q(7) + 13.0*q(6)) - 23.0*q(5)) + 25.0*q(4))*C1_12
	b0 = q(7)*((547.0*q(7) - 3882.0*q(6)) + (4642.0*q(5) - 1854.0*q(4))) + &
			q(6)*((7043.0*q(6) - 17246.0*q(5)) + 7042.0*q(4)) + &
			q(5)*(11003.0*q(5) - 9402.0*q(4)) + 2107.0*(q(4)*q(4))

	P1 = (((q(6) - 5.0*q(5)) + 13.0*q(4)) + 3.0*q(3))*C1_12
	b1 = q(6)*((267.0*q(6) - 1642.0*q(5)) + (1602.0*q(4) - 494.0*q(3))) + &
			q(5)*((2843.0*q(5) - 5966.0*q(4)) + 1922.0*q(3)) + &
			q(4)*(3443.0*q(4) - 2522.0*q(3)) + 547.0*(q(3)*q(3))

	P2 = (((-q(5) + 7.0*q(4)) + 7.0*q(3)) - q(2))*C1_12
	b2 = q(5)*((547.0*q(5) - 2522.0*q(4)) + (1922.0*q(3) - 494.0*q(2))) + &
			q(4)*((3443.0*q(4) - 5966.0*q(3)) + 1602.0*q(2)) + &
			q(3)*(2843.0*q(3) - 1642.0*q(2)) + 267.0*(q(2)*q(2))

	P3 = (((3.0*q(4) + 13.0*q(3)) - 5.0*q(2)) + q(1))*C1_12
	b3 = q(4)*((2107.0*q(4) - 9402.0*q(3)) + (7042.0*q(2) - 1854.0*q(1))) + &
			q(3)*((11003.0*q(3) - 17246.0*q(2)) + 4642.0*q(1)) + &
			q(2)*(7043.0*q(2) - 3882.0*q(1)) + 547.0*(q(1)*q(1))

	! Nonlinear weights
	tau = abs(b0 - b3)
	w0 = d0*weight_fac(tau, b0)
	w1 = d1*weight_fac(tau, b1)
	w2 = d2*weight_fac(tau, b2)
	w3 = d3*weight_fac(tau, b3)

	wnorm = 1.0 / (w0 + w1 + w2 + w3)
	w0 = w0 * wnorm
	w1 = w1 * wnorm
	w2 = w2 * wnorm
	w3 = w3 * wnorm

	disc = (((abs(w0-d0) + abs(w1-d1) + abs(w2-d2) + abs(w3-d3)) > 1.1) &
						.or. (abs(maxval(cfl) - minval(cfl)) > 0.0)) ! check for discontinuity
	wmr = (w0*P0 + w1*P1 + w2*P2 + w3*P3) !* wnorm

	! MP limiter (He et al. 2016)
	qul = q(4) + alpha*(q(4) - q(5))
	qmp = q(4) + minmod2((q(3)-q(4)), (qul-q(4)))

	dm2 = q(7) - 2.0*q(6) + q(5)
	dm1 = q(6) - 2.0*q(5) + q(4)
	dd0 = q(3) - 2.0*q(4) + q(5)
	dd1 = q(4) - 2.0*q(3) + q(2)
	dd2 = q(3) - 2.0*q(2) + q(1)

	dm4p = minmod6( (4.0*dd0 - dd1), (4.0*dd1 - dd0), dd0, dd1, dm1, dd2 )
	dm4m = minmod6( (4.0*dm1 - dd0), (4.0*dd0 - dm1), dm1, dd0, dm2, dd1 )
	qmd = 0.5*((q(3) + q(4)) - dm4p)
	qlc = 0.5*((qul + q(4)) + alpha*dm4m)

	qmin = max(min(q(4), q(3), qmd), min(q(4), qul, qlc))
	qmax = min(max(q(4), q(3), qmd), max(q(4), qul, qlc))
	q0_min = min(q(4), qmp) ; q0_max = max(q(4), qmp)

	if ( (maxval(cfl) > 0.4) .or. ((wmr-qmin)*(wmr-qmax) > 0.0) .or. &
			(((qmax-qmin) > (q0_max-q0_min)) .and. disc ) ) then
			lim = .true.
	endif

	if (lim) then
			wpl = ((-q(3) + 5.0*q(4)) + 2.0*q(5)) / 6.0
			wpl = min(max(q(4), q(5)), wpl) ; wpl = max(min(q(4), q(5)), wpl)
			wmr = ((-q(5) + 5.0*q(4)) + 2.0*q(3)) / 6.0
			wmr = min(max(q(4), q(3)), wmr) ; wmr = max(min(q(4), q(3)), wmr)
			dA = wpl - wmr ; mA = 0.5*( wpl + wmr )
			if ((q(5)-q(4))*(q(4)-q(3)) <= 0.) then
				wmr = q(4) ; wpl = q(4)
			elseif ( dA*(q(4)-mA) > (dA*dA)/6. ) then
				wmr = (3.*q(4)) - 2.*wpl
			elseif ( dA*(q(4)-mA) < - (dA*dA)/6. ) then
				wpl = (3.*q(4)) - 2.*wmr
			endif
	endif

	a6 = 6.*q(4) - 3. * (wpl + wmr) ! Curvature
	if (u >= 0.0) then
			wq = (wpl - 0.5 * cfl(2) * ((wpl - wmr) - a6 * (1. - 2./3. * cfl(2))))
	else
			wq = (wmr + 0.5 * cfl(2) * ((wpl - wmr) + a6 * (1. - 2./3. * cfl(2))))
	endif

end subroutine weno7_reconstruction

!> ppm reconstruction flux
pure subroutine PPM_reconstruction(wq_ppm, qm, q0, qp, u, cfl, qext)
	real, intent(in) :: qm, q0, qp !< tracer concentration for 3-stencil wide [conc]
	real, intent(in) :: u        !< advective flux [H L2 ~> m3 or kg]
	real, intent(in) :: cfl      !< absolute value of the advective upwind-cell CFL number [nondim]
	real, intent(in) :: qext     !< check local extrema
	real, intent(out) :: wq_ppm  !< PPM tracer concentration at the cell interface i+1/2  [conc]

	real :: aL, aR, dA, mA, a6 ! local variables

	wq_ppm = 0.0

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

	r = 0.5 * (sign(1.0, a) + sign(1.0, b)) * min(abs(a), abs(b))
end function minmod2

pure elemental function minmod6(a,b,c,d,e,f) result(r)
	real, intent(in) :: a,b,c,d,e,f
	real :: r
	real :: s

	s = sign(1.0, a)

	if ( (sign(1.0,b)==s) .and. (sign(1.0,c)==s) .and. &
			(sign(1.0,d)==s) .and. (sign(1.0,e)==s) .and. &
			(sign(1.0,f)==s) ) then

			r = s * min( abs(a), abs(b), abs(c), abs(d), abs(e), abs(f) )

	else
			r = 0.0
	end if

end function minmod6

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

	! factor = (1.0 + (tau / (b + 1.0e-40))**2)
	factor = 1.e40; if (abs(b) > 1.0e-20*tau) factor = (1.0 + (tau / b))**2

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

!! Example of WENO stencil:
!!
!!               |___________S0__________|
!!               |                       |
!!               |       |___________S1__________|
!!               |       |                       |
!!               |       |       |___________S2__________|
!!             ..|---o---|---o---|---o---|---o---|---o---|...
!!               | I{i-2}| I{i-1}|  I{i} | I{i+1}| I{i+2}|
!!                               |+     -|
!!                             i-1/2    i+1/2
!!
!! WENO stencil: S{i} = [ I{i-2},...,I{i+2} ]

end module MOM_tracer_advect_weno
