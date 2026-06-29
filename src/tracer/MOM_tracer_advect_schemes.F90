! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!>  This module contains constants for the tracer advection schemes.
module MOM_tracer_advect_schemes

use MOM_error_handler,   only : MOM_error, FATAL

implicit none ; public

! The following are public parameter constants
integer, parameter :: ADVECT_PLM        = 0 !< PLM advection scheme
integer, parameter :: ADVECT_PPMH3      = 1 !< PPM:H3 advection scheme
integer, parameter :: ADVECT_PPM        = 2 !< PPM advection scheme
integer, parameter :: ADVECT_WENO5      = 3 !< WENO5 advection scheme
integer, parameter :: ADVECT_WENO7      = 4 !< WENO7 advection scheme
integer, parameter :: ADVECT_PPMWENO5  = 5 !< PPM with WENO5 edges

!> Documentation for tracer advection schemes
character(len=*), parameter :: TracerAdvectionSchemeDoc = &
  "  PLM    - Piecewise Linear Method\n"//&
  "  PPM    - Piecewise Parabolic Method (Colella-Woodward)\n"//&
  "  PPM:H3 - Piecewise Parabolic Method (Huyhn 3rd order)\n"// &
  "  WENO5  - Weighted Essentially Non-Oscillatory, 5th order\n"//&
  "  WENO7  - Weighted Essentially Non-Oscillatory, 7th order\n"//&
  "  PPM:W5  - Piecewise Parabolic Method (WENO 5th order for the edges)"

contains

!> Numeric value of tracer_advect_scheme corresponding to scheme name
subroutine set_tracer_advect_scheme(scheme_value, advect_scheme_name)
  character(len=*), intent(in) :: advect_scheme_name !< Name of the advection scheme
  integer,         intent(out) :: scheme_value       !< Integer value of the advection scheme

  select case (trim(advect_scheme_name))
    case ("")
      scheme_value = -1
    case ("PLM")
      scheme_value = ADVECT_PLM
    case ("PPM:H3")
      scheme_value = ADVECT_PPMH3
    case ("PPM")
      scheme_value = ADVECT_PPM
    case ("WENO5")
      scheme_value = ADVECT_WENO5
    case ("WENO7")
      scheme_value = ADVECT_WENO7
    case ("PPM:W5")
      scheme_value = ADVECT_PPMWENO5
    case default
      call MOM_error(FATAL, "set_tracer_advect_scheme: "//&
        "Unknown TRACER_ADVECTION_SCHEME = "//trim(advect_scheme_name))
  end select
end subroutine set_tracer_advect_scheme

end module MOM_tracer_advect_schemes
