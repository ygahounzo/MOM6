! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Initial conditions for the idealized baroclinic eddies test case
module baroclinic_eddies_initialization

use MOM_grid, only : ocean_grid_type
use MOM_unit_scaling, only : unit_scale_type
use MOM_verticalGrid, only : verticalGrid_type
use MOM_file_parser, only : param_file_type

implicit none ; private

#include <MOM_memory.h>

! Private (module-wise) parameters
character(len=40) :: mdl = "baroclinic_eddies_initialization" !< This module's name.

public baroclinic_eddies_init_temperature_salinity

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

contains

!> Returns -1, 0, or 1 according to the sign of its argument, matching the behavior of
!! numpy's sign() function.  This is used to build the smoothed step/indicator functions
!! (s_lo, s_hi, and the zonal box function) used to blend the primary front and the
!! localized crest perturbation.
real function sgn(x)
  real, intent(in) :: x !< The value whose sign is returned [arbitrary]
  if (x > 0.0) then
    sgn = 1.0
  elseif (x < 0.0) then
    sgn = -1.0
  else
    sgn = 0.0
  endif
end function sgn

!> Initialization of temperature and salinity for the idealized baroclinic eddies test case of
!! Ilicak et al. (2012) and Petersen et al. (2015).  The domain is a periodic channel that is
!! linearly stratified in the vertical (warmer at the surface), with the northern half of the
!! domain warmer than the southern half.  The two water masses are separated by a narrow
!! meridional transition band whose position wiggles sinusoidally with longitude, and which
!! carries an additional small, spatially localized perturbation over one wave crest to break
!! the exact zonal periodicity and promote the growth of baroclinic instability.
!!
!! All configuration values for this test case are fixed at their default settings, matching
!! Ilicak et al. (2012) / Petersen et al. (2015); no run-time parameters are read.
!!
!! This routine only sets the initial temperature and salinity.  The companion settings that
!! complete the Ilicak et al. (2012) / Petersen et al. (2015) configuration -- a flat bottom
!! (TOPO_CONFIG = "flat"), uniformly spaced initial layers (THICKNESS_CONFIG = "uniform"), a
!! fluid initially at rest (VELOCITY_CONFIG = "zero"), constant planetary rotation
!! (CORIOLIS_SCHEME = "F_PLANE" with F_0 = 1.2e-4 s-1), and a quadratic bottom drag
!! (CDRAG = 0.01) -- are all set with existing, generic run-time parameters.
subroutine baroclinic_eddies_init_temperature_salinity(T, S, h, depth_tot, G, GV, US, &
                                                       param_file, just_read)
  type(ocean_grid_type),   intent(in)  :: G          !< Grid structure
  type(verticalGrid_type), intent(in)  :: GV         !< The ocean's vertical grid structure.
  type(unit_scale_type),   intent(in)  :: US         !< A dimensional unit scaling type
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: T          !< Potential temperature [C ~> degC]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out) :: S          !< Salinity [S ~> ppt]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)  :: h          !< The model thicknesses [Z ~> m]
  real, dimension(SZI_(G),SZJ_(G)), &
                           intent(in)  :: depth_tot  !< The nominal total depth of the ocean [Z ~> m]
  type(param_file_type),   intent(in)  :: param_file !< A structure indicating the open file
                                                     !! to parse for model parameter values.
                                                     !! Retained for interface compatibility;
                                                     !! this routine reads no run-time parameters.
  logical,                 intent(in)  :: just_read  !< If true, this call will only read
                                                     !! parameters without changing T & S.  Since
                                                     !! this routine reads no parameters, it simply
                                                     !! returns without doing any work.

  ! Local variables set to their fixed default values, in place of the run-time parameters
  ! BC_EDDIES_T_BOT, BC_EDDIES_T_SURF, S_REF, BC_EDDIES_FRONT_TEMP_DELTA, BC_EDDIES_FRONT_WIDTH,
  ! BC_EDDIES_FRONT_WAVE_AMP, BC_EDDIES_FRONT_WAVENUMBER, BC_EDDIES_PERTURB_CREST,
  ! BC_EDDIES_CREST_TEMP_DELTA, BC_EDDIES_CREST_X_START, and BC_EDDIES_CREST_X_END that this
  ! test case would otherwise read via MOM_file_parser.
  real :: T_bot       ! Bottom temperature [C ~> degC]
  real :: T_surf      ! Surface temperature [C ~> degC]
  real :: S_ref       ! The uniform background salinity [S ~> ppt]
  real :: front_temp_delta ! Temperature offset south of the front [C ~> degC]
  real :: front_width  ! Meridional width of the temperature transition band, often
                       ! in [km] or [degrees_N], depending on the value of G%y_axis_units
  real :: front_wave_amp ! Amplitude of the front's meridional displacement, in the same
                       ! units as front_width
  integer :: front_wavenumber ! Number of wavelengths of the front's displacement across the
                       ! zonal extent of the domain [nondim]
  logical :: perturb_crest ! If true, apply the single-crest temperature perturbation
  real :: crest_temp_delta ! Temperature offset associated with the localized crest
                       ! perturbation front, applied only within [crest_x_start, crest_x_end] [C ~> degC]
  real :: crest_x_start, crest_x_end ! Zonal range of the single-crest perturbation, in the
                       ! same units as front_width but along the x-axis

  real :: x, y                ! Positions relative to the domain's western and southern edges,
                               ! in the grid's native axis units, often [km] or [degrees]
  real :: y0                  ! The meridional position of the channel center, in the same
                               ! units as y
  real :: yw                  ! The meridional position of the primary front as a function of x,
                               ! in the same units as y
  real :: yw2                 ! The meridional position of the secondary (crest) front as a
                               ! function of x, in the same units as y
  real :: fy                  ! A cross-front position relative to the primary front, normalized
                               ! by front_width [nondim]
  real :: fy2                 ! A cross-front position relative to the crest front, normalized
                               ! by half of front_width [nondim]
  real :: s_lo, s_hi           ! Smoothed step (0 or 1) indicator functions used to ramp and
                               ! then clamp temp_wave between the front's warm and cold sides [nondim]
  real :: s_lo2, s_hi2         ! As s_lo, s_hi, but for the crest front's temp_wave2 [nondim]
  real :: box                 ! An indicator function (0 or 1) that is 1 only for zonal positions
                               ! within [crest_x_start, crest_x_end] [nondim]
  real :: s_wave2              ! An indicator function (0 or 1) used to blend temp_wave2 in place
                               ! of temp_wave within the crest perturbation region [nondim]
  real :: temp_wave           ! The horizontal temperature perturbation associated with the
                               ! primary front [C ~> degC]
  real :: temp_wave2           ! The horizontal temperature perturbation associated with the
                               ! localized crest front [C ~> degC]
  real :: zc, zi              ! Depths in depth units, positive upward [Z ~> m]
  real :: PI                  ! 3.1415926... calculated as 4*atan(1) [nondim]
  integer :: i, j, k, is, ie, js, je, nz

  if (just_read) return ! There are no run-time parameters to read, so return immediately.

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  ! Default values matching Ilicak et al. (2012) / Petersen et al. (2015).
  T_bot            = 10.1 * US%degC_to_C
  T_surf           = 13.1 * US%degC_to_C
  S_ref            = 35.0 * US%ppt_to_S
  front_temp_delta = -1.2 * US%degC_to_C
  front_width      = 40.0
  front_wave_amp   = 40.0
  front_wavenumber = 3
  perturb_crest    = .true.
  if (perturb_crest) then
    crest_temp_delta = -20.0 * US%degC_to_C
    crest_x_start    = 110.0
    crest_x_end      = 130.0
  else
    crest_temp_delta = 0.0 ; crest_x_start = 0.0 ; crest_x_end = 0.0
  endif

  ! T(:,:,:) = 0.0
  ! S(:,:,:) = 0.0
  ! PI = 4.0*atan(1.0)
  ! y0 = 0.5*G%len_lat

  ! do j=js,je ; do i=is,ie
  !   x = G%geoLonT(i,j) - G%west_lon
  !   y = G%geoLatT(i,j) - G%south_lat

  !   ! The primary front's meridional position wiggles sinusoidally with longitude to seed
  !   ! baroclinic instability (Ilicak et al. 2012; Petersen et al. 2015).
  !   yw = y0 - front_wave_amp * sin((2.0*PI*real(front_wavenumber)*x) / G%len_lon)

  !   fy = 1.0 - ((y - yw) / front_width)
  !   s_lo = 0.5*(sgn(fy) + 1.0)
  !   s_hi = 0.5*(sgn(1.0 - fy) + 1.0)
  !   temp_wave = front_temp_delta * (fy*s_lo*s_hi + (1.0 - s_hi))

  !   if (perturb_crest) then
  !     ! A second, narrower and weaker front, offset from the primary front and centered over
  !     ! one wave crest, is blended in only within the zonal band [crest_x_start, crest_x_end],
  !     ! to break the exact zonal periodicity of the primary front and promote the growth of
  !     ! baroclinic instability.
  !     yw2 = y0 - 0.5*front_wave_amp * sin((PI*(x - crest_x_start)) / (crest_x_end - crest_x_start))

  !     fy2 = 1.0 - ((y - yw2) / (0.5*front_width))
  !     s_lo2 = 0.5*(sgn(fy2) + 1.0)
  !     s_hi2 = 0.5*(sgn(2.0 - fy2) + 1.0)
  !     temp_wave2 = crest_temp_delta * (fy2*s_lo2*s_hi2 + (1.0 - s_hi2))

  !     ! box = 1 only for x in [crest_x_start, crest_x_end], 0 otherwise.
  !     box = 0.5*(sgn(x - crest_x_start) * (-1.0) * sgn(x - crest_x_end) + 1.0)
  !     s_wave2 = box * s_hi2

  !     temp_wave = s_wave2*temp_wave2 + (1.0 - s_wave2)*temp_wave
  !   endif

  !   zi = -depth_tot(i,j)
  !   do k=nz,1,-1
  !     zc = zi + 0.5*h(i,j,k)  ! Position of middle of cell
  !     zi = zi + h(i,j,k)      ! Top interface position
  !     T(i,j,k) = T_bot + (T_surf - T_bot) * ((zc + depth_tot(i,j)) / depth_tot(i,j)) & ! Linear
  !                + temp_wave                                                          ! stratification
  !     S(i,j,k) = S_ref
  !   enddo
  ! enddo ; enddo

  T(:,:,:) = 0.0
  S(:,:,:) = 0.0
  PI = 4.0*atan(1.0)
  y0 = 0.5*G%len_lat

  do j=js,je ; do i=is,ie
    x = G%geoLonT(i,j) - G%west_lon
    y = G%geoLatT(i,j) - G%south_lat

    ! The front's meridional position wiggles sinusoidally with longitude to seed
    ! baroclinic instability (Ilicak et al. 2012; Petersen et al. 2015).
    yw = y0 - front_wave_amp * sin((2.0*PI*real(front_wavenumber)*x) / G%len_lon)

    ! An additional, spatially localized bump over one wave crest breaks the exact zonal
    ! periodicity of the front.  The bump is exactly zero at crest_x_start and crest_x_end,
    ! so it blends continuously into the primary wave with no discontinuity in the front.
    if (perturb_crest .and. (x >= crest_x_start) .and. (x <= crest_x_end)) then
      yw = yw + crest_temp_delta * sin((PI*(x - crest_x_start)) / (crest_x_end - crest_x_start))
    endif

    fy = 1.0 - ((y - yw) / front_width)
    if (fy <= 0.0) then
      temp_wave = 0.0
    elseif (fy >= 1.0) then
      temp_wave = front_temp_delta
    else
      temp_wave = front_temp_delta * fy
    endif

    zi = -depth_tot(i,j)
    do k=nz,1,-1
      zc = zi + 0.5*h(i,j,k)  ! Position of middle of cell
      zi = zi + h(i,j,k)      ! Top interface position
      T(i,j,k) = T_bot + (T_surf - T_bot) * ((zc + depth_tot(i,j)) / depth_tot(i,j)) & ! Linear
                 + temp_wave                                                          ! stratification
      S(i,j,k) = S_ref
    enddo
  enddo ; enddo

end subroutine baroclinic_eddies_init_temperature_salinity

!> \namespace baroclinic_eddies_initialization
!!
!! \section section_baroclinic_eddies Description of the baroclinic eddies test case
!!
!! This test case reproduces the idealized eddying channel configuration of \cite ilicak2012,
!! as further documented by \cite petersen2015, an idealization of the Antarctic Circumpolar
!! Current used to evaluate a model's ability to generate baroclinic eddies.  The domain is a
!! horizontally periodic channel of latitudinal (meridional) extent 500 km and longitudinal
!! (zonal) extent 160 km, with a flat bottom 1000 m deep.  The channel is on an f-plane with a
!! constant Coriolis parameter f = 1.2e-4 s-1.  A quadratic bottom drag with a dimensionless
!! drag coefficient of 0.01 is used to promote baroclinic instability.
!!
!! The domain is initially linearly stratified in the vertical, warmer at the surface.  The
!! northern half of the domain is warmer than the southern half, with the two water masses
!! separated by a narrow meridional transition band whose position wiggles sinusoidally with
!! longitude.  Within a localized zonal band centered on one wave crest, a second, narrower
!! and weaker front (half the width and offset amplitude of the primary front, and with its own
!! temperature offset) is blended in in place of the primary front, breaking the exact zonal
!! periodicity of the front and promoting the growth of baroclinic instability rather than a
!! degenerate, perfectly periodic response.  Note that, as formulated, this crest perturbation
!! has a small temperature discontinuity where it saturates at the edge of its meridional
!! transition band, rather than blending smoothly to zero.
!!
!! A linear equation of state that does not depend on salinity is normally used with this test
!! case (EQN_OF_STATE = "LINEAR", RHO_T0_S0 = 1000, dRHO_dT = -0.2, dRHO_dS = 0.0).
!!
!! This module only sets the initial temperature and salinity fields
!! (baroclinic_eddies_init_temperature_salinity), using fixed default values instead of reading
!! run-time parameters. The remaining initial conditions for this test case are all set with
!! existing, generic MOM6 run-time parameters:
!! TOPO_CONFIG = "flat", THICKNESS_CONFIG = "uniform", VELOCITY_CONFIG = "zero", and
!! CORIOLIS_SCHEME = "F_PLANE" with F_0 = 1.2e-4.

end module baroclinic_eddies_initialization