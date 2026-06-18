! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> This tracer package is used to test convergence of advection schemes
!! using a smooth, signed (potentially negative) sine wave initial condition.
module convergence_test_tracer

use MOM_coms,            only : EFP_type
use MOM_coupler_types,   only : set_coupler_type_data, atmos_ocn_coupler_flux
use MOM_diag_mediator,   only : diag_ctrl
use MOM_error_handler,   only : MOM_error, FATAL, WARNING
use MOM_file_parser,     only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type,    only : forcing
use MOM_grid,            only : ocean_grid_type
use MOM_hor_index,       only : hor_index_type
use MOM_io,              only : slasher, vardesc, var_desc, query_vardesc
use MOM_open_boundary,   only : ocean_OBC_type
use MOM_restart,         only : query_initialized, set_initialized, MOM_restart_CS
use MOM_spatial_means,   only : global_mass_int_EFP
use MOM_sponge,          only : set_up_sponge_field, sponge_CS
use MOM_time_manager,    only : time_type
use MOM_tracer_registry, only : register_tracer, tracer_registry_type
use MOM_tracer_diabatic, only : tracer_vertdiff, applyTracerBoundaryFluxesInOut
use MOM_unit_scaling,    only : unit_scale_type
use MOM_variables,       only : surface
use MOM_verticalGrid,    only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public register_convergence_test_tracer, initialize_convergence_test_tracer
public convergence_test_tracer_surface_state, convergence_test_tracer_end
public convergence_test_tracer_column_physics, convergence_test_stock

integer, parameter :: NTR = 1  !< The number of tracers in this module.

!> The control structure for the convergence_test_tracer module
type, public :: convergence_test_tracer_CS ; private
  integer :: ntr = NTR                 !< Number of tracers in this module
  logical :: coupled_tracers = .false. !< These tracers are not offered to the coupler.
  character(len=200) :: tracer_IC_file !< The full path to the IC file, or " " to initialize internally.
  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock.
  type(tracer_registry_type), pointer :: tr_Reg => NULL() !< A pointer to the MOM tracer registry
  real, pointer :: tr(:,:,:,:) => NULL() !< The array of tracers used in this subroutine [conc]
  real :: land_val(NTR) = -1.0 !< The value of tr used where land is masked out [conc]
  logical :: use_sponge    !< If true, sponges may be applied somewhere in the domain.
  logical :: tracers_may_reinit !< If true, the tracers may be set up via the initialization code if
  !! they are not found in the restart files.  Otherwise it is a fatal error
  !! if the tracers are not found in the restart files of a restarted run.
  real :: x_origin !< x-coordinate origin for the sine wave [m] or [km] or [degrees_E]
  real :: y_origin !< y-coordinate origin for the sine wave [m] or [km] or [degrees_N]
  real :: x_wavelength !< Wavelength of the sine wave in x [m] or [km] or [degrees]
  real :: y_wavelength !< Wavelength of the sine wave in y [m] or [km] or [degrees]

  integer, dimension(NTR) :: ind_tr !< Indices returned by atmos_ocn_coupler_flux if it is used and
  !! the surface tracer concentrations are to be provided to the coupler.

  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
  !! regulate the timing of diagnostic output.
  type(MOM_restart_CS), pointer :: restart_CSp => NULL() !< A pointer to the restart control structure.

  type(vardesc) :: tr_desc(NTR) !< Descriptions and metadata for the tracers
end type convergence_test_tracer_CS

contains

!> Register tracer fields and subroutines to be used with MOM.
function register_convergence_test_tracer(G, GV, param_file, CS, tr_Reg, restart_CS)
  type(ocean_grid_type),          intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),        intent(in) :: GV   !< The ocean's vertical grid structure
  type(param_file_type),          intent(in) :: param_file !< A structure to parse for run-time parameters
  type(convergence_test_tracer_CS), pointer  :: CS !< The control structure returned by a previous
  !! call to register_convergence_test_tracer.
  type(tracer_registry_type),     pointer    :: tr_Reg !< A pointer that is set to point to the control
  !! structure for the tracer advection and diffusion module
  type(MOM_restart_CS), target, intent(inout) :: restart_CS !< MOM restart control struct

  ! Local variables
  character(len=80)  :: name, longname
# include "version_variable.h"
  character(len=40)  :: mdl = "convergence_test_tracer"
  character(len=200) :: inputdir
  character(len=48)  :: flux_units
  real, pointer :: tr_ptr(:,:,:) => NULL()
  logical :: register_convergence_test_tracer
  integer :: isd, ied, jsd, jed, nz, m
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed ; nz = GV%ke

  if (associated(CS)) then
      call MOM_error(FATAL, "register_convergence_test_tracer called with an "// &
        "associated control structure.")
  endif
  allocate(CS)

  call log_version(param_file, mdl, version, "")

  call get_param(param_file, mdl, "CONVERGENCE_TEST_X_ORIGIN", CS%x_origin, &
      "The x-coordinate of the phase origin of the sine wave.", &
      units=G%x_ax_unit_short, default=0.)
  call get_param(param_file, mdl, "CONVERGENCE_TEST_Y_ORIGIN", CS%y_origin, &
      "The y-coordinate of the phase origin of the sine wave.", &
      units=G%y_ax_unit_short, default=0.)
  call get_param(param_file, mdl, "CONVERGENCE_TEST_X_WAVELENGTH", CS%x_wavelength, &
      "The wavelength of the sine wave in x.", &
      units=G%x_ax_unit_short, default=1.)
  call get_param(param_file, mdl, "CONVERGENCE_TEST_Y_WAVELENGTH", CS%y_wavelength, &
      "The wavelength of the sine wave in y.", &
      units=G%y_ax_unit_short, default=1.)
  call get_param(param_file, mdl, "CONVERGENCE_TEST_TRACER_IC_FILE", CS%tracer_IC_file, &
      "The name of a file from which to read the initial "//&
      "conditions for the tracers, or blank to initialize "//&
      "them internally.", default=" ")

  if (len_trim(CS%tracer_IC_file) >= 1) then
      call get_param(param_file, mdl, "INPUTDIR", inputdir, default=".")
      CS%tracer_IC_file = trim(slasher(inputdir))//trim(CS%tracer_IC_file)
      call log_param(param_file, mdl, "INPUTDIR/CONVERGENCE_TEST_TRACER_IC_FILE", &
        CS%tracer_IC_file)
  endif
  call get_param(param_file, mdl, "SPONGE", CS%use_sponge, &
      "If true, sponges may be applied anywhere in the domain. "//&
      "The exact location and properties of those sponges are "//&
      "specified from MOM_initialization.F90.", default=.false.)
  call get_param(param_file, mdl, "TRACERS_MAY_REINIT", CS%tracers_may_reinit, &
      "If true, tracers may go through the initialization code "//&
      "if they are not found in the restart files.  Otherwise "//&
      "it is a fatal error if the tracers are not found in the "//&
      "restart files of a restarted run.", default=.false.)

  allocate(CS%tr(isd:ied,jsd:jed,nz,NTR), source=0.0)

  do m=1,NTR
      write(name,'("tr",I0)') m
      write(longname,'("Concentration of Tracer ",I2.2)') m
      CS%tr_desc(m) = var_desc(name, units="kg kg-1", longname=longname, caller=mdl)
      if (GV%Boussinesq) then ; flux_units = "kg kg-1 m3 s-1"
      else ; flux_units = "kg s-1" ; endif

      tr_ptr => CS%tr(:,:,:,m)
      call register_tracer(tr_ptr, tr_Reg, param_file, G%HI, GV, &
        name=name, longname=longname, units="kg kg-1", &
        registry_diags=.true., flux_units=flux_units, &
        restart_CS=restart_CS, mandatory=.not.CS%tracers_may_reinit, &
        nonneg_lim=.false.)

      if (CS%coupled_tracers) &
        CS%ind_tr(m) = atmos_ocn_coupler_flux(trim(name)//'_flux', &
        flux_type=' ', implementation=' ', caller="register_convergence_test_tracer")
  enddo

  CS%tr_Reg => tr_Reg
  CS%restart_CSp => restart_CS
  register_convergence_test_tracer = .true.
end function register_convergence_test_tracer

!> Initializes the convergence test tracer with a sine wave IC.
subroutine initialize_convergence_test_tracer(restart, day, G, GV, h, diag, OBC, CS, &
  sponge_CSp)
  logical,                              intent(in) :: restart !< .true. if the fields have already
  !! been read from a restart file.
  type(time_type),              target, intent(in) :: day  !< Time of the start of the run.
  type(ocean_grid_type),                intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),              intent(in) :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
      intent(in) :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(diag_ctrl),              target, intent(in) :: diag !< A structure that is used to regulate
  !! diagnostic output.
  type(ocean_OBC_type),                 pointer    :: OBC  !< This open boundary condition type specifies
  !! whether, where, and what open boundary conditions are used.
  type(convergence_test_tracer_CS),     pointer    :: CS !< The control structure returned by a previous
  !! call to register_convergence_test_tracer.
  type(sponge_CS),                      pointer    :: sponge_CSp !< Pointer to the control structure for the sponges.

  ! Local variables
  character(len=16) :: name
  real :: locx, locy
  real :: PI        ! 3.14159... [nondim]
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz, m
  integer :: IsdB, IedB, JsdB, JedB

  PI = 4.0*atan(1.0)

  if (.not.associated(CS)) return
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  CS%diag => diag
  CS%ntr = NTR
  do m=1,NTR
    call query_vardesc(CS%tr_desc(m), name=name, &
      caller="initialize_convergence_test_tracer")
    if ((.not.restart) .or. (CS%tracers_may_reinit .and. .not. &
      query_initialized(CS%tr(:,:,:,m), name, CS%restart_CSp))) then
      do k=1,nz ; do j=js,je ; do i=is,ie
          CS%tr(i,j,k,m) = 0.0
      enddo ; enddo ; enddo

      k=1
      do j=js,je ; do i=is,ie
        locx = (G%geoLonT(i,j) - CS%x_origin) / CS%x_wavelength
        locy = (G%geoLatT(i,j) - CS%y_origin) / CS%y_wavelength
        CS%tr(i,j,k,m) = sin(acos(0.0)*2.*locx) !*sin(acos(0.0)*2.*locy)
      enddo ; enddo

      ! k=1 ! Cut cylinder
      ! do j=js,je ; do i=is,ie
      !  locx = (G%geoLonT(i,j)-CS%x_origin)/CS%x_width
      !  locy = (G%geoLatT(i,j)-CS%y_origin)/CS%y_width
      !  if ((locx**2) + (locy**2) <= 1.0) CS%tr(i,j,k,m) = 1.0
      ! !  if (abs(locx) < 0.2 .and. locy < 0.7) CS%tr(i,j,k,m) = 0.0
      ! !  if (locx>0.0 .and. abs(locy)<0.2) CS%tr(i,j,k,m) = 0.0
      ! enddo ; enddo

      call set_initialized(CS%tr(:,:,:,m), name, CS%restart_CSp)
    endif
  enddo

end subroutine initialize_convergence_test_tracer

!> Applies diapycnal diffusion to the convergence test tracers.
subroutine convergence_test_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, G, GV, US, CS, &
  evap_CFL_limit, minimum_forcing_depth)
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h_old !< Layer thickness before entrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h_new !< Layer thickness after entrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: ea   !< Fluid entrained from layer above [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: eb   !< Fluid entrained from layer below [H ~> m or kg m-2].
  type(forcing),           intent(in) :: fluxes !< Thermodynamic and tracer forcing fields.
  real,                    intent(in) :: dt   !< Time step [T ~> s]
  type(unit_scale_type),   intent(in) :: US   !< A dimensional unit scaling type
  type(convergence_test_tracer_CS), pointer :: CS !< The control structure returned by a previous
  !! call to register_convergence_test_tracer.
  real, optional,          intent(in) :: evap_CFL_limit !< Limit on fraction of water fluxed out of top layer [nondim]
  real, optional,          intent(in) :: minimum_forcing_depth !< Smallest depth over which fluxes can be applied [H ~> m or kg m-2]

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: h_work
  integer :: i, j, k, is, ie, js, je, nz, m
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.associated(CS)) return

  if (present(evap_CFL_limit) .and. present(minimum_forcing_depth)) then
    do m=1,CS%ntr
      do k=1,nz ; do j=js,je ; do i=is,ie
        h_work(i,j,k) = h_old(i,j,k)
      enddo ; enddo ; enddo
      call applyTracerBoundaryFluxesInOut(G, GV, CS%tr(:,:,:,m), dt, fluxes, h_work, &
        evap_CFL_limit, minimum_forcing_depth)
      call tracer_vertdiff(h_work, ea, eb, dt, CS%tr(:,:,:,m), G, GV)
      enddo
  else
    do m=1,NTR
      call tracer_vertdiff(h_old, ea, eb, dt, CS%tr(:,:,:,m), G, GV)
    enddo
  endif

end subroutine convergence_test_tracer_column_physics

!> Extracts surface fields from this tracer package for the coupler.
subroutine convergence_test_tracer_surface_state(sfc_state, h, G, GV, CS)
  type(ocean_grid_type),   intent(in)    :: G  !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV !< The ocean's vertical grid structure
  type(surface),           intent(inout) :: sfc_state !< Surface state of the ocean.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h  !< Layer thickness [H ~> m or kg m-2].
  type(convergence_test_tracer_CS), pointer :: CS !< The control structure returned by a previous
  !! call to register_convergence_test_tracer.

  integer :: m, is, ie, js, je, isd, ied, jsd, jed
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  if (.not.associated(CS)) return

  if (CS%coupled_tracers) then
    do m=1,CS%ntr
      call set_coupler_type_data(CS%tr(:,:,1,m), CS%ind_tr(m), sfc_state%tr_fields, &
          idim=(/isd, is, ie, ied/), jdim=(/jsd, js, je, jed/), turns=G%HI%turns)
    enddo
  endif

end subroutine convergence_test_tracer_surface_state

!> Calculate the mass-weighted integral of all tracer stocks.
function convergence_test_stock(h, stocks, G, GV, CS, names, units, stock_index)
  type(ocean_grid_type),            intent(in)    :: G      !< The ocean's grid structure
  type(verticalGrid_type),          intent(in)    :: GV     !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in) :: h   !< Layer thicknesses [H ~> m or kg m-2]
  type(EFP_type), dimension(:),     intent(out)   :: stocks !< Mass-weighted integrated tracer amount [kg conc].
  type(convergence_test_tracer_CS), pointer       :: CS     !< The control structure returned by a previous
  !! call to register_convergence_test_tracer.
  character(len=*), dimension(:),   intent(out)   :: names  !< Names of the stocks calculated.
  character(len=*), dimension(:),   intent(out)   :: units  !< Units of the stocks calculated.
  integer, optional,                intent(in)    :: stock_index !< Coded index of a specific stock being sought.
  integer                                         :: convergence_test_stock !< Number of stocks calculated here.

  ! Local variables
  integer :: is, ie, js, je, nz, m
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  convergence_test_stock = 0
  if (.not.associated(CS)) return
  if (CS%ntr < 1) return

  if (present(stock_index)) then ; if (stock_index > 0) then
    return
  endif ; endif

  do m=1,CS%ntr
    call query_vardesc(CS%tr_desc(m), name=names(m), units=units(m), &
      caller="convergence_test_stock")
    stocks(m) = global_mass_int_EFP(h, G, GV, CS%tr(:,:,:,m), on_PE_only=.true.)
  enddo
  convergence_test_stock = CS%ntr

end function convergence_test_stock

!> Deallocate memory associated with this module
subroutine convergence_test_tracer_end(CS)
  type(convergence_test_tracer_CS), pointer :: CS !< The control structure returned by a previous
  !! call to register_convergence_test_tracer.
  if (associated(CS)) then
    if (associated(CS%tr)) deallocate(CS%tr)
    deallocate(CS)
  endif
end subroutine convergence_test_tracer_end

end module convergence_test_tracer
