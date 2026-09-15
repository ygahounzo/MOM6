!>  This module contains the subroutines of the WENO schemes that advect tracers along coordinate surfaces.
module MOM_tracer_advect_weno

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_cpu_clock,       only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,       only : CLOCK_MODULE, CLOCK_ROUTINE
use MOM_diag_mediator,   only : post_data, query_averaging_enabled, diag_ctrl
use MOM_diag_mediator,   only : register_diag_field, safe_alloc_ptr, time_type
use MOM_domains,         only : max_across_PEs
use MOM_domains,         only : create_group_pass, do_group_pass, group_pass_type, pass_var
use MOM_error_handler,   only : MOM_error, FATAL, WARNING
use MOM_grid,            only : ocean_grid_type
use MOM_open_boundary,   only : ocean_OBC_type, OBC_NONE, OBC_DIRECTION_E
use MOM_open_boundary,   only : OBC_DIRECTION_W, OBC_DIRECTION_N, OBC_DIRECTION_S
use MOM_open_boundary,   only : OBC_segment_type
use MOM_tracer_registry, only : tracer_registry_type, tracer_type
use MOM_unit_scaling,    only : unit_scale_type
use MOM_verticalGrid,    only : verticalGrid_type
use MOM_tracer_advect_schemes, only : ADVECT_WENO5, ADVECT_WENO7

implicit none ; private

#include <MOM_memory.h>

public ppmw5_reconstruction
public PPM_reconstruction
public advect_tracer_RK3

!> Persistent control structure for the WENO/RK3 tracer advection scheme.
!! Holds the stage-1/stage-2 provisional tracer and thickness work arrays and their
!! halo-update group-pass objects, so that create_group_pass only needs to be called
!! once (on first use) rather than being rebuilt on every RK3 substep.
type, public :: weno_advect_CS ; private
  logical :: pass_init = .false. !< True once Ts1_s, Ts2_s, hprev_s1, hprev_s2, and their
                                 !! group passes below have been allocated/created.
  real, allocatable :: Ts1_s(:,:,:,:) !< Stage-1 provisional tracer concentration [conc]
  real, allocatable :: Ts2_s(:,:,:,:) !< Stage-2 provisional tracer concentration [conc]
  real, allocatable :: hprev_s1(:,:,:) !< Stage-1 provisional cell volume [H L2 ~> m3 or kg]
  real, allocatable :: hprev_s2(:,:,:) !< Stage-2 provisional cell volume [H L2 ~> m3 or kg]
  type(group_pass_type) :: pass_Ts1_hprev_s1 !< Halo-update group for Ts1_s/hprev_s1
  type(group_pass_type) :: pass_Ts2_hprev_s2 !< Halo-update group for Ts2_s/hprev_s2
end type weno_advect_CS

contains

!> This routine time steps the tracer concentration using the third-order Runge-Kutta (RK3) method.
!! All tracers in Reg must use WENO5 or WENO7.
subroutine advect_tracer_RK3(h_end, uhtr, vhtr, OBC, dt, G, GV, US, weno_CS, Reg, dt_dyn, &
                         default_advect_scheme, id_clock_advect, id_clock_pass, min_thickness, &
                         x_first_in, vol_prev, max_iter_in, update_vol_prev, uhr_out, vhr_out)
  type(ocean_grid_type),   intent(inout) :: G     !< ocean grid structure
  type(verticalGrid_type), intent(in)    :: GV    !< ocean vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h_end !< Layer thickness after advection [H ~> m or kg m-2]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: uhtr  !< Accumulated volume or mass flux through the
                                                  !! zonal faces [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                           intent(in)    :: vhtr  !< Accumulated volume or mass flux through the
                                                  !! meridional faces [H L2 ~> m3 or kg]
  type(ocean_OBC_type),    pointer       :: OBC   !< specifies whether, where, and what OBCs are used
  real,                    intent(in)    :: dt    !< time increment [T ~> s]
  type(unit_scale_type),   intent(in)    :: US    !< A dimensional unit scaling type
  type(weno_advect_CS),    pointer       :: weno_CS !< Persistent control structure for the WENO/RK3
                                                  !! advection scheme
  type(tracer_registry_type), pointer    :: Reg   !< pointer to tracer registry
  real,                    intent(in)    :: dt_dyn !< The baroclinic dynamics time step [T ~> s]
  integer,                 intent(in)    :: default_advect_scheme !< The default tracer advection
                                                  !! scheme to use when a tracer does not specify one
  integer,                 intent(in)    :: id_clock_advect !< CPU clock id for the whole advection step
  integer,                 intent(in)    :: id_clock_pass   !< CPU clock id for halo updates
  real,                    intent(in)    :: min_thickness !< The minimum layer thickness used to
                                                  !! determine whether a cell is "thin" for CFL-limiting
                                                  !! purposes [H ~> m or kg m-2]
  logical,       optional, intent(in)    :: x_first_in !< If present, indicate whether to update
                                                  !! first in the x- or y-direction.
  ! The remaining optional arguments are only used in offline tracer mode.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(inout) :: vol_prev !< Cell volume before advection [H L2 ~> m3 or kg].
                                                  !! If update_vol_prev is true, the returned value is
                                                  !! the cell volume after the transport that was done
                                                  !! by this call, and if all the transport could be
                                                  !! accommodated it should be close to h_end*G%areaT.
  integer,       optional, intent(in)    :: max_iter_in !< The maximum number of iterations
  logical,       optional, intent(in)    :: update_vol_prev !< If present and true, update vol_prev to
                                                  !! return its value after the tracer have been updated.
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), &
                 optional, intent(out)   :: uhr_out !< Remaining accumulated volume or mass fluxes
                                                  !! through the zonal faces [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), &
                 optional, intent(out)   :: vhr_out !< Remaining accumulated volume or mass fluxes
                                                  !! through the meridional faces [H L2 ~> m3 or kg]

  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)) :: &
    hprev           ! cell volume at the end of previous tracer change [H L2 ~> m3 or kg]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: &
    uhr             ! The remaining zonal thickness flux [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: &
    vhr             ! The remaining meridional thickness fluxes [H L2 ~> m3 or kg]
  real :: uh_neglect(SZIB_(G),SZJ_(G)) ! uh_neglect and vh_neglect are the
  real :: vh_neglect(SZI_(G),SZJB_(G)) ! magnitude of remaining transports that
                                       ! can be simply discarded [H L2 ~> m3 or kg].

  real :: Idt                           ! 1/dt [T-1 ~> s-1].
  integer :: max_iter           ! maximum number of iterations in each layer
  integer :: domore_k(SZK_(GV))
  integer :: stencil            ! stencil of the advection scheme
  integer :: nsten_halo         ! number of stencils that fit in the halos
  integer :: i, j, k, m, is, ie, js, je, isd, ied, jsd, jed, nz, itt, ntr
  integer :: isv, iev, jsv, jev ! The valid range of the indices.
  integer :: IsdB, IedB, JsdB, JedB
  integer :: stencil_local          ! Stencil for the local adection scheme
  integer :: local_advect_scheme(Reg%ntr) ! contains the list of the advection for each tracer
  real :: CFL_max_global  !< global max outflow CFL used to set max_iter [nondim]
  real :: CFL_face        !< per-cell outflow CFL scratch [nondim]
  real, parameter :: CFL_subcycle = 0.4  !< per-subcycle outflow CFL limit, passed through to
                                         !! rk3_substep as the single source of truth [nondim]
  logical :: domore_j(SZJ_(G), SZK_(GV))
  logical :: dump_cfl  !< True on the first subcycle only, to diagnose CFL_scalar_x/CFL_scalar_y
  type(group_pass_type) :: pass_group

  if (.not. associated(Reg)) call MOM_error(FATAL, "MOM_tracer_advect_RK3: "// &
       "register_tracer must be called before advect_tracer.")
  if (Reg%ntr==0) return
  call cpu_clock_begin(id_clock_advect)

  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB
  ntr = Reg%ntr
  Idt = 1.0 / dt

  stencil = 2

  do m=1,ntr
    local_advect_scheme(m) = Reg%Tr(m)%advect_scheme
    if (local_advect_scheme(m) < 0) local_advect_scheme(m) = default_advect_scheme
    if (local_advect_scheme(m) == ADVECT_WENO5) then
      stencil_local = 3
    elseif (local_advect_scheme(m) == ADVECT_WENO7) then
      stencil_local = 4
    else
      call MOM_error(FATAL, "advect_tracer_rk3: all tracers must use WENO5 or WENO7.")
    endif
    stencil = max(stencil, stencil_local)
  enddo

  if (min(is-isd, ied-ie, js-jsd, jed-je) < stencil) &
    call MOM_error(FATAL, "advect_tracer_rk3: stencil wider than halo.")

  max_iter = 2*max(1, INT(CEILING(dt/dt_dyn)))

  ! Set up group pass: uhr, vhr, hprev, and all tracer fields.
  call cpu_clock_begin(id_clock_pass)
  call create_group_pass(pass_group, uhr, vhr, G%Domain)
  call create_group_pass(pass_group, hprev, G%Domain)
  do m=1,ntr
    call create_group_pass(pass_group, Reg%Tr(m)%t, G%Domain)
  enddo
  call cpu_clock_end(id_clock_pass)

  ! Halo rows are never active; initialize once so face-index edge checks are safe.
  domore_j(:,:) = .false.

  !$OMP parallel default(shared)
  !$OMP do
  do k=1,nz
    do j=jsd,jed ; do I=IsdB,IedB ; uhr(I,j,k) = 0.0 ; enddo ; enddo
    do J=JsdB,JedB ; do i=isd,ied ; vhr(i,J,k) = 0.0 ; enddo ; enddo
    do j=jsd,jed ; do i=isd,ied ; hprev(i,j,k) = 0.0 ; enddo ; enddo
    !  Put the remaining (total) thickness fluxes into uhr and vhr.
    do j=js,je ; do I=is-1,ie ; uhr(I,j,k) = uhtr(I,j,k) ; enddo ; enddo
    do J=js-1,je ; do i=is,ie ; vhr(i,J,k) = vhtr(i,J,k) ; enddo ; enddo
    if (.not. present(vol_prev)) then
      !   This loop reconstructs the thickness field the last time that the
      ! tracers were updated, probably just after the diabatic forcing.  A useful
      ! diagnostic could be to compare this reconstruction with that older value.
      do j=js,je ; do i=is,ie
        hprev(i,j,k) = max(0.0, G%areaT(i,j)*h_end(i,j,k) + &
            ((uhtr(I,j,k) - uhtr(I-1,j,k)) + (vhtr(i,J,k) - vhtr(i,J-1,k))))
      ! In the case that the layer is now dramatically thinner than it was previously,
      ! add a bit of mass to avoid truncation errors.  This will lead to
      ! non-conservation of tracers
        hprev(i,j,k) = hprev(i,j,k) + &
            max(0.0, 1.0e-13*hprev(i,j,k) - G%areaT(i,j)*h_end(i,j,k))
      enddo ; enddo
    else
      do j=js,je ; do i=is,ie
        hprev(i,j,k) = vol_prev(i,j,k)
      enddo ; enddo
    endif
  enddo
  !$OMP do
  do j=jsd,jed ; do I=isd,ied-1
    uh_neglect(I,j) = GV%H_subroundoff * MIN(G%areaT(i,j), G%areaT(i+1,j))
  enddo ; enddo
  !$OMP do
  do J=jsd,jed-1 ; do i=isd,ied
    vh_neglect(i,J) = GV%H_subroundoff * MIN(G%areaT(i,j), G%areaT(i,j+1))
  enddo ; enddo
  !$OMP do
  do m=1,ntr
    if (associated(Reg%Tr(m)%ad_x)) Reg%Tr(m)%ad_x(:,:,:) = 0.0
    if (associated(Reg%Tr(m)%ad_y)) Reg%Tr(m)%ad_y(:,:,:) = 0.0
    if (associated(Reg%Tr(m)%advection_xy)) Reg%Tr(m)%advection_xy(:,:,:) = 0.0
    if (associated(Reg%Tr(m)%ad2d_x)) Reg%Tr(m)%ad2d_x(:,:) = 0.0
    if (associated(Reg%Tr(m)%ad2d_y)) Reg%Tr(m)%ad2d_y(:,:) = 0.0
    if (associated(Reg%Tr(1)%cfl_x)) Reg%Tr(1)%cfl_x(:,:,:) = 0.0
    if (associated(Reg%Tr(1)%cfl_y)) Reg%Tr(1)%cfl_y(:,:,:) = 0.0
  enddo
  !$OMP end parallel

  ! Pre-compute the exact number of subcycles from the global max outflow CFL.
  CFL_max_global = 0.0
  do k=1,nz ; do j=js,je ; do i=is,ie
    CFL_face = 0.0
    if (hprev(i,j,k) > G%areaT(i,j)*min_thickness) then
      CFL_face = (max(uhr(I,j,k), 0.0) - min(uhr(I-1,j,k), 0.0) &
                + max(vhr(i,J,k), 0.0) - min(vhr(i,J-1,k), 0.0)) / hprev(i,j,k)
      CFL_max_global = max(CFL_max_global, CFL_face)
    endif
  enddo ; enddo ; enddo
  call max_across_PEs(CFL_max_global)
  max_iter = min(max(ceiling(CFL_max_global / CFL_subcycle), 1), max_iter)

  ! Full domain: fresh halo exchange every iteration makes narrowing unnecessary.
  isv = is ; iev = ie ; jsv = js ; jev = je
  dump_cfl = .true. ! Diagnose CFL_scalar_x/CFL_scalar_y from the first subcycle only.

  do itt=1, max_iter
    ! Exchange uhr, vhr, hprev, and tracers so halos reflect the current residuals.
    call do_group_pass(pass_group, G%Domain, clock=id_clock_pass)

    ! Re-initialize domore_j from current residuals uhr/vhr.
    ! Checking both zonal faces of a row AND the
    ! meridional faces bordering it captures rows that receive inflow from a
    ! CFL-limited neighbour without themselves exceeding CFL.
    !$OMP parallel do default(shared)
    do k=1,nz

      do j=js,je
        domore_j(j,k) = .false.
        do I=is-1,ie
          if (uhr(I,j,k) /= 0.0) then ; domore_j(j,k) = .true. ; exit ; endif
        enddo
        if (.not. domore_j(j,k)) then
          do i=is,ie
            if (vhr(i,j,k) /= 0.0 .or. vhr(i,j-1,k) /= 0.0) then
              domore_j(j,k) = .true. ; exit
            endif
          enddo
        endif
      enddo
      domore_k(k) = 0
      do j=js,je ; if (domore_j(j,k)) then ; domore_k(k) = 1 ; exit ; endif ; enddo
    enddo

    ! SSP-RK3 2D-unsplit step: applies CFL-limited fluxes and subtracts the
    ! consumed portion from uhr/vhr for subsequent iterations.
    call rk3_substep(weno_CS, G, GV, US, OBC, Reg, hprev, uhr, vhr, &
              uh_neglect, vh_neglect, domore_k, domore_j, &
              ntr, nz, isv, iev, jsv, jev, &
              local_advect_scheme, Idt, CFL_subcycle, min_thickness, dump_cfl)

    dump_cfl = .false.

  enddo ! itt

  if (present(uhr_out)) uhr_out(:,:,:) = uhr(:,:,:)
  if (present(vhr_out)) vhr_out(:,:,:) = vhr(:,:,:)
  if (present(vol_prev) .and. present(update_vol_prev)) then
    if (update_vol_prev) vol_prev(:,:,:) = hprev(:,:,:)
  endif

  call cpu_clock_end(id_clock_advect)

end subroutine advect_tracer_RK3

!> One SSP-RK3 2D-unsplit sub-step for use inside advect_tracer_RK3.
subroutine rk3_substep(CS, G, GV, US, OBC, Reg, hprev, uhr, vhr, uh_neglect, vh_neglect, domore_k, &
    domore_j, ntr, nz, isv, iev, jsv, jev, local_advect_scheme, Idt, CFL_subcycle, min_thickness, &
    dump_cfl)
  type(weno_advect_CS),       pointer       :: CS                    !< Persistent WENO/RK3 control structure
  type(ocean_grid_type),      intent(inout) :: G                     !< ocean grid structure
  type(verticalGrid_type),    intent(in)    :: GV                    !< ocean vertical grid structure
  type(unit_scale_type),      intent(in)    :: US                    !< A dimensional unit scaling type
  type(ocean_OBC_type),       pointer       :: OBC             !< specifies whether, where, and what OBCs are used
  type(tracer_registry_type), pointer       :: Reg             !< pointer to tracer registry
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)),  intent(inout) :: hprev !< cell volume at the end of previous
                                                                     !! tracer change [H L2 ~> m3 or kg]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: uhr   !< accumulated volume/mass flux through
                                                                     !! the zonal face [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: vhr   !< accumulated volume/mass flux through
                                                                     !! the meridional face [H L2 ~> m3 or kg]
  real, dimension(SZIB_(G),SZJ_(G)),           intent(in)    :: uh_neglect !< A tiny zonal mass flux that can
                                                                      !! be neglected [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G)),           intent(in)    :: vh_neglect !< A tiny meridional mass flux that can
                                                                  !! be neglected [H L2 ~> m3 or kg]
  integer, dimension(SZK_(GV)),                intent(in)    :: domore_k !< per-layer active flag
  logical, dimension(SZJ_(G),SZK_(GV)),        intent(in)    :: domore_j !< per-row active flag
  integer,                                     intent(in)    :: ntr      !< The number of tracers
  integer,                                     intent(in)    :: nz
  integer,                                     intent(in)    :: isv   !< The starting tracer i-index to work on
  integer,                                     intent(in)    :: iev   !< The ending tracer i-index to work on
  integer,                                     intent(in)    :: jsv   !< The starting tracer j-index to work on
  integer,                                     intent(in)    :: jev   !< The ending tracer j-index to work on
  integer, dimension(ntr),                     intent(in)    :: local_advect_scheme  !< list of advection schemes to use
  real,                                        intent(in)    :: Idt  !< The inverse of dt [T-1 ~> s-1]
  real,                                        intent(in)    :: CFL_subcycle  !< per-subcycle CFL limit [nondim]
  real,                                        intent(in)    :: min_thickness  !< The minimum layer
                                                                               !! thicknesses [H ~> m or kg m-2]
  logical,                                     intent(in)    :: dump_cfl !< If true, write the CFL_scalar_x/CFL_scalar_y
                                                                     !! diagnostics from this substep's combined
                                                                     !! zonal+meridional outflow CFL (should only be
                                                                     !! true on the first subcycle, itt==1)

  real :: uhh(SZIB_(G),SZJ_(G),SZK_(GV))
  real :: vhh(SZI_(G),SZJB_(G),SZK_(GV))
  real :: flux_x(SZIB_(G),SZJ_(G),ntr,nz)
  real :: flux_y(SZI_(G),SZJB_(G),ntr,nz)
  real :: flux_xs(SZIB_(G),SZJ_(G),ntr,nz)
  real :: flux_ys(SZI_(G),SZJB_(G),ntr,nz)
  real :: h_old, h_new, Ihnew_ij, dh_ij, h_neglect
  logical :: do_ij
  integer :: i, j, k, m, n
  real :: CFL_cell_ij                     !< per-cell combined 2D outflow CFL [nondim]
  real :: scale(SZI_(G),SZJ_(G),nz)      !< outflow scale factor: min(1, CFL_subcycle/CFL_cell) [nondim]
  real :: dh(SZI_(G),SZJ_(G),nz)            !< precomputed flux divergence per cell [H]
  logical :: no_flux(SZI_(G),SZJ_(G),nz)    !< .true. if all four face fluxes are zero
  logical :: apply_lim_zs(ntr)  !< per-tracer Zhang-Shu/MPP limiter flag
  type(OBC_segment_type), pointer :: segment => null()
  integer :: m_zs

  h_neglect = GV%H_subroundoff

  if (.not. CS%pass_init) then
    allocate(CS%Ts1_s(SZI_(G),SZJ_(G),SZK_(GV),ntr), source=0.0)
    allocate(CS%Ts2_s(SZI_(G),SZJ_(G),SZK_(GV),ntr), source=0.0)
    allocate(CS%hprev_s1(SZI_(G),SZJ_(G),nz), source=0.0)
    allocate(CS%hprev_s2(SZI_(G),SZJ_(G),nz), source=0.0)
    do m=1,ntr
      call create_group_pass(CS%pass_Ts1_hprev_s1, CS%Ts1_s(:,:,:,m), G%Domain)
    enddo
    call create_group_pass(CS%pass_Ts1_hprev_s1, CS%hprev_s1, G%Domain)
    do m=1,ntr
      call create_group_pass(CS%pass_Ts2_hprev_s2, CS%Ts2_s(:,:,:,m), G%Domain)
    enddo
    call create_group_pass(CS%pass_Ts2_hprev_s2, CS%hprev_s2, G%Domain)
    CS%pass_init = .true.
  endif

  do m_zs=1,ntr ; apply_lim_zs(m_zs) = Reg%Tr(m_zs)%nonneg_lim ; enddo

  flux_x = 0.0 ; flux_y = 0.0

  ! Stage 1: initialize T^n snapshot, compute mass fluxes,
  ! reconstruct stage-1 fluxes, and compute T* and h*.
  do k=1,nz
    ! CS%Ts1_s/CS%Ts2_s/CS%hprev_s1/CS%hprev_s2 persist across calls and are halo-exchanged
    ! in full below, regardless of domore_k.
    do j=G%jsd,G%jed ; do i=G%isd,G%ied
      CS%hprev_s1(i,j,k) = hprev(i,j,k)
      CS%hprev_s2(i,j,k) = hprev(i,j,k)
    enddo ; enddo
    do m=1,ntr
      do j=G%jsd,G%jed ; do i=G%isd,G%ied
        CS%Ts1_s(i,j,k,m) = Reg%Tr(m)%t(i,j,k)
        CS%Ts2_s(i,j,k,m) = Reg%Tr(m)%t(i,j,k)
      enddo ; enddo
    enddo

    if (domore_k(k) > 0) then

      do j=G%jsd,G%jed ; do i=G%isd,G%ied
        scale(i,j,k) = 1.0
      enddo ; enddo

      ! Compute per-cell outflow CFL from uhr/vhr and derive scale in one pass.
      do j=G%jsd,G%jed
        do i=G%isd,G%ied
          CFL_cell_ij = 0.0
          if (hprev(i,j,k) > G%areaT(i,j)*min_thickness) then
            CFL_cell_ij = (max(G%mask2dCu(I,j)*uhr(I,j,k), 0.0) &
                          -min(G%mask2dCu(I-1,j)*uhr(I-1,j,k), 0.0) &
                          +max(G%mask2dCv(i,J)*vhr(i,J,k), 0.0) &
                          -min(G%mask2dCv(i,J-1)*vhr(i,J-1,k), 0.0)) / hprev(i,j,k)
            if (CFL_cell_ij > CFL_subcycle) scale(i,j,k) = CFL_subcycle / CFL_cell_ij
          endif
          ! Diagnose CFL_scalar_x/CFL_scalar_y: face I's value is this cell's own combined 2D
          ! outflow CFL.
          if (dump_cfl) then
            if ((Reg%Tr(1)%id_cflx > 0) .and. &
                (I >= isv-1 .and. I <= iev) .and. (j >= jsv .and. j <= jev)) &
              Reg%Tr(1)%cfl_x(I,j,k) = scale(i,j,k) * CFL_cell_ij
            if ((Reg%Tr(1)%id_cfly > 0) .and. &
                (J >= jsv-1 .and. J <= jev) .and. (i >= isv .and. i <= iev)) &
              Reg%Tr(1)%cfl_y(i,J,k) = scale(i,j,k) * CFL_cell_ij
          endif
        enddo
      enddo

      ! Compute mass fluxes with thin-cell zeroing and CFL scale in a single pass.
      ! uhh/vhh are written exactly once.
      do j=jsv,jev ; do I=isv-1,iev
        if (.not. domore_j(j,k)) then ; uhh(I,j,k) = 0.0 ; cycle ; endif
        if ((uhr(I,j,k) == 0.0) .or. &
            ((uhr(I,j,k) < 0.0) .and. (hprev(i+1,j,k) <= G%areaT(i+1,j)*min_thickness)) .or. &
            ((uhr(I,j,k) > 0.0) .and. (hprev(i,j,k) <= G%areaT(i,j)*min_thickness))) then
          uhh(I,j,k) = 0.0
        elseif (uhr(I,j,k) > 0.0) then
          uhh(I,j,k) = uhr(I,j,k) * scale(I,j,k)
        else
          uhh(I,j,k) = uhr(I,j,k) * scale(I+1,j,k)
        endif
      enddo ; enddo
      do J=jsv-1,jev ; do i=isv,iev
        if (.not. (domore_j(max(J,G%jsd),k) .or. domore_j(min(J+1,G%jed),k))) then
          vhh(i,J,k) = 0.0 ; cycle
        endif
        if ((G%mask2dCv(i,J)*vhr(i,J,k) == 0.0) .or. &
            ((G%mask2dCv(i,J)*vhr(i,J,k) < 0.0) .and. &
            (hprev(i,j+1,k) <= G%areaT(i+1,j)*min_thickness)) .or. &
            ((G%mask2dCv(i,J)*vhr(i,J,k) > 0.0) .and. &
            (hprev(i,j,k) <= G%areaT(i,j)*min_thickness))) then
          vhh(i,J,k) = 0.0
        elseif (G%mask2dCv(i,J)*vhr(i,J,k) > 0.0) then
          vhh(i,J,k) = G%mask2dCv(i,J)*vhr(i,J,k) * scale(i,J,k)
        else
          vhh(i,J,k) = G%mask2dCv(i,J)*vhr(i,J,k) * scale(i,J+1,k)
        endif
      enddo ; enddo
    endif ! mass fluxes; compute_flux_2d is a collective (pass_var) and cannot sit behind per-PE domore_k

    call compute_flux_2d(CS%Ts1_s, uhh, vhh, hprev(:,:,k), OBC, ntr, &
        isv, iev, jsv, jev, k, G, GV, local_advect_scheme, &
        flux_x(:,:,:,k), flux_y(:,:,:,k), domore_j, uhr, vhr, apply_lim_zs, &
        (domore_k(k) > 0))

    if (domore_k(k) > 0) then
      ! Compute T* = (h^n T^n - F1) / h*
      ! h* = h^n - dh  (stored in hprev_s1 for stage-2 reconstruction)
      do j=jsv,jev ; do i=isv,iev
        if (.not. domore_j(j,k)) cycle
        no_flux(i,j,k) = (uhh(I,j,k) == 0.0) .and. (uhh(I-1,j,k) == 0.0) .and. &
                          (vhh(i,J,k) == 0.0) .and. (vhh(i,J-1,k) == 0.0)
        dh(i,j,k) = (uhh(I,j,k) - uhh(I-1,j,k)) + (vhh(i,J,k) - vhh(i,J-1,k))
        if (no_flux(i,j,k)) then
          do_ij = .false.
        else
          do_ij = .true.
          h_old = hprev(i,j,k)
          h_new = h_old - dh(i,j,k)
          if (h_new <= 0.0) then
            Ihnew_ij = 0.0 ; do_ij = .false.
          elseif (h_new < h_neglect*G%areaT(i,j)) then
            h_old = h_old + (h_neglect*G%areaT(i,j) - h_new)
            Ihnew_ij = 1.0 / (h_neglect*G%areaT(i,j))
          else
            Ihnew_ij = 1.0 / h_new
          endif
          CS%hprev_s1(i,j,k) = h_new
        endif
        if (associated(OBC)) then
          if ((.not.OBC%exterior_OBC_bug) .and. (OBC%OBC_pe)) then
            if (((OBC%specified_u_BCs_exist_globally .or. OBC%open_u_BCs_exist_globally) .and. &
                  ((OBC%segnum_u(I-1,j) > 0) .or. (OBC%segnum_u(I,j) < 0))) .or. &
                ((OBC%specified_v_BCs_exist_globally .or. OBC%open_v_BCs_exist_globally) .and. &
                  ((OBC%segnum_v(i,J-1) > 0) .or. (OBC%segnum_v(i,J) < 0)))) do_ij = .false.
          endif
        endif
        do m=1,ntr
          if (do_ij) &
            CS%Ts1_s(i,j,k,m) = (Reg%Tr(m)%t(i,j,k)*h_old &
                            - (flux_x(I,j,m,k) - flux_x(I-1,j,m,k)) &
                            - (flux_y(i,J,m,k) - flux_y(i,J-1,m,k))) * Ihnew_ij
        enddo
      enddo ; enddo
    endif ! domore_k
  enddo ! Stage 1 k-loop

  ! T* halo exchange: Ts1_s/hprev_s1's group pass was created once, on first use, above.
  call do_group_pass(CS%pass_Ts1_hprev_s1, G%Domain)

  ! Stage 2: reconstruct from T*, compute T** and h**.
  do k=1,nz
    call compute_flux_2d(CS%Ts1_s, uhh, vhh, CS%hprev_s1(:,:,k), OBC, ntr, &
        isv, iev, jsv, jev, k, G, GV, local_advect_scheme, &
        flux_xs(:,:,:,k), flux_ys(:,:,:,k), domore_j, uhr, vhr, apply_lim_zs, &
        (domore_k(k) > 0))

    if (domore_k(k) > 0) then
      ! Accumulate F1+F2.
      do m=1,ntr
        do j=jsv,jev ; do I=isv-1,iev
          flux_x(I,j,m,k) = flux_x(I,j,m,k) + flux_xs(I,j,m,k)
        enddo ; enddo
        do J=jsv-1,jev ; do i=isv,iev
          flux_y(i,J,m,k) = flux_y(i,J,m,k) + flux_ys(i,J,m,k)
        enddo ; enddo
      enddo

      ! Compute T** = (h^n T^n - (1/4)*(F1+F2)) / h**
      ! T^n is still Reg%Tr(m)%t — it is not updated until the end of stage 3.
      ! h** = h^n - (1/2)*dh  (stored in hprev_s2 for stage-3 reconstruction)
      do j=jsv,jev ; do i=isv,iev
        if (.not. domore_j(j,k)) cycle
        if (no_flux(i,j,k)) then
          do_ij = .false.
        else
          do_ij = .true.
          h_old = hprev(i,j,k)
          h_new = h_old - 0.5*dh(i,j,k)
          if (h_new <= 0.0) then
            Ihnew_ij = 0.0 ; do_ij = .false.
          elseif (h_new < h_neglect*G%areaT(i,j)) then
            h_old = h_old + (h_neglect*G%areaT(i,j) - h_new)
            Ihnew_ij = 1.0 / (h_neglect*G%areaT(i,j))
          else
            Ihnew_ij = 1.0 / h_new
          endif
          CS%hprev_s2(i,j,k) = h_new
        endif
        if (associated(OBC)) then
          if ((.not.OBC%exterior_OBC_bug) .and. (OBC%OBC_pe)) then
            if (((OBC%specified_u_BCs_exist_globally .or. OBC%open_u_BCs_exist_globally) .and. &
                  ((OBC%segnum_u(I-1,j) > 0) .or. (OBC%segnum_u(I,j) < 0))) .or. &
                ((OBC%specified_v_BCs_exist_globally .or. OBC%open_v_BCs_exist_globally) .and. &
                  ((OBC%segnum_v(i,J-1) > 0) .or. (OBC%segnum_v(i,J) < 0)))) do_ij = .false.
          endif
        endif
        do m=1,ntr
          if (do_ij) &
            CS%Ts2_s(i,j,k,m) = (Reg%Tr(m)%t(i,j,k)*h_old &
                            - 0.25*((flux_x(I,j,m,k) - flux_x(I-1,j,m,k)) &
                            +       (flux_y(i,J,m,k) - flux_y(i,J-1,m,k)))) * Ihnew_ij
        enddo
      enddo ; enddo
    endif ! domore_k
  enddo ! Stage 2 k-loop

  ! T** halo exchange: Ts2_s/hprev_s2's group pass was created once, on first use, above.
  call do_group_pass(CS%pass_Ts2_hprev_s2, G%Domain)

  ! Stage 3: reconstruct from T**, combine fluxes, final update.
  do k=1,nz
    call compute_flux_2d(CS%Ts2_s, uhh, vhh, CS%hprev_s2(:,:,k), OBC, ntr, &
        isv, iev, jsv, jev, k, G, GV, local_advect_scheme, &
        flux_xs(:,:,:,k), flux_ys(:,:,:,k), domore_j, uhr, vhr, apply_lim_zs, &
        (domore_k(k) > 0))

    if (domore_k(k) > 0) then
      ! Combine: RK3 tracers get (1/6)*(F1+F2) + (2/3)*F3.
      do m=1,ntr
        do j=jsv,jev ; do I=isv-1,iev
          flux_x(I,j,m,k) = (1.0/6.0)*flux_x(I,j,m,k) + (2.0/3.0)*flux_xs(I,j,m,k)
        enddo ; enddo
        do J=jsv-1,jev ; do i=isv,iev
          flux_y(i,J,m,k) = (1.0/6.0)*flux_y(i,J,m,k) + (2.0/3.0)*flux_ys(i,J,m,k)
        enddo ; enddo
      enddo

      ! Final tracer and thickness update — each stage's flux was already passed through
      ! the Zhang-Shu limiter inside compute_flux_2d, which
      ! guarantees positivity per stage; SSP-RK3's convex combination carries that to T^{n+1}.
      do j=jsv,jev ; do i=isv,iev
        if (.not. domore_j(j,k)) cycle
        if (no_flux(i,j,k)) then
          do_ij = .false.
        else
          do_ij = .true.
          h_old = hprev(i,j,k)
          h_new = h_old - dh(i,j,k)
          if (h_new <= 0.0) then
            Ihnew_ij = 0.0 ; do_ij = .false.
          elseif (h_new < h_neglect*G%areaT(i,j)) then
            h_old = h_old + (h_neglect*G%areaT(i,j) - h_new)
            Ihnew_ij = 1.0 / (h_neglect*G%areaT(i,j))
          else
            Ihnew_ij = 1.0 / h_new
          endif
          hprev(i,j,k) = max(h_new, 0.0)
        end if
        if (associated(OBC)) then
          if ((.not.OBC%exterior_OBC_bug) .and. (OBC%OBC_pe)) then
            if (((OBC%specified_u_BCs_exist_globally .or. OBC%open_u_BCs_exist_globally) .and. &
                  ((OBC%segnum_u(I-1,j) > 0) .or. (OBC%segnum_u(I,j) < 0))) .or. &
                ((OBC%specified_v_BCs_exist_globally .or. OBC%open_v_BCs_exist_globally) .and. &
                  ((OBC%segnum_v(i,J-1) > 0) .or. (OBC%segnum_v(i,J) < 0)))) do_ij = .false.
          endif
        endif
        do m=1,ntr
          if (do_ij) &
            Reg%Tr(m)%t(i,j,k) = (Reg%Tr(m)%t(i,j,k)*h_old &
                - (flux_x(I,j,m,k) - flux_x(I-1,j,m,k)) &
                - (flux_y(i,J,m,k) - flux_y(i,J-1,m,k))) * Ihnew_ij
        enddo
      enddo ; enddo

      ! Consume the applied flux (scaled uhh/vhh) from the remaining uhr/vhr.
      ! Residuals below neglect are zeroed so future iterations don't spin on rounding noise.
      do j=jsv,jev ; do I=isv-1,iev
        if (.not. domore_j(j,k)) cycle
        uhr(I,j,k) = uhr(I,j,k) - uhh(I,j,k)
        if (abs(uhr(I,j,k)) < uh_neglect(I,j)) uhr(I,j,k) = 0.0
      enddo ; enddo
      do J=jsv-1,jev ; do i=isv,iev
        if (.not. (domore_j(max(J,G%jsd),k) .or. domore_j(min(J+1,G%jed),k))) cycle
        vhr(i,J,k) = vhr(i,J,k) - vhh(i,J,k)
        if (abs(vhr(i,J,k)) < vh_neglect(i,J)) vhr(i,J,k) = 0.0
      enddo ; enddo

    endif ! domore_k

    do m=1,ntr
      if (associated(Reg%Tr(m)%ad_x)) then
        do j=jsv,jev ; do I=isv-1,iev
          Reg%Tr(m)%ad_x(I,j,k) = Reg%Tr(m)%ad_x(I,j,k) + flux_x(I,j,m,k)*Idt
        enddo ; enddo
      endif
      if (associated(Reg%Tr(m)%ad_y)) then
        do J=jsv-1,jev ; do i=isv,iev
          Reg%Tr(m)%ad_y(i,J,k) = Reg%Tr(m)%ad_y(i,J,k) + flux_y(i,J,m,k)*Idt
        enddo ; enddo
      endif
      if (associated(Reg%Tr(m)%advection_xy)) then
        do j=jsv,jev ; do i=isv,iev
          Reg%Tr(m)%advection_xy(i,j,k) = Reg%Tr(m)%advection_xy(i,j,k) - &
              ((flux_x(I,j,m,k) - flux_x(I-1,j,m,k)) + &
                (flux_y(i,J,m,k) - flux_y(i,J-1,m,k))) * Idt * G%IareaT(i,j)
        enddo ; enddo
      endif
      if (associated(Reg%Tr(m)%ad2d_x)) then
        do j=jsv,jev ; do I=isv-1,iev
          Reg%Tr(m)%ad2d_x(I,j) = Reg%Tr(m)%ad2d_x(I,j) + flux_x(I,j,m,k)*Idt
        enddo ; enddo
      endif
      if (associated(Reg%Tr(m)%ad2d_y)) then
        do J=jsv-1,jev ; do i=isv,iev
          Reg%Tr(m)%ad2d_y(i,J) = Reg%Tr(m)%ad2d_y(i,J) + flux_y(i,J,m,k)*Idt
        enddo ; enddo
      endif
      if (Reg%Tr(m)%conc_underflow > 0.0) then
        do j=jsv,jev ; do i=isv,iev
          if (abs(Reg%Tr(m)%t(i,j,k)) < Reg%Tr(m)%conc_underflow) &
              Reg%Tr(m)%t(i,j,k) = 0.0
        enddo ; enddo
      endif
    enddo

  enddo ! Stage 3 k-loop

end subroutine rk3_substep

!> Compute zonal and meridional tracer fluxes and apply the Zhang-Shu limiter in one call.
subroutine compute_flux_2d(tk, uhh_in, vhh_in, h_k, OBC, ntr, &
  is, ie, js, je, k, G, GV, advect_schemes, flux_x_out, flux_y_out, domore_j_k, &
  uhr, vhr, apply_lim, do_recon)
  type(ocean_grid_type),                          intent(inout) :: G   !< Ocean grid structure
  type(verticalGrid_type),                        intent(in)    :: GV  !< Ocean vertical grid structure
  integer,                                        intent(in)    :: ntr !< Number of tracers
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV),ntr), intent(in)    :: tk  !< Tracer concentrations [conc]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)),     intent(inout) :: uhh_in   !< Zonal mass flux [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)),     intent(inout) :: vhh_in   !< Meridional mass flux [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJ_(G)),                intent(in)    :: h_k     !< Cell volume at
                                                                            ! current step [H L2 ~> m3 or kg]
  type(ocean_OBC_type),                           pointer       :: OBC  !< Open boundary condition structure
  integer,                                        intent(in)    :: is  !< Start of i-index computational domain
  integer,                                        intent(in)    :: ie  !< End of i-index computational domain
  integer,                                        intent(in)    :: js  !< Start of j-index computational domain
  integer,                                        intent(in)    :: je  !< End of j-index computational domain
  integer,                                        intent(in)    :: k   !< Vertical layer index
  integer, dimension(ntr),                        intent(in)    :: advect_schemes !< Per-tracer advection scheme
                                                                        !! identifier
  real, dimension(SZIB_(G),SZJ_(G),ntr),          intent(out)   :: flux_x_out !< Zonal tracer
                                                                              ! flux [conc H L2 ~> conc m3]
  real, dimension(SZI_(G),SZJB_(G),ntr),          intent(out)   :: flux_y_out !< Meridional tracer
                                                                              ! flux [conc H L2 ~> conc m3]
  logical, dimension(SZJ_(G),SZK_(GV)),           intent(in)    :: domore_j_k !< Per-row per-layer
                                                                              ! activity flag [nondim]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)),     intent(inout) :: uhr   !< accumulated zonal mass flux [H L2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)),     intent(inout) :: vhr   !< accumulated meridional mass flux [H L2]
  logical, dimension(ntr),                        intent(in)    :: apply_lim  !< per-tracer Zhang-Shu enable
  logical,                                        intent(in)    :: do_recon !< If false, skip reconstruction
                                                                          !! (this PE is inactive at k) but still
                                                                          !! participate in pass_var(theta)

  real, dimension(SZI_(G),ntr)          :: Ts2x  !< x-reconstruction stencil (per row)
  real, dimension(SZI_(G),ntr,SZJB_(G)) :: Ts2y  !< y-reconstruction stencil (all rows)
  real, dimension(SZI_(G),SZJ_(G))      :: theta  !< positivity limiter rescaling factor [nondim]
  real :: order3, order5, order7
  real :: T7(7), wq, qext
  real, dimension(SZIB_(G))      :: cfl_row !< CFL at every zonal face along the current j-row
  real, dimension(SZI_(G),SZJB_(G)) :: cfl_v !< CFL at every meridional face in the full column, computed once
  integer :: i, j, m, n, i_up, j_up, ntr_id
  type(OBC_segment_type), pointer :: segment=>NULL()
  real :: hT                  !< tracer mass available in the cell, h_k*max(tk,0) [conc H L2 ~> conc m3]
  real :: outgoing            !< total attempted outgoing flux across all 4 faces [conc H L2 ~> conc m3]

  flux_x_out(:,:,:) = 0.0
  flux_y_out(:,:,:) = 0.0

  if (do_recon) then

    ! Y pre-init: fill Ts2y for all j
    do j=G%jsd,G%jed ; do m=1,ntr ; do i=G%isd,G%ied
      Ts2y(i,m,j) = tk(i,j,k,m)
    enddo ; enddo ; enddo
    if (associated(OBC)) then ; if (OBC%OBC_pe) then
      do n=1,OBC%number_of_segments
        segment=>OBC%segment(n)
        if (.not. segment%on_pe) cycle
        if (.not. associated(segment%tr_Reg)) cycle
        do i=is,ie
          if (segment%is_N_or_S .and. i>=segment%HI%isd .and. i<=segment%HI%ied) then
            J = segment%HI%JsdB
            do m=1,segment%tr_Reg%ntseg
              ntr_id = segment%tr_reg%Tr(m)%ntr_index
              if (segment%direction == OBC_DIRECTION_S) then
                Ts2y(i,ntr_id,j) = segment%tr_Reg%Tr(m)%tres(i,j,k)
              else
                Ts2y(i,ntr_id,j+1) = segment%tr_Reg%Tr(m)%tres(i,j,k)
              endif
            enddo
          endif
        enddo
      enddo
    endif ; endif

    ! X reconstruction
    do j=js,je ; if (domore_j_k(j,k)) then
      ! CFL at every zonal face this row needs, computed once per row.
      do I=is-1,ie
        if (uhh_in(I,j,k) >= 0.0) then ; i_up = I ; else ; i_up = I+1 ; endif
        if (h_k(i_up,j) > 0.0) then
          cfl_row(I) = abs(uhh_in(I,j,k)) / h_k(i_up,j)
        else
          cfl_row(I) = 0.0
        endif
      enddo
      do m=1,ntr ; do i=G%isd,G%ied ; Ts2x(i,m) = tk(i,j,k,m) ; enddo ; enddo
      if (associated(OBC)) then ; if (OBC%OBC_pe) then
        do n=1,OBC%number_of_segments
          segment=>OBC%segment(n)
          if (.not. segment%on_pe) cycle
          if (.not. associated(segment%tr_Reg)) cycle
          if (segment%is_E_or_W .and. j>=segment%HI%jsd .and. j<=segment%HI%jed) then
            I = segment%HI%IsdB
            do m=1,segment%tr_Reg%ntseg
              ntr_id = segment%tr_reg%Tr(m)%ntr_index
              if (segment%direction == OBC_DIRECTION_W) then
                Ts2x(i,ntr_id) = segment%tr_Reg%Tr(m)%tres(i,j,k)
              else
                Ts2x(i+1,ntr_id) = segment%tr_Reg%Tr(m)%tres(i,j,k)
              endif
            enddo
          endif
        enddo
      endif ; endif
      do m=1,ntr
        if ((advect_schemes(m) == ADVECT_WENO5) .or. (advect_schemes(m) == ADVECT_WENO7)) then
          do I=is-1,ie
            if (uhh_in(I,j,k) >= 0.0) then ; i_up = i ; else ; i_up = i+1 ; endif
            T7(:) = Ts2x(i_up-3:i_up+3,m)
            order3 = G%mask2dCu(I_up-2,j)*G%mask2dCu(I_up-1,j)*G%mask2dCu(I_up,j)*&
                      G%mask2dCu(I_up+1,j)
            order5 = order3*G%mask2dCu(I_up-3,j)*G%mask2dCu(I_up+2,j)
            order7 = 0.0
            if (advect_schemes(m) == ADVECT_WENO7) &
              order7 = order5*G%mask2dCu(I_up-4,j)*G%mask2dCu(I_up+3,j)
            if (order7 == 1.0) then
              if (uhh_in(I,j,k) >= 0.0) then
                call weno7_face(wq, T7, cfl_row(I))
              else
                call weno7_face(wq, T7(7:1:-1), cfl_row(I))
              endif
            elseif (order5 == 1.0) then
              if (uhh_in(I,j,k) >= 0.0) then
                call weno5_face(wq, T7, cfl_row(I))
              else
                call weno5_face(wq, T7(7:1:-1), cfl_row(I))
              endif
            else
              qext = G%mask2dCu(I_up,j)*G%mask2dCu(I_up-1,j)
              call PPM_reconstruction(wq, T7(3), T7(4), T7(5), uhh_in(I,j,k), cfl_row(I), qext)
            endif
            flux_x_out(I,j,m) = uhh_in(I,j,k)*wq
          enddo
        endif
      enddo
      if (associated(OBC)) then ; if (OBC%OBC_pe) then
        if (OBC%specified_u_BCs_exist_globally .or. OBC%open_u_BCs_exist_globally) then
          do n=1,OBC%number_of_segments
            segment=>OBC%segment(n)
            if (.not. segment%on_pe) cycle
            if (.not. associated(segment%tr_Reg)) cycle
            if (segment%is_E_or_W .and. j>=segment%HI%jsd .and. j<=segment%HI%jed) then
              I = segment%HI%IsdB
              ! Tracer fluxes are set to prescribed values only for inflows from masked areas.
              ! Now changing to simply fixed inflows.
              if ((uhr(I,j,k) > 0.0 .and. segment%direction == OBC_DIRECTION_W) .or. &
                (uhr(I,j,k) < 0.0 .and. segment%direction == OBC_DIRECTION_E)) then
                uhh_in(I,j,k) = uhr(I,j,k)
                do m=1,segment%tr_Reg%ntseg
                  ntr_id = segment%tr_reg%Tr(m)%ntr_index
                  flux_x_out(I,j,ntr_id) = uhh_in(I,j,k)*segment%tr_Reg%Tr(m)%tres(I,j,k)
                enddo
              endif
            endif
          enddo
        endif
        if (OBC%open_u_BCs_exist_globally) then
          do n=1,OBC%number_of_segments
            segment=>OBC%segment(n)
            if (.not. segment%on_pe) cycle
            I = segment%HI%IsdB
            if (segment%is_E_or_W .and. j>=segment%HI%jsd .and. j<=segment%HI%jed) then
              if (segment%specified) cycle
              if (.not. associated(segment%tr_Reg)) cycle

              ! Tracer fluxes are set to prescribed values only for inflows from masked areas.
              if ((uhr(I,j,k) > 0.0 .and. G%mask2dT(i,j) < 0.5) .or. &
                (uhr(I,j,k) < 0.0 .and. G%mask2dT(i+1,j) < 0.5)) then
                uhh_in(I,j,k) = uhr(I,j,k)
                do m=1,segment%tr_Reg%ntseg
                  ntr_id = segment%tr_reg%Tr(m)%ntr_index
                  flux_x_out(I,j,ntr_id) = uhh_in(I,j,k)*segment%tr_Reg%Tr(m)%tres(I,j,k)
                enddo
              endif
            endif
          enddo
        endif
      endif ; endif
    endif ; enddo ! j-loop

    ! Y reconstruction and flux assembly (face-centric), analogous to the X-direction above.
    do J=js-1,je
      if (.not. (domore_j_k(J,k) .or. domore_j_k(J+1,k))) cycle
      ! CFL at every meridional face this column needs, computed once per J (shared across tracers).
      do i=is,ie
        if (vhh_in(i,J,k) >= 0.0) then ; j_up = J ; else ; j_up = J+1 ; endif
        if (h_k(i,j_up) > 0.0) then
          cfl_v(i,J) = abs(vhh_in(i,J,k)) / h_k(i,j_up)
        else
          cfl_v(i,J) = 0.0
        endif
      enddo
      do m=1,ntr
        if ((advect_schemes(m) == ADVECT_WENO5) .or. (advect_schemes(m) == ADVECT_WENO7)) then
          do i=is,ie
            if (vhh_in(i,J,k) >= 0.0) then ; j_up = j ; else ; j_up = j+1 ; endif
            T7(:) = Ts2y(i,m,j_up-3:j_up+3)
            order3 = G%mask2dCv(i,J_up-2)*G%mask2dCv(i,J_up-1)*G%mask2dCv(i,J_up)*&
                      G%mask2dCv(i,J_up+1)
            order5 = order3*G%mask2dCv(i,J_up-3)*G%mask2dCv(i,J_up+2)
            order7 = 0.0
            if (advect_schemes(m) == ADVECT_WENO7) &
              order7 = order5*G%mask2dCv(i,J_up-4)*G%mask2dCv(i,J_up+3)
            if (order7 == 1.0) then
              if (vhh_in(i,J,k) >= 0.0) then
                call weno7_face(wq, T7, cfl_v(i,J))
              else
                call weno7_face(wq, T7(7:1:-1), cfl_v(i,J))
              endif
            elseif (order5 == 1.0) then
              if (vhh_in(i,J,k) >= 0.0) then
                call weno5_face(wq, T7, cfl_v(i,J))
              else
                call weno5_face(wq, T7(7:1:-1), cfl_v(i,J))
              endif
            else
              qext = G%mask2dCv(i,J_up)*G%mask2dCv(i,J_up-1)
              call PPM_reconstruction(wq, T7(3), T7(4), T7(5), vhh_in(i,J,k), cfl_v(i,J), qext)
            endif
            flux_y_out(i,J,m) = vhh_in(i,J,k)*wq
          enddo
        endif
      enddo
      if (associated(OBC)) then ; if (OBC%OBC_pe) then
        if (OBC%specified_v_BCs_exist_globally .or. OBC%open_v_BCs_exist_globally) then
          do n=1,OBC%number_of_segments
            segment=>OBC%segment(n)
            if (.not. segment%on_pe) cycle
            if (.not. segment%specified) cycle
            if (.not. associated(segment%tr_Reg)) cycle
            if (OBC%segment(n)%is_N_or_S) then
              if (J >= segment%HI%JsdB .and. J<= segment%HI%JedB) then
                do i=segment%HI%isd,segment%HI%ied
                  ! Tracer fluxes are set to prescribed values only for inflows from masked areas.
                  ! Now changing to simply fixed inflows.
                  if ((vhr(i,J,k) > 0.0) .and. (segment%direction == OBC_DIRECTION_S) .or. &
                      (vhr(i,J,k) < 0.0) .and. (segment%direction == OBC_DIRECTION_N)) then
                    vhh_in(i,J,k) = vhr(i,J,k)
                    do m=1,segment%tr_Reg%ntseg
                      ntr_id = segment%tr_reg%Tr(m)%ntr_index
                      flux_y_out(i,J,ntr_id) = vhh_in(i,J,k)*OBC%segment(n)%tr_Reg%Tr(m)%tres(i,J,k)
                    enddo
                  endif
                enddo
              endif
            endif
          enddo
        endif

        if (OBC%open_v_BCs_exist_globally) then
          do n=1,OBC%number_of_segments
            segment=>OBC%segment(n)
            if (.not. segment%on_pe) cycle
            if (segment%specified) cycle
            if (.not. associated(segment%tr_Reg)) cycle
            if (segment%is_N_or_S .and. J>=segment%HI%JsdB .and. J<=segment%HI%JedB) then
              do i=segment%HI%isd,segment%HI%ied
                if ((vhr(i,J,k) > 0.0 .and. G%mask2dT(i,j) < 0.5) .or. &
                  (vhr(i,J,k) < 0.0 .and. G%mask2dT(i,j+1) < 0.5)) then
                  vhh_in(i,J,k) = vhr(i,J,k)
                  do m=1,segment%tr_Reg%ntseg
                    ntr_id = segment%tr_reg%Tr(m)%ntr_index
                    flux_y_out(i,J,ntr_id) = vhh_in(i,J,k)*segment%tr_Reg%Tr(m)%tres(i,J,k)
                  enddo
                endif
              enddo
            endif
          enddo
        endif
      endif ; endif
    enddo ! J-loop

  endif ! do_recon

  ! Zhang & Shu (2010) positivity limiter: rescale each cell's outgoing fluxes, all together,
  ! so it cannot draw down more tracer mass in one step than it has available.
  ! pass_var(theta) fills halo theta from the neighbor PE so a shared PE-boundary face
  ! is scaled by the same factor on both sides.
  ! Every PE must call pass_var, including those with do_recon=.false.
  do m = 1, ntr
    if (.not. apply_lim(m)) cycle

    do j = G%jsd, G%jed ; do i = G%isd, G%ied ; theta(i,j) = 1.0 ; enddo ; enddo
    if (do_recon) then
      do j = js, je ; if (.not. domore_j_k(j,k)) cycle
        do i = is, ie
          if (h_k(i,j) <= 0.0) cycle
          hT = h_k(i,j) * max(tk(i,j,k,m), 0.0)
          outgoing = max(flux_x_out(I,j,m), 0.0) - min(flux_x_out(I-1,j,m), 0.0) &
                    + max(flux_y_out(i,J,m), 0.0) - min(flux_y_out(i,J-1,m), 0.0)
          if (outgoing > hT) theta(i,j) = hT / outgoing
        enddo
      enddo
    endif
    call pass_var(theta, G%Domain)
    if (do_recon) then
      do j = js, je ; if (.not. domore_j_k(j,k)) cycle
        do I = is-1, ie
          if (flux_x_out(I,j,m) > 0.0) then
            flux_x_out(I,j,m) = flux_x_out(I,j,m) * theta(i,j)
          elseif (flux_x_out(I,j,m) < 0.0) then
            flux_x_out(I,j,m) = flux_x_out(I,j,m) * theta(i+1,j)
          endif
        enddo
      enddo
      do J = js-1, je ; do i = is, ie
        if (.not. (domore_j_k(max(J,G%jsd),k) .or. domore_j_k(min(J+1,G%jed),k))) cycle
        if (flux_y_out(i,J,m) > 0.0) then
          flux_y_out(i,J,m) = flux_y_out(i,J,m) * theta(i,J)
        elseif (flux_y_out(i,J,m) < 0.0) then
          flux_y_out(i,J,m) = flux_y_out(i,J,m) * theta(i,J+1)
        endif
      enddo ; enddo
    endif
  enddo ! m

end subroutine compute_flux_2d

!> WENO5-Z + MP reconstruction at the right face of cell q(4), given a 7-point upwind-ordered stencil.
pure subroutine weno5_face(wf, q, cfl)
  real, intent(in)  :: q(7)      !< stencil ordered upwind to downwind [conc]
  real, intent(in)  :: cfl       !< absolute value of the advective CFL number at this face [nondim]
  real, intent(out) :: wf        !< reconstructed value at the right face of q(4) [conc]

  real :: P0, P1, P2                   ! sub-stencil polynomial reconstructions
  real :: b0, b1, b2                   ! smoothness indicators
  real :: w0, w1, w2                   ! nonlinear weights
  real :: tau, wnorm, bmin             ! WENO-Z indicators
  real :: dm2, dm1, dd0, dd1, dd2      ! second differences
  real :: dm4p, dm4m                   ! 4th-order minmod combinations
  real :: qul, qmp, qmd, qlc           ! MP limiter reference values
  real :: qmin, qmax, bM               ! monotone range
  real :: q0_min, q0_max               ! range of {q(4), qmp}, for the discontinuity flag
  real :: wpl, wmr                     ! PPM-fallback left/right edge states
  real :: dA, mA                       ! PPM-fallback edge difference/mean, for shape correction
  real :: a6                           ! PPM-fallback curvature
  real, parameter :: C1_6 = 1.0/6.0
  real, parameter :: d0 = 1.0/10.0, d1 = 6.0/10.0, d2 = 3.0/10.0
  real, parameter :: alpha = 2.0

  ! WENO5-Z sub-stencil reconstructions
  P0 = ((2.0*q(2) - 7.0*q(3)) + 11.0*q(4))*C1_6
  P1 = ((-q(3) + 5.0*q(4)) + 2.0*q(5))*C1_6
  P2 = ((2.0*q(4) + 5.0*q(5)) - q(6))*C1_6

  ! Smoothness indicators (Jiang & Shu 1996)
  b0 = (13.0/12.0)*(q(2) - 2.0*q(3) + q(4))**2 &
        + ( 1.0/ 4.0)*(q(2) - 4.0*q(3) + 3.0*q(4))**2
  b1 = (13.0/12.0)*(q(3) - 2.0*q(4) + q(5))**2 &
        + ( 1.0/ 4.0)*(q(3) - q(5))**2
  b2 = (13.0/12.0)*(q(4) - 2.0*q(5) + q(6))**2 &
        + ( 1.0/ 4.0)*(3.0*q(4) - 4.0*q(5) + q(6))**2

  bM = min(b0, b1, b2)

  ! WENO-Z nonlinear weights
  tau = abs(b0 - 2.0*b1 + b2)
  w0 = d0 * weight_fac(tau, b0)
  w1 = d1 * weight_fac(tau, b1)
  w2 = d2 * weight_fac(tau, b2)
  wnorm = 1.0 / (w0 + w1 + w2)
  wf = (w0*P0 + w1*P1 + w2*P2) * wnorm

  ! MP limiter (Suresh & Huynh 1997)
  qul  = q(4) + alpha * (q(4) - q(3))
  qmp  = q(4) + minmod2(q(5) - q(4), qul - q(4))

  dm2  = q(3) - 2.0*q(2) + q(1)
  dm1  = q(2) - 2.0*q(3) + q(4)
  dd0  = q(5) - 2.0*q(4) + q(3)
  dd1  = q(4) - 2.0*q(5) + q(6)
  dd2  = q(5) - 2.0*q(6) + q(7)

  dm4p = minmod6(4.0*dd0 - dd1, 4.0*dd1 - dd0, dd0, dd1, dm1, dd2)
  dm4m = minmod6(4.0*dm1 - dd0, 4.0*dd0 - dm1, dm1, dd0, dm2, dd1)

  qmd  = 0.5*(q(5) + q(4)) - 0.5*dm4p
  qlc  = 0.5*(3.0*q(4) - q(3)) + (4.0/3.0)*dm4m

  qmin = max(min(q(4), q(5), qmd), min(q(4), qul, qlc))
  qmax = min(max(q(4), q(5), qmd), max(q(4), qul, qlc))

  wf = min(max(qmin, qmax), wf) ; wf = max(min(qmin, qmax), wf)
  q0_min = min(q(4), qmp) ; q0_max = max(q(4), qmp)

  if (((qmax-qmin) > (q0_max-q0_min) .and. (tau > bM) ) ) then

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

    a6 = 6.*q(4) - 3. * (wpl + wmr) ! Curvature
    wf = (wpl - 0.5 * cfl * ((wpl - wmr) - a6 * (1. - 2./3. * cfl)))
  endif

end subroutine weno5_face

!> WENO7-Z + MP reconstruction at the right face of cell q(4), given a 7-point upwind-ordered stencil.
pure subroutine weno7_face(wf, q, cfl)
  real, intent(in)  :: q(7)     !< stencil ordered upwind to downwind [conc]
  real, intent(in)  :: cfl      !< absolute value of the advective CFL number at this face [nondim]
  real, intent(out) :: wf       !< reconstructed value at the right face of q(4) [conc]

  real :: P0, P1, P2, P3               ! sub-stencil polynomial reconstructions
  real :: b0, b1, b2, b3               ! smoothness indicators
  real :: bM                           ! smoothness indicator mean, for the discontinuity flag
  real :: w0, w1, w2, w3               ! nonlinear weights
  real :: tau, wnorm                   ! WENO-Z indicators
  real :: dm2, dm1, dd0, dd1, dd2      ! second differences
  real :: dm4p, dm4m                   ! 4th-order minmod combinations
  real :: qul, qmp, qmd, qlc           ! MP limiter reference values
  real :: qmin, qmax                   ! monotone range
  real :: q0_min, q0_max               ! range of {q(4), qmp}, for the discontinuity flag
  real :: wpl, wmr                     ! PPM-fallback left/right edge states
  real :: dA, mA                       ! PPM-fallback edge difference/mean, for shape correction
  real :: a6                           ! PPM-fallback curvature
  real, parameter :: C1_12 = 1.0/12.0
  real, parameter :: d0 = 1.0/35.0, d1 = 12.0/35.0, d2 = 18.0/35.0, d3 = 4.0/35.0
  real, parameter :: alpha = 2.0

  ! WENO7-Z sub-stencil reconstructions
  P0 = (-3.0*q(1) + 13.0*q(2) - 23.0*q(3) + 25.0*q(4)) * C1_12
  P1 = (q(2) - 5.0*q(3) + 13.0*q(4) + 3.0*q(5)) * C1_12
  P2 = (-q(3) + 7.0*q(4) + 7.0*q(5) - q(6)) * C1_12
  P3 = (3.0*q(4) + 13.0*q(5) - 5.0*q(6) + q(7)) * C1_12

  ! Smoothness indicators, in the compact expanded form of Balsara & Shu (2000, JCP),
  ! (also reproduced in Balsara, Garain & Shu 2016). Coefficients pre-divided by 1000
  ! (normalized out by weight_fac)
  b0 = ( q(1)*((0.547*q(1) - 3.882*q(2)) + (4.642*q(3) - 1.854*q(4))) &
       + q(2)*((7.043*q(2) - 17.246*q(3)) + 7.042*q(4))) &
       + ( q(3)*(11.003*q(3) - 9.402*q(4)) &
       + 2.107*q(4)**2 )
  b1 = ( q(2)*((0.267*q(2) - 1.642*q(3)) + (1.602*q(4) - 0.494*q(5))) &
       + q(3)*((2.843*q(3) - 5.966*q(4)) + 1.922*q(5))) &
       + ( q(4)*(3.443*q(4) - 2.522*q(5)) &
       + 0.547*q(5)**2 )
  b2 = ( q(3)*((0.547*q(3) - 2.522*q(4)) + (1.922*q(5) - 0.494*q(6))) &
       + q(4)*((3.443*q(4) - 5.966*q(5)) + 1.602*q(6))) &
       + ( q(5)*(2.843*q(5) - 1.642*q(6)) &
       + 0.267*q(6)**2 )
  b3 = ( q(4)*((2.107*q(4) - 9.402*q(5)) + (7.042*q(6) - 1.854*q(7))) &
       + q(5)*((11.003*q(5) - 17.246*q(6)) + 4.642*q(7))) &
       + ( q(6)*(7.043*q(6) - 3.882*q(7)) &
       + 0.547*q(7)**2 )

  bM = min(b0, b1, b2, b3)

  ! WENO7-Z nonlinear weights
  tau = abs((b0 - b3) + 3*(b1 - b2))
  w0 = d0 * weight_fac(tau, b0)
  w1 = d1 * weight_fac(tau, b1)
  w2 = d2 * weight_fac(tau, b2)
  w3 = d3 * weight_fac(tau, b3)
  wnorm = 1.0 / ((w0 + w1) + (w2 + w3))
  wf = ((w0*P0 + w1*P1) + (w2*P2 + w3*P3)) * wnorm

  ! MP limiter (Suresh & Huynh 1997)
  qul  = q(4) + alpha * (q(4) - q(3))
  qmp  = q(4) + minmod2(q(5) - q(4), qul - q(4))

  dm2  = q(3) - 2.0*q(2) + q(1)
  dm1  = q(2) - 2.0*q(3) + q(4)
  dd0  = q(5) - 2.0*q(4) + q(3)
  dd1  = q(4) - 2.0*q(5) + q(6)
  dd2  = q(5) - 2.0*q(6) + q(7)

  dm4p = minmod6(4.0*dd0 - dd1, 4.0*dd1 - dd0, dd0, dd1, dm1, dd2)
  dm4m = minmod6(4.0*dm1 - dd0, 4.0*dd0 - dm1, dm1, dd0, dm2, dd1)

  qmd  = 0.5*(q(5) + q(4)) - 0.5*dm4p
  qlc  = 0.5*(3.0*q(4) - q(3)) + (4.0/3.0)*dm4m

  qmin = max(min(q(4), q(5), qmd), min(q(4), qul, qlc))
  qmax = min(max(q(4), q(5), qmd), max(q(4), qul, qlc))

  wf = min(max(qmin, qmax), wf) ; wf = max(min(qmin, qmax), wf)
  q0_min = min(q(4), qmp) ; q0_max = max(q(4), qmp)

  if (((qmax-qmin) > (q0_max-q0_min) .and. (tau > bM) ) ) then

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

    a6 = 6.*q(4) - 3. * (wpl + wmr) ! Curvature
    wf = (wpl - 0.5 * cfl * ((wpl - wmr) - a6 * (1. - 2./3. * cfl)))
  endif

end subroutine weno7_face

pure subroutine ppmw5_reconstruction(wq, q, u, cfl, pos_def)
	real, intent(in) :: q(5)   !< tracer concentration from cell i-2 to i+2 [conc]
	real, intent(in) :: u      !< advective flux [H L2 ~> m3 or kg]
	real, intent(in) :: cfl(3) !< absolute value of the advective upwind-cell CFL number [nondim]
	logical, intent(in) :: pos_def !< If true, fall back to PPM:H3 when a WENO edge is negative
	real, intent(out) :: wq    !< weno flux  [conc]

	real :: P0, P1, P2         ! reconstructed polynomials
	real :: b0, b1, b2         ! smoothness indicator
	real :: w0, w1, w2         ! nonlinear weights
	real :: tau                ! Difference of smoothness indicators
	real, parameter :: C1_6 = 1.0/6.0             ! The ratio of 1/6 [nondim]
	real, parameter :: d0 = 1.0/10.0              ! The ratio of 1/10 [nondim]
	real, parameter :: d1 = 6.0/10.0              ! The ratio of 3/5 [nondim]
	real, parameter :: d2 = 3.0/10.0              ! The ratio of 3/10 [nondim]
	real, parameter :: cfl_disc_tol = 1.0e-2      ! Relative CFL-spread threshold used as an extra
	                                               ! discontinuity signal [nondim].
	real :: wnorm                                 ! Temporary variable
	real :: dm1, dd0, dd1, dm4p, dm4m             ! Temporary variables
	real :: qul, qmd, qlc, qmin, qmax, alpha      ! Temporary variables
	real :: qmp, q0_min, q0_max                   ! MP bounds
	real :: wpl, wmr, dA, mA, a6
	logical :: lim, disc_L, disc_R, disc_cfl, wide_L, wide_R, colella

	! Discontinuity signal from CFL variation across the local stencil (cfl(1:3)).
	! Relative to |cfl(2)| so it doesn't just fire on floating-point noise at tiny CFL.
	disc_cfl = (abs(maxval(cfl) - minval(cfl)) > cfl_disc_tol * max(abs(cfl(2)), 1.0e-6))

	! Left state at i+1/2
	P0 = ((2.0*q(1) - 7.0*q(2)) + 11.0*q(3))*C1_6
	b0 = (13.0/12.0)*(q(1) - 2.0*q(2) + q(3))**2 &
	      + ( 1.0/ 4.0)*(q(1) - 4.0*q(2) + 3.0*q(3))**2

	P1 = ((-q(2) + 5.0*q(3)) + 2.0*q(4))*C1_6
	b1 = (13.0/12.0)*(q(2) - 2.0*q(3) + q(4))**2 &
	      + ( 1.0/ 4.0)*(q(2) - q(4))**2

	P2 = ((2.0*q(3) + 5.0*q(4)) - q(5))*C1_6
	b2 = (13.0/12.0)*(q(3) - 2.0*q(4) + q(5))**2 &
	      + ( 1.0/ 4.0)*(3.0*q(3) - 4.0*q(4) + q(5))**2

	disc_L = disc_cfl .or. (abs(b2-b0) >= min(b0, b1, b2))

	! Nonlinear weights
	tau = abs(b0 - 2.0*b1 + b2)
	w0 = d0*weight_fac(tau, b0)
	w1 = d1*weight_fac(tau, b1)
	w2 = d2*weight_fac(tau, b2)

	wnorm = 1.0 / (w0 + w1 + w2)
	wpl = (w0*P0 + w1*P1 + w2*P2) * wnorm

	! MP limiter (Suresh & Huynh 1997, He et al. 2016)
	alpha = 2.0
	qul = q(3) + alpha*(q(3) - q(2))
	qmp = q(3) + minmod2((q(4)-q(3)), (qul-q(3)))

	dm1 = q(1) - 2.0*q(2) + q(3)
	dd0 = q(4) - 2.0*q(3) + q(2)
	dd1 = q(3) - 2.0*q(4) + q(5)

	dm4p = minmod4( (4.0*dd0 - dd1), (4.0*dd1 - dd0), dd0, dd1 )
	dm4m = minmod4( (4.0*dm1 - dd0), (4.0*dd0 - dm1), dm1, dd0 )
	qmd = 0.5*((q(4) + q(3)) - dm4p)
	qlc = 0.5*(3.0*q(3) - q(2)) + (4.0/3.0)*dm4m

	qmin = max(min(q(3), q(4), qmd), min(q(3), qul, qlc))
	qmax = min(max(q(3), q(4), qmd), max(q(3), qul, qlc))
	q0_min = min(q(3), qmp) ; q0_max = max(q(3), qmp)
	wpl = min(max(qmin, qmax), wpl) ; wpl = max(min(qmin, qmax), wpl)
	wide_L = ((qmax-qmin) > (q0_max-q0_min))

	! Right state at i-1/2
	P0 = ((2.0*q(5) - 7.0*q(4)) + 11.0*q(3))*C1_6
	b0 = (13.0/12.0)*(q(5) - 2.0*q(4) + q(3))**2 &
	      + ( 1.0/ 4.0)*(q(5) - 4.0*q(4) + 3.0*q(3))**2

	P1 = ((-q(4) + 5.0*q(3)) + 2.0*q(2))*C1_6
	b1 = (13.0/12.0)*(q(4) - 2.0*q(3) + q(2))**2 &
	      + ( 1.0/ 4.0)*(q(4) - q(2))**2

	P2 = ((2.0*q(3) + 5.0*q(2)) - q(1))*C1_6
	b2 = (13.0/12.0)*(q(3) - 2.0*q(2) + q(1))**2 &
	      + ( 1.0/ 4.0)*(3.0*q(3) - 4.0*q(2) + q(1))**2

	disc_R = disc_cfl .or. (abs(b2-b0) >= min(b0, b1, b2))

	! Nonlinear weights
	tau = abs(b0 - 2.0*b1 + b2)
	w0 = d0*weight_fac(tau, b0)
	w1 = d1*weight_fac(tau, b1)
	w2 = d2*weight_fac(tau, b2)

	wnorm = 1.0 / (w0 + w1 + w2)
	wmr = (w0*P0 + w1*P1 + w2*P2) * wnorm

	! MP limiter (Suresh & Huynh 1997, He et al. 2016)
	qul = q(3) + alpha*(q(3) - q(4))
	qmp = q(3) + minmod2((q(2)-q(3)), (qul-q(3)))

	dm1 = q(5) - 2.0*q(4) + q(3)
	dd0 = q(2) - 2.0*q(3) + q(4)
	dd1 = q(3) - 2.0*q(2) + q(1)

	dm4p = minmod4( (4.0*dd0 - dd1), (4.0*dd1 - dd0), dd0, dd1 )
	dm4m = minmod4( (4.0*dm1 - dd0), (4.0*dd0 - dm1), dm1, dd0 )

	qmd = 0.5*((q(2) + q(3)) - dm4p)
	qlc = 0.5*(3.0*q(3) - q(4)) + (4.0/3.0)*dm4m

	qmin = max(min(q(3), q(2), qmd), min(q(3), qul, qlc))
	qmax = min(max(q(3), q(2), qmd), max(q(3), qul, qlc))
	q0_min = min(q(3), qmp) ; q0_max = max(q(3), qmp)
	wmr = min(max(qmin, qmax), wmr) ; wmr = max(min(qmin, qmax), wmr)
	wide_R = ((qmax-qmin) > (q0_max-q0_min))

	! H3 fallback trigger
  lim = ((wide_L .and. disc_L) .or. (wide_R .and. disc_R)) .or. &
	      (pos_def .and. ((wpl < 0.0) .or. (wmr < 0.0)))

	if (lim) then
		wpl = ((-q(2) + 5.0*q(3)) + 2.0*q(4))*C1_6
		wpl = min(max(q(3), q(4)), wpl) ; wpl = max(min(q(3), q(4)), wpl)
		wmr = ((-q(4) + 5.0*q(3)) + 2.0*q(2))*C1_6
		wmr = min(max(q(3), q(2)), wmr) ; wmr = max(min(q(3), q(2)), wmr)
		dA = wpl - wmr ; mA = 0.5*( wpl + wmr )
		if ((q(4)-q(3))*(q(3)-q(2)) <= 0.0) then
			wmr = q(3) ; wpl = q(3)
		elseif ( dA*(q(3)-mA) > (dA*dA)/6. ) then
			wmr = (3.*q(3)) - 2.*wpl
		elseif ( dA*(q(3)-mA) < - (dA*dA)/6. ) then
			wpl = (3.*q(3)) - 2.*wmr
		endif
	endif

	a6 = 6.*q(3) - 3. * (wpl + wmr) ! Curvature
	if (u >= 0.0) then
		wq = (wpl - 0.5 * cfl(2) * ((wpl - wmr) - a6 * (1. - 2./3. * cfl(2))))
	else
		wq = (wmr + 0.5 * cfl(2) * ((wpl - wmr) + a6 * (1. - 2./3. * cfl(2))))
	endif

end subroutine ppmw5_reconstruction

!> PPM reconstruction at the upwind face of the donor cell.
pure subroutine PPM_reconstruction(wq_ppm, qm, q0, qp, u, cfl, qext)
  real, intent(in)  :: qm, q0, qp !< tracer concentration for 3-stencil wide [conc]
  real, intent(in)  :: u        !< advective flux [H L2 ~> m3 or kg]
  real, intent(in)  :: cfl      !< absolute value of the advective upwind-cell CFL number [nondim]
  real, intent(in)  :: qext     !< check local extrema
  real, intent(out) :: wq_ppm  !< PPM tracer concentration at the cell interface i+1/2  [conc]

  real :: aL, aR, dA, mA, a6

  wq_ppm = 0.0

  aL = (5.*q0 + (2.*qm - qp)) / 6.
  aL = max(min(q0, qm), aL) ; aL = min(max(q0, qm), aL)
  aR = (5.*q0 + (2.*qp - qm)) / 6.
  aR = max(min(q0, qp), aR) ; aR = min(max(q0, qp), aR)

  dA = aR - aL ; mA = 0.5*(aR + aL)
  if (qext*(qp - q0)*(q0 - qm) <= 0.) then
    aL = q0 ; aR = q0
  elseif (dA*(q0 - mA) > (dA*dA)/6.) then
    aL = (3.*q0) - 2.*aR
  elseif (dA*(q0 - mA) < -(dA*dA)/6.) then
    aR = (3.*q0) - 2.*aL
  endif

  a6 = 6.*q0 - 3.*(aR + aL)

  if (u >= 0.0) then
    wq_ppm = aR - 0.5*cfl*((aR - aL) - a6*(1. - 2./3.*cfl))
  else
    wq_ppm = aL + 0.5*cfl*((aR - aL) + a6*(1. - 2./3.*cfl))
  endif

end subroutine PPM_reconstruction

pure elemental function minmod2(a, b) result(r)
  real, intent(in) :: a, b
  real :: r
  r = 0.5 * (sign(1.0, a) + sign(1.0, b)) * min(abs(a), abs(b))
end function minmod2

pure elemental function minmod4(a, b, c, d) result(r)
  real, intent(in) :: a, b, c, d
  real :: r, s
  s = sign(1.0, a)
  if ((sign(1.0,b)==s) .and. (sign(1.0,c)==s) .and. (sign(1.0,d)==s)) then
    r = s * min(abs(a), abs(b), abs(c), abs(d))
  else
    r = 0.0
  endif
end function minmod4

pure elemental function minmod6(a, b, c, d, e, f) result(r)
  real, intent(in) :: a, b, c, d, e, f
  real :: r, s
  s = sign(1.0, a)
  if ((sign(1.0,b)==s) .and. (sign(1.0,c)==s) .and. &
      (sign(1.0,d)==s) .and. (sign(1.0,e)==s) .and. (sign(1.0,f)==s)) then
    r = s * min(abs(a), abs(b), abs(c), abs(d), abs(e), abs(f))
  else
    r = 0.0
  endif
end function minmod6

!> Compute the WENO-Z weight factor.
pure function weight_fac(tau, b) result(factor)
  real, intent(in) :: tau  !< Difference of the smoothness indicator [A ~> a]
  real, intent(in) :: b    !< The smoothness indicator [A ~> a]
  real :: factor

  factor = 1.0e20 ; if (abs(b) > 1.0e-20*tau) factor = (1 + tau / b)

end function weight_fac

!> \namespace mom_tracer_advect_weno
!!
!!  WENO5-Z and WENO7-Z tracer reconstruction schemes (Balsara et al. 2016, Borges et al. 2008)
!!  with MP monotonicity-preserving limiting (Suresh & Huynh 1997, He et al. 2016)
!!  and using RK3 time stepping.
!!  A near-discontinuity fallback uses a van Leer harmonic-mean reconstruction when the
!!  smoothness indicators detect strong variation (flag condition).
!!
!! WENO5 stencil (q(4) is the donor cell):
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
!! The advection scheme of some tracers can be set to be different
!! to that used by active tracers.  However, either all or none of
!! the tracers must use RK3 and WENO5-Z or WENO7-Z.

end module MOM_tracer_advect_weno
