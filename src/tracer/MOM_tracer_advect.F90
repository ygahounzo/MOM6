!>  This module contains the subroutines that advect tracers along coordinate surfaces.
module MOM_tracer_advect

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_cpu_clock,       only : cpu_clock_id, cpu_clock_begin, cpu_clock_end
use MOM_cpu_clock,       only : CLOCK_MODULE, CLOCK_ROUTINE
use MOM_diag_mediator,   only : post_data, query_averaging_enabled, diag_ctrl
use MOM_diag_mediator,   only : register_diag_field, safe_alloc_ptr, time_type
use MOM_domains,         only : sum_across_PEs, max_across_PEs
use MOM_domains,         only : create_group_pass, do_group_pass, group_pass_type, pass_var
use MOM_error_handler,   only : MOM_error, FATAL, WARNING, MOM_mesg, is_root_pe
use MOM_file_parser,     only : get_param, log_version, param_file_type
use MOM_grid,            only : ocean_grid_type
use MOM_open_boundary,   only : ocean_OBC_type, OBC_NONE, OBC_DIRECTION_E
use MOM_open_boundary,   only : OBC_DIRECTION_W, OBC_DIRECTION_N, OBC_DIRECTION_S
use MOM_open_boundary,   only : OBC_segment_type
use MOM_tracer_registry, only : tracer_registry_type, tracer_type
use MOM_unit_scaling,    only : unit_scale_type
use MOM_verticalGrid,    only : verticalGrid_type
implicit none ; private

#include <MOM_memory.h>

public advect_tracer
public tracer_advect_init
public tracer_advect_end

!> Control structure for this module
type, public :: tracer_advect_CS ; private
  real    :: dt                    !< The baroclinic dynamics time step [T ~> s].
  type(diag_ctrl), pointer :: diag !< A structure that is used to regulate the
                                   !< timing of diagnostic output.
  logical :: debug                 !< If true, write verbose checksums for debugging purposes.
  logical :: usePPM                !< If true, use PPM instead of PLM
  logical :: useHuynh              !< If true, use the Huynh scheme for PPM interface values
  logical :: useWENO5              !< if true, use WENO5 instead of PPM or PLM
  logical :: useWENO7              !< if true, use WENO7 instead of WENO5, PPM or PLM
  logical :: useHuynhStencilBug = .false. !< If true, use the incorrect stencil width.
                                   !! This is provided for compatibility with legacy simuations.
  type(group_pass_type) :: pass_uhr_vhr_t_hprev !< A structure used for group passes

  !>@{ Diagnostic IDs
  integer :: id_ppmTr_x  = -1
  integer :: id_ppmTr_y  = -1
  integer :: id_wenoTr_x  = -1
  integer :: id_wenoTr_y  = -1
  !>@}

end type tracer_advect_CS

!>@{ CPU time clocks
integer :: id_clock_advect
integer :: id_clock_pass
integer :: id_clock_sync
!>@}

contains

!> This routine time steps the tracer concentration using a
!! monotonic, conservative, weakly diffusive scheme.
subroutine advect_tracer(h_end, uhtr, vhtr, OBC, dt, G, GV, US, CS, Reg, x_first_in, &
                         vol_prev, max_iter_in, update_vol_prev, uhr_out, vhr_out)
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
  type(tracer_advect_CS),  pointer       :: CS    !< control structure for module
  type(tracer_registry_type), pointer    :: Reg   !< pointer to tracer registry
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

  real :: landvolfill                   ! An arbitrary? nonzero cell volume [H L2 ~> m3 or kg].
  logical :: use_PPM_stencil            ! If true, use the correct PPM stencil width.
  real :: Idt                           ! 1/dt [T-1 ~> s-1].
  logical :: domore_u(SZJ_(G),SZK_(GV))  ! domore_u and domore_v indicate whether there is more
  logical :: domore_v(SZJB_(G),SZK_(GV)) ! advection to be done in the corresponding row or column.
  logical :: domore_v2(SZI_(G),SZK_(GV)) ! advection to be done in the corresponding row or column.
  logical :: x_first            ! If true, advect in the x-direction first.
  integer :: max_iter           ! maximum number of iterations in each layer
  integer :: domore_k(SZK_(GV))
  integer :: stencil            ! stencil of the advection scheme
  integer :: nsten_halo         ! number of stencils that fit in the halos
  integer :: i, j, k, m, is, ie, js, je, isd, ied, jsd, jed, nz, itt, ntr, do_any
  integer :: isv, iev, jsv, jev ! The valid range of the indices.
  integer :: IsdB, IedB, JsdB, JedB

  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: ppm_y
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: weno_y
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: ppm_x
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: weno_x

  domore_u(:,:) = .false.
  domore_v(:,:) = .false.
  domore_v2(:,:) = .false.
  is  = G%isc ; ie  = G%iec ; js  = G%jsc ; je  = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB
  landvolfill = 1.0e-20         ! This is arbitrary, but must be positive.
  stencil = 2                   ! The scheme's stencil; 2 for PLM

  if (.not. associated(CS)) call MOM_error(FATAL, "MOM_tracer_advect: "// &
       "tracer_advect_init must be called before advect_tracer.")
  if (.not. associated(Reg)) call MOM_error(FATAL, "MOM_tracer_advect: "// &
       "register_tracer must be called before advect_tracer.")
  if (Reg%ntr==0) return
  call cpu_clock_begin(id_clock_advect)
  x_first = (MOD(G%first_direction,2) == 0)

  ! increase stencil size for Colella & Woodward PPM
  use_PPM_stencil = CS%usePPM .and. .not. CS%useHuynhStencilBug
  if (use_PPM_stencil) stencil = 3
  if (CS%useWENO5) stencil = 5
  if (CS%useWENO7) stencil = 7 

  ntr = Reg%ntr
  Idt = 1.0 / dt

  max_iter = 2*INT(CEILING(dt/CS%dt)) + 1

  if (present(max_iter_in)) max_iter = max_iter_in
  if (present(x_first_in))  x_first = x_first_in
  call cpu_clock_begin(id_clock_pass)
  call create_group_pass(CS%pass_uhr_vhr_t_hprev, uhr, vhr, G%Domain)
  call create_group_pass(CS%pass_uhr_vhr_t_hprev, hprev, G%Domain)
  do m=1,ntr
    call create_group_pass(CS%pass_uhr_vhr_t_hprev, Reg%Tr(m)%t, G%Domain)
  enddo
  call cpu_clock_end(id_clock_pass)

  !$OMP parallel default(shared)

  ! This initializes the halos of uhr and vhr because pass_vector might do
  ! calculations on them, even though they are never used.
  !$OMP do
  do k=1,nz
    do j=jsd,jed ; do I=IsdB,IedB ; uhr(I,j,k) = 0.0 ; enddo ; enddo
    do J=jsdB,jedB ; do i=Isd,Ied ; vhr(i,J,k) = 0.0 ; enddo ; enddo
    do j=jsd,jed ; do i=Isd,Ied ; hprev(i,j,k) = 0.0 ; enddo ; enddo
    domore_k(k)=1
    !  Put the remaining (total) thickness fluxes into uhr and vhr.
    do j=js,je ; do I=is-1,ie ; uhr(I,j,k) = uhtr(I,j,k) ; enddo ; enddo
    do J=js-1,je ; do i=is,ie ; vhr(i,J,k) = vhtr(i,J,k) ; enddo ; enddo
    if (.not. present(vol_prev)) then
    !   This loop reconstructs the thickness field the last time that the
    ! tracers were updated, probably just after the diabatic forcing.  A useful
    ! diagnostic could be to compare this reconstruction with that older value.
      do j=js,je ; do i=is,ie
        hprev(i,j,k) = max(0.0, G%areaT(i,j)*h_end(i,j,k) + &
             ((uhr(I,j,k) - uhr(I-1,j,k)) + (vhr(i,J,k) - vhr(i,J-1,k))))
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

  ! initialize diagnostic fluxes and tendencies
  !$OMP do
  do m=1,ntr
    if (associated(Reg%Tr(m)%ad_x)) Reg%Tr(m)%ad_x(:,:,:) = 0.0
    if (associated(Reg%Tr(m)%ad_y)) Reg%Tr(m)%ad_y(:,:,:) = 0.0
    if (associated(Reg%Tr(m)%advection_xy)) Reg%Tr(m)%advection_xy(:,:,:) = 0.0
    if (associated(Reg%Tr(m)%ad2d_x)) Reg%Tr(m)%ad2d_x(:,:) = 0.0
    if (associated(Reg%Tr(m)%ad2d_y)) Reg%Tr(m)%ad2d_y(:,:) = 0.0
  enddo
  !$OMP end parallel

  isv = is ; iev = ie ; jsv = js ; jev = je

  do itt=1,max_iter

    if (isv > is-stencil) then
      call do_group_pass(CS%pass_uhr_vhr_t_hprev, G%Domain, clock=id_clock_pass)

      nsten_halo = min(is-isd,ied-ie,js-jsd,jed-je)/stencil
      isv = is-nsten_halo*stencil ; jsv = js-nsten_halo*stencil
      iev = ie+nsten_halo*stencil ; jev = je+nsten_halo*stencil
      ! Reevaluate domore_u & domore_v unless the valid range is the same size as
      ! before.  Also, do this if there is Strang splitting.
      if ((nsten_halo > 1) .or. (itt==1)) then
        !$OMP parallel do default(shared)
        do k=1,nz ; if (domore_k(k) > 0) then
          do j=jsv,jev ; if (.not.domore_u(j,k)) then
            do i=isv+stencil-1,iev-stencil ; if (uhr(I,j,k) /= 0.0) then
              domore_u(j,k) = .true. ; exit
            endif ; enddo ! i-loop
          endif ; enddo

          !do i=isv,iev ; if (.not.domore_v2(i,k)) then
          !  do j=jsv+stencil-1,jev-stencil ; if (vhr(i,J,k) /= 0.0) then
          !    domore_v2(i,k) = .true. ; exit
          !  endif ; enddo ! j-loop
          !endif ; enddo

          do J=jsv+stencil-1,jev-stencil ; if (.not.domore_v(J,k)) then
            do i=isv+stencil,iev-stencil ; if (vhr(i,J,k) /= 0.0) then
              domore_v(J,k) = .true. ; exit
            endif ; enddo ! i-loop
          endif ; enddo

          !   At this point, domore_k is global.  Change it so that it indicates
          ! whether any work is needed on a layer on this processor.
          domore_k(k) = 0
          do j=jsv,jev ; if (domore_u(j,k)) domore_k(k) = 1 ; enddo
          do J=jsv+stencil-1,jev-stencil ; if (domore_v(J,k)) domore_k(k) = 1 ; enddo

        endif ; enddo ! k-loop
      endif
    endif
    
    ! Set the range of valid points after this iteration.
    isv = isv + stencil ; iev = iev - stencil
    jsv = jsv + stencil ; jev = jev - stencil

    !  To ensure positive definiteness of the thickness at each iteration, the
    !  mass fluxes out of each layer are checked each step, and limited to keep
    !  the thicknesses positive.  This means that several iterations may be required
    !  for all the transport to happen.  The sum over domore_k keeps the processors
    !  synchronized.  This may not be very efficient, but it should be reliable.

    !$OMP parallel default(shared)

    if (x_first) then

      !$OMP do ordered
      do k=1,nz ; if (domore_k(k) > 0) then
        ! First, advect zonally.
        call advect_x(Reg%Tr, hprev, uhr, uh_neglect, OBC, domore_u, ntr, Idt, &
                      isv, iev, jsv-stencil, jev+stencil, k, G, GV, US, CS%usePPM, &
                      CS%useHuynh, CS%useWENO5, CS%useWENO7, ppm_x, weno_x)
      endif ; enddo

      if(CS%id_ppmTr_x > 0) call post_data(CS%id_ppmTr_x, ppm_x, CS%diag) 

      !$OMP do ordered
      do k=1,nz ; if (domore_k(k) > 0) then
        !  Next, advect meridionally.
        call advect_y(Reg%Tr, hprev, vhr, vh_neglect, OBC, domore_v, ntr, Idt, &
                      isv, iev, jsv, jev, k, G, GV, US, CS%usePPM, CS%useHuynh, &
                      CS%useWENO5, CS%useWENO7, ppm_y, weno_y)

        !call advect_y2(Reg%Tr, hprev, vhr, vh_neglect, OBC, domore_v2, ntr, Idt, &
        !              isv, iev, jsv, jev, k, G, GV, US, CS%usePPM, CS%useHuynh, &
        !              CS%useWENO5, CS%useWENO7, ppm_y, weno_y)

        ! Update domore_k(k) for the next iteration
        domore_k(k) = 0
        do j=jsv-stencil,jev+stencil ; if (domore_u(j,k)) domore_k(k) = 1 ; enddo
        do J=jsv-1,jev ; if (domore_v(J,k)) domore_k(k) = 1 ; enddo
        !do i=isv,iev ; if (domore_v2(i,k)) domore_k(k) = 1 ; enddo

      endif ; enddo

      if(CS%id_ppmTr_y > 0) call post_data(CS%id_ppmTr_y, ppm_y, CS%diag) 

    else

      !$OMP do ordered
      do k=1,nz ; if (domore_k(k) > 0) then
        ! First, advect meridionally.
        call advect_y(Reg%Tr, hprev, vhr, vh_neglect, OBC, domore_v, ntr, Idt, &
                      isv-stencil, iev+stencil, jsv, jev, k, G, GV, US, CS%usePPM, &
                      CS%useHuynh, CS%useWENO5, CS%useWENO7, ppm_y, weno_y)

        !call advect_y2(Reg%Tr, hprev, vhr, vh_neglect, OBC, domore_v2, ntr, Idt, &
        !              isv-stencil, iev+stencil, jsv, jev, k, G, GV, US, CS%usePPM, &
        !              CS%useHuynh, CS%useWENO5, CS%useWENO7, ppm_y, weno_y)

      endif ; enddo

      if(CS%id_ppmTr_y > 0) call post_data(CS%id_ppmTr_y, ppm_y, CS%diag) 
      !if(CS%id_wenoTr_y > 0) call post_data(CS%id_wenoTr_y, weno_y, CS%diag) 

      !$OMP do ordered
      do k=1,nz ; if (domore_k(k) > 0) then
        ! Next, advect zonally.
        call advect_x(Reg%Tr, hprev, uhr, uh_neglect, OBC, domore_u, ntr, Idt, &
                      isv, iev, jsv, jev, k, G, GV, US, CS%usePPM, CS%useHuynh, &
                      CS%useWENO5, CS%useWENO7, ppm_x, weno_x)


        ! Update domore_k(k) for the next iteration
        domore_k(k) = 0
        do j=jsv,jev ; if (domore_u(j,k)) domore_k(k) = 1 ; enddo
        do J=jsv-1,jev ; if (domore_v(J,k)) domore_k(k) = 1 ; enddo
        !do i=isv-stencil,iev+stencil ; if (domore_v2(i,k)) domore_k(k) = 1 ; enddo
      endif ; enddo

      if(CS%id_ppmTr_x > 0) call post_data(CS%id_ppmTr_x, ppm_x, CS%diag) 

    endif ! x_first

    !$OMP end parallel

    ! If the advection just isn't finishing after max_iter, move on.
    if (itt >= max_iter) then
      exit
    endif

    ! Exit if there are no layers that need more iterations.
    if (isv > is-stencil) then
      do_any = 0
      call cpu_clock_begin(id_clock_sync)
      call sum_across_PEs(domore_k(:), nz)
      call cpu_clock_end(id_clock_sync)
      do k=1,nz ; do_any = do_any + domore_k(k) ; enddo
      if (do_any == 0) then
        exit
      endif

    endif

  enddo ! Iterations loop

  if (present(uhr_out)) uhr_out(:,:,:) = uhr(:,:,:)
  if (present(vhr_out)) vhr_out(:,:,:) = vhr(:,:,:)
  if (present(vol_prev) .and. present(update_vol_prev)) then
    if (update_vol_prev) vol_prev(:,:,:) = hprev(:,:,:)
  endif

  call cpu_clock_end(id_clock_advect)

end subroutine advect_tracer

!> This subroutine does 1-d flux-form advection in the zonal direction using
!! a monotonic piecewise linear scheme.
subroutine advect_x(Tr, hprev, uhr, uh_neglect, OBC, domore_u, ntr, Idt, &
                    is, ie, js, je, k, G, GV, US, usePPM, useHuynh, useWENO5, useWENO7, ppm_x, weno_x)
  type(ocean_grid_type),                     intent(inout) :: G    !< The ocean's grid structure
  type(verticalGrid_type),                   intent(in)    :: GV   !< The ocean's vertical grid structure
  integer,                                   intent(in)    :: ntr  !< The number of tracers
  type(tracer_type), dimension(ntr),         intent(inout) :: Tr   !< The array of registered tracers to work on
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: hprev !< cell volume at the end of previous
                                                                  !! tracer change [H L2 ~> m3 or kg]
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)), intent(inout) :: uhr !< accumulated volume/mass flux through
                                                                  !! the zonal face [H L2 ~> m3 or kg]
  real, dimension(SZIB_(G),SZJ_(G)),         intent(in)    :: uh_neglect !< A tiny zonal mass flux that can
                                                                  !! be neglected [H L2 ~> m3 or kg]
  type(ocean_OBC_type),                      pointer       :: OBC !< specifies whether, where, and what OBCs are used
  logical, dimension(SZJ_(G),SZK_(GV)),      intent(inout) :: domore_u !< If true, there is more advection to be
                                                                  !! done in this u-row
  real,                                      intent(in)    :: Idt !< The inverse of dt [T-1 ~> s-1]
  integer,                                   intent(in)    :: is  !< The starting tracer i-index to work on
  integer,                                   intent(in)    :: ie  !< The ending tracer i-index to work on
  integer,                                   intent(in)    :: js  !< The starting tracer j-index to work on
  integer,                                   intent(in)    :: je  !< The ending tracer j-index to work on
  integer,                                   intent(in)    :: k   !< The k-level to work on
  type(unit_scale_type),                     intent(in)    :: US  !< A dimensional unit scaling type
  logical,                                   intent(in)    :: usePPM !< If true, use PPM instead of PLM
  logical,                                   intent(in)    :: useHuynh !< If true, use the Huynh scheme
                                                                     !! for PPM interface values
  logical,                                   intent(in)    :: useWENO5 !< If true, use WENO5 instead of PPM or PLM
  logical,                                   intent(in)    :: useWENO7 !< If true, use WENO7 instead of WENO5, PPM or PLM

  real, intent(inout), dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: ppm_x
  real, intent(inout), dimension(SZIB_(G),SZJ_(G),SZK_(GV)) :: weno_x

  real, dimension(SZI_(G),ntr) :: &
    slope_x             ! The concentration slope per grid point [conc].
  real, dimension(SZIB_(G),SZJ_(G),ntr) :: &
    flux_x              ! The tracer flux across a boundary [H L2 conc ~> m3 conc or kg conc].
  real, dimension(SZI_(G),ntr) :: &
    T_tmp               ! The copy of the tracer concentration at constant i,k [conc].

  real :: hup, hlos     ! hup is the upwind volume, hlos is the
                        ! part of that volume that might be lost
                        ! due to advection out the other side of
                        ! the grid box, both in [H L2 ~> m3 or kg].
  real :: uhh(SZIB_(G)) ! The zonal flux that occurs during the
                        ! current iteration [H L2 ~> m3 or kg].
  real, dimension(SZIB_(G)) :: &
    hlst, &             ! Work variable [H L2 ~> m3 or kg].
    Ihnew, &            ! Work variable [H-1 L-2 ~> m-3 or kg-1].
    CFL                 ! The absolute value of the advective upwind-cell CFL number [nondim].
  real :: min_h         ! The minimum thickness that can be realized during
                        ! any of the passes [H ~> m or kg m-2].
  real :: tiny_h        ! The smallest numerically invertible thickness [H ~> m or kg m-2].
  real :: h_neglect     ! A thickness that is so small it is usually lost
                        ! in roundoff and can be neglected [H ~> m or kg m-2].
  real :: aR, aL        ! Reconstructed tracer concentrations at the right and left edges [conc]
  real :: dMx           ! Difference between the maximum of the surrounding cell concentrations and
                        ! the value in the cell whose reconstruction is being found [conc]
  real :: dMn           ! Difference between the tracer concentration in the cell whose reconstruction
                        ! is being found and the minimum of the surrounding values [conc]
  real :: Tp, Tc, Tm    ! Tracer concentrations around the upstream cell [conc]
  real :: dA            ! Difference between the reconstruction tracer edge values [conc]
  real :: mA            ! Average of the reconstruction tracer edge values [conc]
  real :: a6            ! Curvature of the reconstruction tracer values [conc]
  logical :: do_i(SZI_(G),SZJ_(G))     ! If true, work on given points.
  logical :: usePLMslope
  integer :: i, j, m, n, i_up, stencil, ntr_id
  type(OBC_segment_type), pointer :: segment=>NULL()
  logical, dimension(SZJ_(G),SZK_(GV)) :: domore_u_initial
  real :: wq, Tpp, Tmm, Tppp, Tmmm
  real :: do_plm, cst, bcv, D2w, D2wL, D2wR, D2wlim
  real :: order3, order5, order7
  real :: Tm3, Tm2, Tm1, Tp1, Tp2, Tp3, Tp4, ddx
  real :: wml,wpl,wpr, wmr, wppl, minT, maxT, mu, Tm4, Tm5, Tp5
  real :: Tmin(ntr), Tmax(ntr), u, pcm, appm, aweno, Ts, Te, dx, pcm1


  do_plm = 2.0

  ! keep a local copy of the initial values of domore_u, which is to be used when computing ad2d_x
  ! diagnostic at the end of this subroutine.
  domore_u_initial = domore_u

  usePLMslope = .not. (usePPM .and. useHuynh)
  ! stencil for calculating slope values
  stencil = 1
  if (usePPM .and. .not. useHuynh) stencil = 2

  min_h = 0.1*GV%Angstrom_H
  tiny_h = tiny(min_h)
  h_neglect = GV%H_subroundoff

  do I=is-1,ie ; CFL(I) = 0.0 ; enddo
  
  do m = 1,ntr
    Tmin(m) = 0.0 ; Tmax(m) = 1.0 
    !Tmin(m) = minval(Tr(m)%t(:,:,k)) ; Tmax(m) = maxval(Tr(m)%t(:,:,k))
  enddo

  do j=js,je ; if (domore_u(j,k)) then
    domore_u(j,k) = .false.

    ! Calculate the i-direction profiles (slopes) of each tracer that is being advected.
    if (usePLMslope) then
      do m=1,ntr ; do i=is-stencil,ie+stencil
       !if (ABS(Tr(m)%t(i+1,j,k)-Tr(m)%t(i,j,k)) < &
       !    ABS(Tr(m)%t(i,j,k)-Tr(m)%t(i-1,j,k))) then
       !  maxslope = 4.0*(Tr(m)%t(i+1,j,k)-Tr(m)%t(i,j,k))
       !else
       !  maxslope = 4.0*(Tr(m)%t(i,j,k)-Tr(m)%t(i-1,j,k))
       !endif
       !if ((Tr(m)%t(i+1,j,k)-Tr(m)%t(i,j,k)) * (Tr(m)%t(i,j,k)-Tr(m)%t(i-1,j,k)) < 0.0) then
       !  slope_x(i,m) = 0.0
       !elseif (ABS(Tr(m)%t(i+1,j,k)-Tr(m)%t(i-1,j,k))<ABS(maxslope)) then
       !  slope_x(i,m) = G%mask2dCu(I,j)*G%mask2dCu(I-1,j) * &
       !                 0.5*(Tr(m)%t(i+1,j,k)-Tr(m)%t(i-1,j,k))
       !else
       !  slope_x(i,m) = G%mask2dCu(I,j)*G%mask2dCu(I-1,j) * 0.5*maxslope
       !endif
        Tp = Tr(m)%t(i+1,j,k) ; Tc = Tr(m)%t(i,j,k) ; Tm = Tr(m)%t(i-1,j,k)
        dMx = max( Tp, Tc, Tm ) - Tc
        dMn= Tc - min( Tp, Tc, Tm )
        slope_x(i,m) = G%mask2dCu(I,j)*G%mask2dCu(I-1,j) * &
            sign( min(0.5*abs(Tp-Tm), 2.0*dMx, 2.0*dMn), Tp-Tm )
      enddo ; enddo
    endif ! usePLMslope

    ! make a copy of the tracers in case values need to be overridden for OBCs
    do m = 1,ntr
      do i=G%isd,G%ied
        T_tmp(i,m) = Tr(m)%t(i,j,k)
      enddo
    enddo
    ! loop through open boundaries and recalculate flux terms
    if (associated(OBC)) then ; if (OBC%OBC_pe) then
      do n=1,OBC%number_of_segments
        segment=>OBC%segment(n)
        if (.not. associated(segment%tr_Reg)) cycle
        if (segment%is_E_or_W) then
          if (j>=segment%HI%jsd .and. j<=segment%HI%jed) then
            I = segment%HI%IsdB
            do m = 1,segment%tr_Reg%ntseg ! replace tracers with OBC values
              ntr_id = segment%tr_reg%Tr(m)%ntr_index
              if (allocated(segment%tr_Reg%Tr(m)%tres)) then
                if (segment%direction == OBC_DIRECTION_W) then
                  T_tmp(i,ntr_id) = segment%tr_Reg%Tr(m)%tres(i,j,k)
                else
                  T_tmp(i+1,ntr_id) = segment%tr_Reg%Tr(m)%tres(i,j,k)
                endif
              else
                if (segment%direction == OBC_DIRECTION_W) then
                  T_tmp(i,ntr_id) = segment%tr_Reg%Tr(m)%OBC_inflow_conc
                else
                  T_tmp(i+1,ntr_id) = segment%tr_Reg%Tr(m)%OBC_inflow_conc
                endif
              endif
            enddo
            do m = 1,ntr ! Apply update tracer values for slope calculation
              do i=segment%HI%IsdB-1,segment%HI%IsdB+1
                Tp = T_tmp(i+1,m) ; Tc = T_tmp(i,m) ; Tm = T_tmp(i-1,m)
                dMx = max( Tp, Tc, Tm ) - Tc
                dMn= Tc - min( Tp, Tc, Tm )
                slope_x(i,m) = G%mask2dCu(I,j)*G%mask2dCu(I-1,j) * &
                     sign( min(0.5*abs(Tp-Tm), 2.0*dMx, 2.0*dMn), Tp-Tm )

              enddo
            enddo

          endif
        endif
      enddo
    endif ; endif


    ! Calculate the i-direction fluxes of each tracer, using as much
    ! the minimum of the remaining mass flux (uhr) and the half the mass
    ! in the cell plus whatever part of its half of the mass flux that
    ! the flux through the other side does not require.
    do I=is-1,ie
      if ((uhr(I,j,k) == 0.0) .or. &
          ((uhr(I,j,k) < 0.0) .and. (hprev(i+1,j,k) <= tiny_h)) .or. &
          ((uhr(I,j,k) > 0.0) .and. (hprev(i,j,k) <= tiny_h)) ) then
        uhh(I) = 0.0
        CFL(I) = 0.0
      elseif (uhr(I,j,k) < 0.0) then
        hup = hprev(i+1,j,k) - G%areaT(i+1,j)*min_h
        hlos = MAX(0.0, uhr(I+1,j,k))
        if ((((hup - hlos) + uhr(I,j,k)) < 0.0) .and. &
            ((0.5*hup + uhr(I,j,k)) < 0.0)) then
          uhh(I) = MIN(-0.5*hup, -hup+hlos, 0.0)
          domore_u(j,k) = .true.
        else
          uhh(I) = uhr(I,j,k)
        endif
        CFL(I) = - uhh(I) / (hprev(i+1,j,k))  ! CFL is positive
      else
        hup = hprev(i,j,k) - G%areaT(i,j)*min_h
        hlos = MAX(0.0, -uhr(I-1,j,k))
        if ((((hup - hlos) - uhr(I,j,k)) < 0.0) .and. &
            ((0.5*hup - uhr(I,j,k)) < 0.0)) then
          uhh(I) = MAX(0.5*hup, hup-hlos, 0.0)
          domore_u(j,k) = .true.
        else
          uhh(I) = uhr(I,j,k)
        endif
        CFL(I) = uhh(I) / (hprev(i,j,k))  ! CFL is positive
      endif

    enddo

    if (usePPM) then
      do m=1,ntr ; do I=is-1,ie
        ! centre cell depending on upstream direction
        if (uhh(I) >= 0.0) then
          i_up = i
        else
          i_up = i+1
        endif

        ! Implementation of PPM-H3
        Tp = T_tmp(i_up+1,m) ; Tc = T_tmp(i_up,m) ; Tm = T_tmp(i_up-1,m)

        if (useHuynh) then
          aL = ( 5.*Tc + ( 2.*Tm - Tp ) )/6. ! H3 estimate
          aL = max( min(Tc,Tm), aL) ; aL = min( max(Tc,Tm), aL) ! Bound
          aR = ( 5.*Tc + ( 2.*Tp - Tm ) )/6. ! H3 estimate
          aR = max( min(Tc,Tp), aR) ; aR = min( max(Tc,Tp), aR) ! Bound
        else
          aL = 0.5 * ((Tm + Tc) + (slope_x(i_up-1,m) - slope_x(i_up,m)) / 3.)
          aR = 0.5 * ((Tc + Tp) + (slope_x(i_up,m) - slope_x(i_up+1,m)) / 3.)
        endif
        
        dA = aR - aL ; mA = 0.5*( aR + aL )
        if (G%mask2dCu(I_up,j)*G%mask2dCu(I_up-1,j)*(Tp-Tc)*(Tc-Tm) <= 0.) then
          aL = Tc ; aR = Tc ! PCM for local extrema and boundary cells
        elseif ( dA*(Tc-mA) > (dA*dA)/6. ) then
          aL = (3.*Tc) - 2.*aR
        elseif ( dA*(Tc-mA) < - (dA*dA)/6. ) then
          aR = (3.*Tc) - 2.*aL
        endif

        a6 = 6.*Tc - 3. * (aR + aL) ! Curvature

        if (uhh(I) >= 0.0) then
          flux_x(I,j,m) = uhh(I)*( aR - 0.5 * CFL(I) * ( &
               ( aR - aL ) - a6 * ( 1. - 2./3. * CFL(I) ) ) )
        else
          flux_x(I,j,m) = uhh(I)*( aL + 0.5 * CFL(I) * ( &
               ( aR - aL ) + a6 * ( 1. - 2./3. * CFL(I) ) ) )
        endif
      enddo ; enddo

    elseif(useWENO5) then
      do m=1,ntr ; do I=is-1,ie
        
        i_up = i

        order3 = G%mask2dCu(I_up,j)*G%mask2dCu(I_up-1,j)*G%mask2dCu(I_up+1,j)*G%mask2dCu(I_up-2,j)*G%mask2dCu(I_up+2,j)
        order5 = order3*G%mask2dCu(I_up-3,j)*G%mask2dCu(I_up+3,j)
        
        Tm2 = T_tmp(i_up-2,m); Tm1 = T_tmp(i_up-1,m); Tc = T_tmp(i_up,m) ;
        Tp1 = T_tmp(i_up+1,m); Tp2 = T_tmp(i_up+2,m); Tp3 = T_tmp(i_up+3,m)

        u = uhh(I)
        mu = CFL(I)
        pcm = G%mask2dCu(I_up,j)*G%mask2dCu(I_up-1,j)*(Tp1-Tc)*(Tc-Tm1)
        pcm1 = G%mask2dCu(I_up+1,j)*G%mask2dCu(I_up,j)*(Tp2-Tp1)*(Tp1-Tc)
        appm = 0.0

        dx = G%dxT(i,j)

        if(order5 == 1.0) then
            call weno5_reconstruction_MPP(wq, Tm2, Tm1, Tc, Tp1, Tp2, Tp3, u, pcm, pcm1, mu, appm)
        elseif(order3 == 1.0) then
            call weno3_reconstruction_MPP(wq, Tm2, Tm1, Tc, Tp1, Tp2, Tp3, u, pcm, pcm1, mu, appm)
        else
            if(u >= 0.0) then
               wq = Tc
            else
               wq = Tp1
            endif 
        endif        

        flux_x(I,j,m) = uhh(I)*wq 
        ppm_x(I,j,k) = appm

      enddo ; enddo

    elseif(useWENO7) then
      do m=1,ntr ; do I=is-1,ie
        
        i_up = i
            
        order3 = G%mask2dCu(I_up,j)*G%mask2dCu(I_up-1,j)*G%mask2dCu(I_up+1,j)*G%mask2dCu(I_up-2,j)*G%mask2dCu(I_up+2,j)
        order5 = order3*G%mask2dCu(I_up-3,j)*G%mask2dCu(I_up+3,j)
        order7 = order5*G%mask2dCu(I_up-4,j)*G%mask2dCu(I_up+4,j)

        Tm3 = T_tmp(i_up-3,m); Tm2 = T_tmp(i_up-2,m); Tm1 = T_tmp(i_up-1,m); Tc = T_tmp(i_up,m) ;
        Tp1 = T_tmp(i_up+1,m); Tp2 = T_tmp(i_up+2,m); Tp3 = T_tmp(i_up+3,m); Tp4 = T_tmp(i_up+4,m)
        Tm5 = T_tmp(i_up-5,m); Tm4 = T_tmp(i_up-4,m); Tp5 = T_tmp(i_up+5,m)

        u = uhh(I)
        mu = CFL(I)
        pcm = G%mask2dCu(I_up,j)*G%mask2dCu(I_up-1,j)*(Tp1-Tc)*(Tc-Tm1)
        pcm1 = G%mask2dCu(I_up+1,j)*G%mask2dCu(I_up,j)*(Tp2-Tp1)*(Tp1-Tc)
        appm = 0.0

        if(order7 == 1.0) then
            call weno7_reconstruction_MPP(wq, Tm3, Tm2, Tm1, Tc, Tp1, Tp2, Tp3, Tp4, u, pcm, pcm1, mu, appm)
        elseif(order5 == 1.0) then
            call weno5_reconstruction_MPP(wq, Tm2, Tm1, Tc, Tp1, Tp2, Tp3, u, pcm, pcm1, mu, appm)
        elseif(order3 == 1.0) then
            call weno3_reconstruction_MPP(wq, Tm2, Tm1, Tc, Tp1, Tp2, Tp3, u, pcm, pcm1, mu, appm)
        else
            if(u > 0.0) then
               wq = Tc
            else
               wq = Tp1
            endif 
        endif        

        flux_x(I,j,m) = uhh(I)*wq
        ppm_x(I,j,k) = appm

      enddo ; enddo
    else ! PLM
      do m=1,ntr ; do I=is-1,ie
        if (uhh(I) >= 0.0) then
          ! Indirect implementation of PLM
         !aL = Tr(m)%t(i,j,k) - 0.5 * slope_x(i,m)
         !aR = Tr(m)%t(i,j,k) + 0.5 * slope_x(i,m)
         !flux_x(I,j,m) = uhh(I)*( aR - 0.5 * (aR-aL) * CFL(I) )
          ! Alternative implementation of PLM
          Tc = T_tmp(i,m)
          flux_x(I,j,m) = uhh(I)* ( Tc + 0.5 * slope_x(i,m) * ( 1. - CFL(I) ) )
        else
          ! Indirect implementation of PLM
         !aL = Tr(m)%t(i+1,j,k) - 0.5 * slope_x(i+1,m)
         !aR = Tr(m)%t(i+1,j,k) + 0.5 * slope_x(i+1,m)
         !flux_x(I,j,m) = uhh(I)*( aL + 0.5 * (aR-aL) * CFL(I) )
          ! Alternative implementation of PLM
          Tc = T_tmp(i+1,m)
          flux_x(I,j,m) = uhh(I)*( Tc - 0.5 * slope_x(i+1,m) * ( 1. - CFL(I) ) )
        endif
      enddo ; enddo
    endif ! usePPM
    
    if (associated(OBC)) then ; if (OBC%OBC_pe) then
      if (OBC%specified_u_BCs_exist_globally .or. OBC%open_u_BCs_exist_globally) then
        do n=1,OBC%number_of_segments
          segment=>OBC%segment(n)
          if (.not. associated(segment%tr_Reg)) cycle
          if (segment%is_E_or_W) then
            if (j>=segment%HI%jsd .and. j<=segment%HI%jed) then
              I = segment%HI%IsdB
              ! Tracer fluxes are set to prescribed values only for inflows from masked areas.
              ! Now changing to simply fixed inflows.
              if ((uhr(I,j,k) > 0.0) .and. (segment%direction == OBC_DIRECTION_W) .or. &
                  (uhr(I,j,k) < 0.0) .and. (segment%direction == OBC_DIRECTION_E)) then
                uhh(I) = uhr(I,j,k)
              ! should the reservoir evolve for this case Kate ?? - Nope
                do m=1,segment%tr_Reg%ntseg
                  ntr_id = segment%tr_reg%Tr(m)%ntr_index
                  if (allocated(segment%tr_Reg%Tr(m)%tres)) then
                    flux_x(I,j,ntr_id) = uhh(I)*segment%tr_Reg%Tr(m)%tres(I,j,k)
                  else ; flux_x(I,j,ntr_id) = uhh(I)*segment%tr_Reg%Tr(m)%OBC_inflow_conc ; endif
                enddo
              endif
            endif
          endif
        enddo
      endif

      if (OBC%open_u_BCs_exist_globally) then
        do n=1,OBC%number_of_segments
          segment=>OBC%segment(n)
          I = segment%HI%IsdB
          if (segment%is_E_or_W .and. (j >= segment%HI%jsd .and. j<= segment%HI%jed)) then
            if (segment%specified) cycle
            if (.not. associated(segment%tr_Reg)) cycle

            ! Tracer fluxes are set to prescribed values only for inflows from masked areas.
            if ((uhr(I,j,k) > 0.0) .and. (G%mask2dT(i,j) < 0.5) .or. &
                (uhr(I,j,k) < 0.0) .and. (G%mask2dT(i+1,j) < 0.5)) then
              uhh(I) = uhr(I,j,k)
              do m=1,segment%tr_Reg%ntseg
                ntr_id = segment%tr_reg%Tr(m)%ntr_index
                if (allocated(segment%tr_Reg%Tr(m)%tres)) then
                  flux_x(I,j,ntr_id) = uhh(I)*segment%tr_Reg%Tr(m)%tres(I,j,k)
                else; flux_x(I,j,ntr_id) = uhh(I)*segment%tr_Reg%Tr(m)%OBC_inflow_conc; endif
              enddo
            endif
          endif
        enddo
      endif
    endif ; endif

    ! Calculate new tracer concentration in each cell after accounting
    ! for the i-direction fluxes.
    do I=is-1,ie
      uhr(I,j,k) = uhr(I,j,k) - uhh(I)
      if (abs(uhr(I,j,k)) < uh_neglect(I,j)) uhr(I,j,k) = 0.0
    enddo
    do i=is,ie
      if ((uhh(I) /= 0.0) .or. (uhh(I-1) /= 0.0)) then
        do_i(i,j) = .true.
        hlst(i) = hprev(i,j,k)
        hprev(i,j,k) = hprev(i,j,k) - (uhh(I) - uhh(I-1))
        if (hprev(i,j,k) <= 0.0) then ; do_i(i,j) = .false.
        elseif (hprev(i,j,k) < h_neglect*G%areaT(i,j)) then
          hlst(i) = hlst(i) + (h_neglect*G%areaT(i,j) - hprev(i,j,k))
          Ihnew(i) = 1.0 / (h_neglect*G%areaT(i,j))
        else ;  Ihnew(i) = 1.0 / hprev(i,j,k) ; endif
      else
        do_i(i,j) = .false.
      endif
    enddo

    ! update tracer concentration from i-flux and save some diagnostics
    do m=1,ntr

      ! update tracer
      do i=is,ie
        if (do_i(i,j)) then
          if (Ihnew(i) > 0.0) then
            Tr(m)%t(i,j,k) = (Tr(m)%t(i,j,k) * hlst(i) - &
                              (flux_x(I,j,m) - flux_x(I-1,j,m))) * Ihnew(i)
          endif
        endif
      enddo

      ! diagnostics
      if (associated(Tr(m)%ad_x)) then ; do I=is-1,ie ; if (do_i(i,j) .or. do_i(i+1,j)) then
        Tr(m)%ad_x(I,j,k) = Tr(m)%ad_x(I,j,k) + flux_x(I,j,m)*Idt
      endif ; enddo ; endif

      ! diagnose convergence of flux_x (do not use the Ihnew(i) part of the logic).
      ! division by areaT to get into W/m2 for heat and kg/(s*m2) for salt.
      if (associated(Tr(m)%advection_xy)) then
        do i=is,ie ; if (do_i(i,j)) then
          Tr(m)%advection_xy(i,j,k) = Tr(m)%advection_xy(i,j,k) - (flux_x(I,j,m) - flux_x(I-1,j,m)) * &
                                          Idt * G%IareaT(i,j)
        endif ; enddo
      endif

    enddo

  endif ; enddo ! End of j-loop.

  ! Do user controlled underflow of the tracer concentrations.
  do m=1,ntr ; if (Tr(m)%conc_underflow > 0.0) then
    do j=js,je ; do i=is,ie
      if (abs(Tr(m)%t(i,j,k)) < Tr(m)%conc_underflow) Tr(m)%t(i,j,k) = 0.0
    enddo ; enddo
  endif ; enddo

  ! compute ad2d_x diagnostic outside above j-loop so as to make the summation ordered when OMP is active.

  !$OMP ordered
  do m=1,ntr ; if (associated(Tr(m)%ad2d_x)) then
    do j=js,je ; if (domore_u_initial(j,k)) then
      do I=is-1,ie ; if (do_i(i,j) .or. do_i(i+1,j)) then
        Tr(m)%ad2d_x(I,j) = Tr(m)%ad2d_x(I,j) + flux_x(I,j,m)*Idt
      endif ; enddo
    endif ; enddo
  endif ; enddo ! End of m-loop.
  !$OMP end ordered

end subroutine advect_x

!> This subroutine does 1-d flux-form advection using a monotonic piecewise
!! linear scheme.
subroutine advect_y(Tr, hprev, vhr, vh_neglect, OBC, domore_v, ntr, Idt, &
                    is, ie, js, je, k, G, GV, US, usePPM, useHuynh, useWENO5, useWENO7, ppm_y, weno_y)
  type(ocean_grid_type),                     intent(inout) :: G    !< The ocean's grid structure
  type(verticalGrid_type),                   intent(in)    :: GV   !< The ocean's vertical grid structure
  integer,                                   intent(in)    :: ntr !< The number of tracers
  type(tracer_type), dimension(ntr),         intent(inout) :: Tr   !< The array of registered tracers to work on
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: hprev !< cell volume at the end of previous
                                                                  !! tracer change [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)), intent(inout) :: vhr !< accumulated volume/mass flux through
                                                                  !! the meridional face [H L2 ~> m3 or kg]
  real, dimension(SZI_(G),SZJB_(G)),         intent(inout) :: vh_neglect !< A tiny meridional mass flux that can
                                                                  !! be neglected [H L2 ~> m3 or kg]
  type(ocean_OBC_type),                      pointer       :: OBC !< specifies whether, where, and what OBCs are used
  logical, dimension(SZJB_(G),SZK_(GV)),     intent(inout) :: domore_v !< If true, there is more advection to be
                                                                  !! done in this v-row
  real,                                      intent(in)    :: Idt !< The inverse of dt [T-1 ~> s-1]
  integer,                                   intent(in)    :: is  !< The starting tracer i-index to work on
  integer,                                   intent(in)    :: ie  !< The ending tracer i-index to work on
  integer,                                   intent(in)    :: js  !< The starting tracer j-index to work on
  integer,                                   intent(in)    :: je  !< The ending tracer j-index to work on
  integer,                                   intent(in)    :: k   !< The k-level to work on
  type(unit_scale_type),                     intent(in)    :: US  !< A dimensional unit scaling type
  logical,                                   intent(in)    :: usePPM !< If true, use PPM instead of PLM
  logical,                                   intent(in)    :: useHuynh !< If true, use the Huynh scheme
                                                                     !! for PPM interface values
  logical,                                   intent(in)    :: useWENO5 !< If true, use WENO5 instead of PPM or PLM
  logical,                                   intent(in)    :: useWENO7 !< If true, use WENO7 instead of WENO5, PPM or PLM

  real, intent(inout), dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: ppm_y
  real, intent(inout), dimension(SZI_(G),SZJB_(G),SZK_(GV)) :: weno_y

  real, dimension(SZI_(G),ntr,SZJ_(G)) :: &
    slope_y                     ! The concentration slope per grid point [conc].
  real, dimension(SZI_(G),ntr,SZJB_(G)) :: &
    flux_y                      ! The tracer flux across a boundary [H L2 conc ~> m3 conc or kg conc].
  real, dimension(SZI_(G),ntr,SZJB_(G)) :: &
    T_tmp               ! The copy of the tracer concentration at constant i,k [conc].
  real :: vhh(SZI_(G),SZJB_(G)) ! The meridional flux that occurs during the
                                ! current iteration [H L2 ~> m3 or kg].
  real :: hup, hlos             ! hup is the upwind volume, hlos is the
                                ! part of that volume that might be lost
                                ! due to advection out the other side of
                                ! the grid box, both in  [H L2 ~> m3 or kg].
  real, dimension(SZIB_(G)) :: &
    hlst, &             ! Work variable [H L2 ~> m3 or kg].
    Ihnew, &            ! Work variable [H-1 L-2 ~> m-3 or kg-1].
    CFL                 ! The absolute value of the advective upwind-cell CFL number [nondim].
  real :: min_h         ! The minimum thickness that can be realized during
                        ! any of the passes [H ~> m or kg m-2].
  real :: tiny_h        ! The smallest numerically invertible thickness [H ~> m or kg m-2].
  real :: h_neglect     ! A thickness that is so small it is usually lost
                        ! in roundoff and can be neglected [H ~> m or kg m-2].
  real :: aR, aL        ! Reconstructed tracer concentrations at the right and left edges [conc]
  real :: dMx           ! Difference between the maximum of the surrounding cell concentrations and
                        ! the value in the cell whose reconstruction is being found [conc]
  real :: dMn           ! Difference between the tracer average in the cell whose reconstruction
                        ! is being found and the minimum of the surrounding values [conc]
  real :: Tp, Tc, Tm    ! Tracer concentrations around the upstream cell [conc]
  real :: dA            ! Difference between the reconstruction tracer edge values [conc]
  real :: mA            ! Average of the reconstruction tracer edge values [conc]
  real :: a6            ! Curvature of the reconstruction tracer values [conc]
  logical :: do_j_tr(SZJ_(G))   ! If true, calculate the tracer profiles.
  logical :: do_i(SZI_(G), SZJ_(G))     ! If true, work on given points.
  logical :: usePLMslope
  integer :: i, j, j2, m, n, j_up, stencil, ntr_id
  type(OBC_segment_type), pointer :: segment=>NULL()
  logical :: domore_v_initial(SZJB_(G)) ! Initial state of domore_v
  real :: wq, Tpp, Tmm, Tppp, Tmmm, bcv, cst, D2w, D2wL, D2wR, D2wlim
  real :: order3, order5, order7
  real :: Tm3, Tm2, Tm1, Tp1, Tp2, Tp3, Tp4, ddx
  real :: wpl, wpr, wml, wmr, wppl, minT, maxT, mu, Tm4, Tm5, Tp5
  real :: Tmin(ntr), Tmax(ntr), v, pcm, appm, aweno, Ts, Te, dy, pcm1

  usePLMslope = .not. (usePPM .and. useHuynh)
  ! stencil for calculating slope values
  stencil = 1
  if (usePPM .and. .not. useHuynh) stencil = 2
  
  do m = 1,ntr
    Tmin(m) = 0.0 ; Tmax(m) = 1.0 
    !Tmin(m) = minval(Tr(m)%t(:,:,k)) ; Tmax(m) = maxval(Tr(m)%t(:,:,k))
  enddo

  min_h = 0.1*GV%Angstrom_H
  tiny_h = tiny(min_h)
  h_neglect = GV%H_subroundoff

  ! We conditionally perform work on tracer points: calculating the PLM slope,
  ! and updating tracer concentration within a cell
  ! this depends on whether there is a flux which would affect this tracer point,
  ! as indicated by domore_v. In the case of PPM reconstruction, a flux requires
  ! slope calculations at the two tracer points on either side (as indicated by
  ! the stencil variable), so we account for this with the do_j_tr flag array
  !
  ! Note: this does lead to unnecessary work in updating tracer concentrations,
  ! since that doesn't need a wider stencil with the PPM advection scheme, but
  ! this would require an additional loop, etc.
  do_j_tr(:) = .false.
  do J=js-1,je ; if (domore_v(J,k)) then ; do j2=1-stencil,stencil ; do_j_tr(j+j2) = .true. ; enddo ; endif ; enddo
  domore_v_initial(:) = domore_v(:,k)

  ! Calculate the j-direction profiles (slopes) of each tracer that
  ! is being advected.
  if (usePLMslope) then
    do j=js-stencil,je+stencil ; if (do_j_tr(j)) then ; do m=1,ntr ; do i=is,ie
      !if (ABS(Tr(m)%t(i,j+1,k)-Tr(m)%t(i,j,k)) < &
      !    ABS(Tr(m)%t(i,j,k)-Tr(m)%t(i,j-1,k))) then
      !  maxslope = 4.0*(Tr(m)%t(i,j+1,k)-Tr(m)%t(i,j,k))
      !else
      !  maxslope = 4.0*(Tr(m)%t(i,j,k)-Tr(m)%t(i,j-1,k))
      !endif
      !if ((Tr(m)%t(i,j+1,k)-Tr(m)%t(i,j,k))*(Tr(m)%t(i,j,k)-Tr(m)%t(i,j-1,k)) < 0.0) then
      !  slope_y(i,m,j) = 0.0
      !elseif (ABS(Tr(m)%t(i,j+1,k)-Tr(m)%t(i,j-1,k))<ABS(maxslope)) then
      !  slope_y(i,m,j) = G%mask2dCv(i,J) * G%mask2dCv(i,J-1) * &
      !                 0.5*(Tr(m)%t(i,j+1,k)-Tr(m)%t(i,j-1,k))
      !else
      !  slope_y(i,m,j) = G%mask2dCv(i,J) * G%mask2dCv(i,J-1) * 0.5*maxslope
      !endif
      Tp = Tr(m)%t(i,j+1,k) ; Tc = Tr(m)%t(i,j,k) ; Tm = Tr(m)%t(i,j-1,k)
      dMx = max( Tp, Tc, Tm ) - Tc
      dMn = Tc - min( Tp, Tc, Tm )
      slope_y(i,m,j) = G%mask2dCv(i,J)*G%mask2dCv(i,J-1) * &
           sign( min(0.5*abs(Tp-Tm), 2.0*dMx, 2.0*dMn), Tp-Tm )
    enddo ; enddo ; endif ; enddo ! End of i-, m-, & j- loops.
  endif ! usePLMslope


  ! make a copy of the tracers in case values need to be overridden for OBCs

  do j=G%jsd,G%jed ; do m=1,ntr ; do i=G%isd,G%ied
    T_tmp(i,m,j) = Tr(m)%t(i,j,k)
  enddo ; enddo ; enddo

  ! loop through open boundaries and recalculate flux terms
  if (associated(OBC)) then ; if (OBC%OBC_pe) then
    do n=1,OBC%number_of_segments
      segment=>OBC%segment(n)
      if (.not. associated(segment%tr_Reg)) cycle
      do i=is,ie
        if (segment%is_N_or_S) then
          if (i>=segment%HI%isd .and. i<=segment%HI%ied) then
            J = segment%HI%JsdB
            do m = 1,segment%tr_Reg%ntseg ! replace tracers with OBC values
              ntr_id = segment%tr_reg%Tr(m)%ntr_index
              if (allocated(segment%tr_Reg%Tr(m)%tres)) then
                if (segment%direction == OBC_DIRECTION_S) then
                  T_tmp(i,ntr_id,j) = segment%tr_Reg%Tr(m)%tres(i,j,k)
                else
                  T_tmp(i,ntr_id,j+1) = segment%tr_Reg%Tr(m)%tres(i,j,k)
                endif
              else
                if (segment%direction == OBC_DIRECTION_S) then
                  T_tmp(i,ntr_id,j) = segment%tr_Reg%Tr(m)%OBC_inflow_conc
                else
                  T_tmp(i,ntr_id,j+1) = segment%tr_Reg%Tr(m)%OBC_inflow_conc
                endif
              endif
            enddo
            do m = 1,ntr ! Apply update tracer values for slope calculation
              do j=segment%HI%JsdB-1,segment%HI%JsdB+1
                Tp = T_tmp(i,m,j+1) ; Tc = T_tmp(i,m,j) ; Tm = T_tmp(i,m,j-1)
                dMx = max( Tp, Tc, Tm ) - Tc
                dMn= Tc - min( Tp, Tc, Tm )
                slope_y(i,m,j) = G%mask2dCv(i,J)*G%mask2dCv(i,J-1) * &
                     sign( min(0.5*abs(Tp-Tm), 2.0*dMx, 2.0*dMn), Tp-Tm )
              enddo
            enddo
          endif
        endif ! is_N_S
      enddo ! i-loop
    enddo ! segment loop
  endif ; endif

  ! Calculate the j-direction fluxes of each tracer, using as much
  ! the minimum of the remaining mass flux (vhr) and the half the mass
  ! in the cell plus whatever part of its half of the mass flux that
  ! the flux through the other side does not require.
  do J=js-1,je ; if (domore_v(J,k)) then
    domore_v(J,k) = .false.

    do i=is,ie
      if ((vhr(i,J,k) == 0.0) .or. &
          ((vhr(i,J,k) < 0.0) .and. (hprev(i,j+1,k) <= tiny_h)) .or. &
          ((vhr(i,J,k) > 0.0) .and. (hprev(i,j,k) <= tiny_h)) ) then
        vhh(i,J) = 0.0
        CFL(i) = 0.0
      elseif (vhr(i,J,k) < 0.0) then
        hup = hprev(i,j+1,k) - G%areaT(i,j+1)*min_h
        hlos = MAX(0.0, vhr(i,J+1,k))
        if ((((hup - hlos) + vhr(i,J,k)) < 0.0) .and. &
            ((0.5*hup + vhr(i,J,k)) < 0.0)) then
          vhh(i,J) = MIN(-0.5*hup, -hup+hlos, 0.0)
          domore_v(J,k) = .true.
        else
          vhh(i,J) = vhr(i,J,k)
        endif
        CFL(i) = - vhh(i,J) / hprev(i,j+1,k)  ! CFL is positive
      else
        hup = hprev(i,j,k) - G%areaT(i,j)*min_h
        hlos = MAX(0.0, -vhr(i,J-1,k))
        if ((((hup - hlos) - vhr(i,J,k)) < 0.0) .and. &
            ((0.5*hup - vhr(i,J,k)) < 0.0)) then
          vhh(i,J) = MAX(0.5*hup, hup-hlos, 0.0)
          domore_v(J,k) = .true.
        else
          vhh(i,J) = vhr(i,J,k)
        endif
        CFL(i) = vhh(i,J) / hprev(i,j,k)  ! CFL is positive
      endif

    enddo

    if (usePPM) then
      do m=1,ntr ; do i=is,ie
        ! centre cell depending on upstream direction
        if (vhh(i,J) >= 0.0) then
          j_up = j
        else
          j_up = j + 1
        endif

        ! Implementation of PPM-H3
        Tp = T_tmp(i,m,j_up+1) ; Tc = T_tmp(i,m,j_up) ; Tm = T_tmp(i,m,j_up-1)

        if (useHuynh) then
          aL = ( 5.*Tc + ( 2.*Tm - Tp ) )/6. ! H3 estimate
          aL = max( min(Tc,Tm), aL) ; aL = min( max(Tc,Tm), aL) ! Bound
          aR = ( 5.*Tc + ( 2.*Tp - Tm ) )/6. ! H3 estimate
          aR = max( min(Tc,Tp), aR) ; aR = min( max(Tc,Tp), aR) ! Bound
        else
          aL = 0.5 * ((Tm + Tc) + (slope_y(i,m,j_up-1) - slope_y(i,m,j_up)) / 3.)
          aR = 0.5 * ((Tc + Tp) + (slope_y(i,m,j_up) - slope_y(i,m,j_up+1)) / 3.)
        endif
        
        dA = aR - aL ; mA = 0.5*( aR + aL )
        if (G%mask2dCv(i,J_up)*G%mask2dCv(i,J_up-1)*(Tp-Tc)*(Tc-Tm) <= 0.) then
          aL = Tc ; aR = Tc ! PCM for local extrema and boundary cells
        elseif ( dA*(Tc-mA) > (dA*dA)/6. ) then
          aL = (3.*Tc) - 2.*aR
        elseif ( dA*(Tc-mA) < - (dA*dA)/6. ) then
          aR = (3.*Tc) - 2.*aL
        endif

        a6 = 6.*Tc - 3. * (aR + aL) ! Curvature

        if (vhh(i,J) >= 0.0) then
          flux_y(i,m,J) = vhh(i,J)*( aR - 0.5 * CFL(i) * ( &
               ( aR - aL ) - a6 * ( 1. - 2./3. * CFL(i) ) ) )
        else
          flux_y(i,m,J) = vhh(i,J)*( aL + 0.5 * CFL(i) * ( &
               ( aR - aL ) + a6 * ( 1. - 2./3. * CFL(i) ) ) )
        endif
      enddo ; enddo

    elseif(useWENO5) then ! WENO5
      do m=1,ntr ; do i=is,ie

        j_up = j

        order3 = G%mask2dCv(i,J_up)*G%mask2dCv(i,J_up-1)*G%mask2dCv(i,J_up+1)*G%mask2dCv(i,J_up-2)*G%mask2dCv(i,J_up+2)
        order5 = order3*G%mask2dCv(i,J_up-3)*G%mask2dCv(i,J_up+3)

        Tm2 = T_tmp(i,m,j_up-2); Tm1 = T_tmp(i,m,j_up-1); Tc = T_tmp(i,m,j_up) ;
        Tp1 = T_tmp(i,m,j_up+1); Tp2 = T_tmp(i,m,j_up+2); Tp3 = T_tmp(i,m,j_up+3)

        v = vhh(i,J)
        mu = CFL(i)
        pcm = G%mask2dCv(i,J_up)*G%mask2dCv(i,J_up-1)*(Tp1-Tc)*(Tc-Tm1)
        pcm1 = G%mask2dCv(i,J_up+1)*G%mask2dCv(i,J_up)*(Tp2-Tp1)*(Tp1-Tc)
        appm = 0.0

        dy = G%dyT(i,j)

        if(order5 == 1.0) then
            call weno5_reconstruction_MPP(wq, Tm2, Tm1, Tc, Tp1, Tp2, Tp3, v, pcm, pcm1, mu, appm)
        elseif(order3 == 1.0) then
            call weno3_reconstruction_MPP(wq, Tm2, Tm1, Tc, Tp1, Tp2, Tp3, v, pcm, pcm1, mu, appm)
        else
            if(v >= 0.0) then
               wq = Tc
            else
               wq = Tp1
            endif 
        endif        

        flux_y(i,m,J) = vhh(i,J)*wq    

        ppm_y(i,J,k) = appm

      enddo ; enddo

    elseif(useWENO7) then ! WENO7
      do m=1,ntr ; do i=is,ie

        j_up = j

        order3 = G%mask2dCv(i,J_up)*G%mask2dCv(i,J_up-1)*G%mask2dCv(i,J_up+1)*G%mask2dCv(i,J_up-2)*G%mask2dCv(i,J_up+2)
        order5 = order3*G%mask2dCv(i,J_up-3)*G%mask2dCv(i,J_up+3)
        order7 = order5*G%mask2dCv(i,J_up-4)*G%mask2dCv(i,J_up+4)

        Tm3 = T_tmp(i,m,j_up-3); Tm2 = T_tmp(i,m,j_up-2); Tm1 = T_tmp(i,m,j_up-1); Tc = T_tmp(i,m,j_up) ;
        Tp1 = T_tmp(i,m,j_up+1); Tp2 = T_tmp(i,m,j_up+2); Tp3 = T_tmp(i,m,j_up+3); Tp4 = T_tmp(i,m,j_up+4)
        Tm5 = T_tmp(i,m,j_up-5); Tm4 = T_tmp(i,m,j_up-4); Tp5 = T_tmp(i,m,j_up+5)

        v = vhh(i,J)
        mu = CFL(i)
        pcm = G%mask2dCv(i,J_up)*G%mask2dCv(i,J_up-1)*(Tp1-Tc)*(Tc-Tm1)
        pcm1 = G%mask2dCv(i,J_up+1)*G%mask2dCv(i,J_up)*(Tp2-Tp1)*(Tp1-Tc)
        appm = 0.0

        if(order7 == 1.0) then
            call weno7_reconstruction_MPP(wq, Tm3, Tm2, Tm1, Tc, Tp1, Tp2, Tp3, Tp4, v, pcm, pcm1, mu, appm)
        elseif(order5 == 1.0) then
            call weno5_reconstruction_MPP(wq, Tm2, Tm1, Tc, Tp1, Tp2, Tp3, v, pcm, pcm1, mu, appm)
        elseif(order3 == 1.0) then
            call weno3_reconstruction_MPP(wq, Tm2, Tm1, Tc, Tp1, Tp2, Tp3, v, pcm, pcm1, mu, appm)
        else
            if(v >= 0.0) then
               wq = Tc
            else
               wq = Tp1
            endif 
        endif        

        flux_y(i,m,J) = vhh(i,J)*wq
        
        ppm_y(i,J,k) = appm

      enddo ; enddo
    else ! PLM
      do m=1,ntr ; do i=is,ie
        if (vhh(i,J) >= 0.0) then
          ! Indirect implementation of PLM
         !aL = Tr(m)%t(i,j,k) - 0.5 * slope_y(i,m,j)
         !aR = Tr(m)%t(i,j,k) + 0.5 * slope_y(i,m,j)
         !flux_y(i,m,J) = vhh(i,J)*( aR - 0.5 * (aR-aL) * CFL(i) )
          ! Alternative implementation of PLM
          Tc = T_tmp(i,m,j)
          flux_y(i,m,J) = vhh(i,J)*( Tc + 0.5 * slope_y(i,m,j) * ( 1. - CFL(i) ) )
        else
          ! Indirect implementation of PLM
         !aL = Tr(m)%t(i,j+1,k) - 0.5 * slope_y(i,m,j+1)
         !aR = Tr(m)%t(i,j+1,k) + 0.5 * slope_y(i,m,j+1)
         !flux_y(i,m,J) = vhh(i,J)*( aL + 0.5 * (aR-aL) * CFL(i) )
          ! Alternative implementation of PLM
          Tc = T_tmp(i,m,j+1)
          flux_y(i,m,J) = vhh(i,J)*( Tc - 0.5 * slope_y(i,m,j+1) * ( 1. - CFL(i) ) )
        endif
      enddo ; enddo
    endif ! usePPM

    if (associated(OBC)) then ; if (OBC%OBC_pe) then
      if (OBC%specified_v_BCs_exist_globally .or. OBC%open_v_BCs_exist_globally) then
        do n=1,OBC%number_of_segments
          segment=>OBC%segment(n)
          if (.not. segment%specified) cycle
          if (.not. associated(segment%tr_Reg)) cycle
          if (OBC%segment(n)%is_N_or_S) then
            if (J >= segment%HI%JsdB .and. J<= segment%HI%JedB) then
              do i=segment%HI%isd,segment%HI%ied
                ! Tracer fluxes are set to prescribed values only for inflows from masked areas.
                ! Now changing to simply fixed inflows.
                if ((vhr(i,J,k) > 0.0) .and. (segment%direction == OBC_DIRECTION_S) .or. &
                    (vhr(i,J,k) < 0.0) .and. (segment%direction == OBC_DIRECTION_N)) then
                  vhh(i,J) = vhr(i,J,k)
                  do m=1,segment%tr_Reg%ntseg
                    ntr_id = segment%tr_reg%Tr(m)%ntr_index
                    if (allocated(segment%tr_Reg%Tr(m)%tres)) then
                      flux_y(i,ntr_id,J) = vhh(i,J)*OBC%segment(n)%tr_Reg%Tr(m)%tres(i,J,k)
                    else ; flux_y(i,ntr_id,J) = vhh(i,J)*OBC%segment(n)%tr_Reg%Tr(m)%OBC_inflow_conc ; endif
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
          if (segment%specified) cycle
          if (.not. associated(segment%tr_Reg)) cycle
          if (segment%is_N_or_S .and. (J >= segment%HI%JsdB .and. J<= segment%HI%JedB)) then
            do i=segment%HI%isd,segment%HI%ied
              ! Tracer fluxes are set to prescribed values only for inflows from masked areas.
              if ((vhr(i,J,k) > 0.0) .and. (G%mask2dT(i,j) < 0.5) .or. &
                  (vhr(i,J,k) < 0.0) .and. (G%mask2dT(i,j+1) < 0.5)) then
                vhh(i,J) = vhr(i,J,k)
                do m=1,segment%tr_Reg%ntseg
                  ntr_id = segment%tr_reg%Tr(m)%ntr_index
                  if (allocated(segment%tr_Reg%Tr(m)%tres)) then
                    flux_y(i,ntr_id,J) = vhh(i,J)*segment%tr_Reg%Tr(m)%tres(i,J,k)
                  else ; flux_y(i,ntr_id,J) = vhh(i,J)*segment%tr_Reg%Tr(m)%OBC_inflow_conc ; endif
                enddo
              endif
            enddo
          endif
        enddo
      endif
    endif ; endif

  else ! not domore_v.
    do i=is,ie ; vhh(i,J) = 0.0 ; enddo
    do m=1,ntr ; do i=is,ie ; flux_y(i,m,J) = 0.0 ; enddo ; enddo
  endif ; enddo ! End of j-loop

  do J=js-1,je ; do i=is,ie
    vhr(i,J,k) = vhr(i,J,k) - vhh(i,J)
    if (abs(vhr(i,J,k)) < vh_neglect(i,J)) vhr(i,J,k) = 0.0
  enddo ; enddo

  ! Calculate new tracer concentration in each cell after accounting
  ! for the j-direction fluxes.
  do j=js,je ; if (do_j_tr(j)) then
    do i=is,ie
      if ((vhh(i,J) /= 0.0) .or. (vhh(i,J-1) /= 0.0)) then
        do_i(i,j) = .true.
        hlst(i) = hprev(i,j,k)
        hprev(i,j,k) = max(hprev(i,j,k) - (vhh(i,J) - vhh(i,J-1)), 0.0)
        if (hprev(i,j,k) <= 0.0) then ; do_i(i,j) = .false.
        elseif (hprev(i,j,k) < h_neglect*G%areaT(i,j)) then
          hlst(i) = hlst(i) + (h_neglect*G%areaT(i,j) - hprev(i,j,k))
          Ihnew(i) = 1.0 / (h_neglect*G%areaT(i,j))
        else ;  Ihnew(i) = 1.0 / hprev(i,j,k) ; endif
      else ; do_i(i,j) = .false. ; endif
    enddo

    ! update tracer and save some diagnostics
    do m=1,ntr
      do i=is,ie ; if (do_i(i,j)) then
        Tr(m)%t(i,j,k) = (Tr(m)%t(i,j,k) * hlst(i) - &
                          (flux_y(i,m,J) - flux_y(i,m,J-1))) * Ihnew(i)
      endif ; enddo

      ! diagnose convergence of flux_y and add to convergence of flux_x.
      ! division by areaT to get into W/m2 for heat and kg/(s*m2) for salt.
      if (associated(Tr(m)%advection_xy)) then
        do i=is,ie ; if (do_i(i,j)) then
          Tr(m)%advection_xy(i,j,k) = Tr(m)%advection_xy(i,j,k) - (flux_y(i,m,J) - flux_y(i,m,J-1))* Idt * &
                                          G%IareaT(i,j)
        endif ; enddo
      endif

    enddo
  endif ; enddo ! End of j-loop.

  ! Do user controlled underflow of the tracer concentrations.
  do m=1,ntr ; if (Tr(m)%conc_underflow > 0.0) then
    do j=js,je ; do i=is,ie
      if (abs(Tr(m)%t(i,j,k)) < Tr(m)%conc_underflow) Tr(m)%t(i,j,k) = 0.0
    enddo ; enddo
  endif ; enddo

  ! compute ad_y and ad2d_y diagnostic outside above j-loop so as to make the summation ordered when OMP is active.
  !$OMP ordered
  do m=1,ntr ; if (associated(Tr(m)%ad_y)) then
    do J=js-1,je ; if (domore_v_initial(J)) then
      do i=is,ie ; if (do_i(i,j) .or. do_i(i,j+1)) then
        Tr(m)%ad_y(i,J,k) = Tr(m)%ad_y(i,J,k) + flux_y(i,m,J)*Idt
      endif ; enddo
    endif ; enddo
  endif ; enddo ! End of m-loop.

  do m=1,ntr ; if (associated(Tr(m)%ad2d_y)) then
    do J=js-1,je ; if (domore_v_initial(J)) then
      do i=is,ie ; if (do_i(i,j) .or. do_i(i,j+1)) then
        Tr(m)%ad2d_y(i,J) = Tr(m)%ad2d_y(i,J) + flux_y(i,m,J)*Idt
      endif ; enddo
    endif ; enddo
  endif ; enddo ! End of m-loop.
  !$OMP end ordered

end subroutine advect_y

subroutine tracer_advect_init(Time, G, US, param_file, diag, CS)
  type(time_type), target, intent(in)    :: Time        !< current model time
  type(ocean_grid_type),   intent(in)    :: G           !< ocean grid structure
  type(unit_scale_type),   intent(in)    :: US          !< A dimensional unit scaling type
  type(param_file_type),   intent(in)    :: param_file  !< open file to parse for model parameters
  type(diag_ctrl), target, intent(inout) :: diag        !< regulates diagnostic output
  type(tracer_advect_CS),  pointer       :: CS          !< module control structure

  ! This include declares and sets the variable "version".
# include "version_variable.h"
  character(len=40)  :: mdl = "MOM_tracer_advect" ! This module's name.
  character(len=256) :: mesg    ! Message for error messages.

  if (associated(CS)) then
    call MOM_error(WARNING, "tracer_advect_init called with associated control structure.")
    return
  endif
  allocate(CS)

  CS%diag => diag

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "DT", CS%dt, fail_if_missing=.true., &
          desc="The (baroclinic) dynamics time step.", units="s", scale=US%s_to_T)
  call get_param(param_file, mdl, "DEBUG", CS%debug, default=.false.)
  call get_param(param_file, mdl, "TRACER_ADVECTION_SCHEME", mesg, &
          desc="The horizontal transport scheme for tracers:\n"//&
          "  PLM    - Piecewise Linear Method\n"//&
          "  PPM:H3 - Piecewise Parabolic Method (Huyhn 3rd order)\n"// &
          "  PPM    - Piecewise Parabolic Method (Colella-Woodward)\n"// &
          "  WENO5  - Weighted Essentially Non-Oscillatory, 5th order"&
          "  WENO7  - Weighted Essentially Non-Oscillatory, 7th order"&
          , default='PLM')
  select case (trim(mesg))
    case ("PLM")
      CS%usePPM = .false.
    case ("PPM:H3")
      CS%usePPM = .true.
      CS%useHuynh = .true.
    case ("PPM")
      CS%usePPM = .true.
      CS%useHuynh = .false.
    case ("WENO5")
      CS%useWENO5 = .true.
    case ("WENO7")
      CS%useWENO7 = .true.
    case default
      call MOM_error(FATAL, "MOM_tracer_advect, tracer_advect_init: "//&
           "Unknown TRACER_ADVECTION_SCHEME = "//trim(mesg))
  end select

  if (CS%useHuynh) then
    call get_param(param_file, mdl, "USE_HUYNH_STENCIL_BUG", &
        CS%useHuynhStencilBug, &
        desc="If true, use a stencil width of 2 in PPM:H3 tracer advection. " &
        // "This is incorrect and will produce regressions in certain " &
        // "configurations, but may be required to reproduce results in " &
        // "legacy simulations.", &
        default=.false.)
  endif

  id_clock_advect = cpu_clock_id('(Ocean advect tracer)', grain=CLOCK_MODULE)
  id_clock_pass = cpu_clock_id('(Ocean tracer halo updates)', grain=CLOCK_ROUTINE)
  id_clock_sync = cpu_clock_id('(Ocean tracer global synch)', grain=CLOCK_ROUTINE)

  CS%id_ppmTr_x  = -1
  CS%id_ppmTr_y  = -1
  CS%id_wenoTr_x  = -1
  CS%id_wenoTr_y  = -1

  CS%id_ppmTr_x = register_diag_field('ocean_model', 'tr_ppm_loc_x', diag%axesCuL, Time,&
          'Location where PPM is used along with WENO x-direction', 'nondim')
  CS%id_ppmTr_y = register_diag_field('ocean_model', 'tr_ppm_loc_y', diag%axesCvL, Time,&
          'Location where PPM is used along with WENO y-direction', 'nondim')
  CS%id_wenoTr_x = register_diag_field('ocean_model', 'tr_weno_loc_x', diag%axesCuL, Time,&
          'Location where WENO is used x-direction', 'nondim')
  CS%id_wenoTr_y = register_diag_field('ocean_model', 'tr_weno_loc_y', diag%axesCvL, Time,&
          'Location where WENO is used y-direction', 'nondim')
end subroutine tracer_advect_init

!> Close the tracer advection module
subroutine tracer_advect_end(CS)
  type(tracer_advect_CS), pointer :: CS  !< module control structure

  if (associated(CS)) deallocate(CS)

end subroutine tracer_advect_end

! weno reconstruction subroutines

subroutine weno3_reconstruction_MPP(wq,  qm2, qm, q0, qp, qp2, qp3, u, pcm, pcm1, mu, appm)

   real, intent(in) :: qm2, qm, q0, qp, qp2, qp3, u, pcm, mu, pcm1
   real, intent(out) :: wq, appm

   real :: wmr, wpl, w0

   w0 = 1.0

   if(u >= 0.0) then 
      call weno3_reconstruction(wmr, qp2, qp, q0, qm, qm2, pcm, mu, 1.0, appm)
      call weno3_reconstruction(wpl, qm2, qm, q0, qp, qp2, pcm, mu, 1.0, appm)
      call PP_limiter(wq, qm, q0, qp, wmr, wpl, w0)
   else
      call weno3_reconstruction(wpl, qp3, qp2, qp, q0, qm, pcm1, mu, -1.0, appm)
      call weno3_reconstruction(wmr, qm, q0, qp, qp2, qp3, pcm1, mu, -1.0, appm)
      call PP_limiter(wq, qp2, qp, q0, wmr, wpl, w0)
   endif

end subroutine weno3_reconstruction_MPP

subroutine weno3_reconstruction(wq, qmm, qm, q0, qp, qpp, pcm, mu, sig, appm)

   real, intent(in) :: qmm, qm, q0, qp, qpp, pcm, mu, sig
   real, intent(out) :: wq, appm

   real :: a1, a2, b1, b2, h1, h2, nu
   real :: eps, wnorm, w_1, w_2, P1, P2, tau
   
      call weno3_weights(b1, b2, qm, q0, qp)
      call weno3_poly(P1, P2, qm, q0, qp)

      wq = 0.0
      h1 = 1.0/3.0
      h2 = 2.0/3.0

      ! Alpha values
      eps = 1.0e-2
      tau = abs(b2-b1)
      a1 = h1*(1.0 + (tau/(b1+eps))**2)
      a2 = h2*(1.0 + (tau/(b2+eps))**2)

      ! Normalization
      wnorm = a1+a2
      w_1 = a1 / wnorm
      w_2 = a2 / wnorm

      wq = w_1*P1 + w_2*P2

      ! Monotonicity Preserving
      !call apply_MP(wq, qmm, qm, q0, qp, qpp, pcm, mu, sig, appm)

end subroutine weno3_reconstruction

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

subroutine weno5_reconstruction_MPP(wq, qm2, qm, q0, qp, qp2, qp3, u, pcm, pcm1, mu, appm)

   real, intent(in) :: qm2, qm, q0, qp, qp2, qp3, u, pcm, mu, pcm1
   real, intent(out) :: wq, appm

   real :: wmr, wpl, w0

   w0 = 5.0/18.0

   if(u >= 0.0) then
      call weno5_reconstruction(wmr, qp2, qp, q0, qm, qm2, pcm, mu, 1.0, appm)
      call weno5_reconstruction(wpl, qm2, qm, q0, qp, qp2, pcm, mu, 1.0, appm)
      call PP_limiter(wq, qm, q0, qp, wmr, wpl, w0)
   else
      call weno5_reconstruction(wpl, qp3, qp2, qp, q0, qm, pcm1, mu, -1.0, appm)
      call weno5_reconstruction(wmr, qm, q0, qp, qp2, qp3, pcm1, mu, -1.0, appm)
      call PP_limiter(wq, qp2, qp, q0, wmr, wpl, w0)
   endif

end subroutine weno5_reconstruction_MPP

subroutine weno5_reconstruction_MPP2(wq, qm2, qm, q0, qp, qp2, qp3, u, pcm, pcm1, mu, appm)

   real, intent(in) :: qm2, qm, q0, qp, qp2, qp3, u, pcm, mu, pcm1
   real, intent(out) :: wq, appm

   if(u >= 0.0) then
      call weno5_reconstruction(wq, qm2, qm, q0, qp, qp2, pcm, mu, 1.0, appm)
   else
      call weno5_reconstruction(wq, qp3, qp2, qp, q0, qm, pcm1, mu, -1.0, appm)
   endif

end subroutine weno5_reconstruction_MPP2

subroutine weno5_reconstruction(wq, qmm, qm, q0, qp, qpp, pcm, mu, sig, appm)

   real, intent(in) :: qmm, qm, q0, qp, qpp, pcm, mu, sig
   real, intent(out) :: wq, appm

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
   !call apply_MP(wq, qmm, qm, q0, qp, qpp, pcm, mu, sig, appm)

end subroutine weno5_reconstruction

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

subroutine weno7_reconstruction_MPP(wq, qm3, qm2, qm1, q0, qp1, qp2, qp3, qp4, u, pcm, pcm1, mu, appm)

   real, intent(in) :: qm3, qm2, qm1, q0, qp1, qp2, qp3, qp4, u, pcm, mu, pcm1
   real, intent(out) :: wq, appm
   real :: wmr, wpl, w0

   w0 = (322.0-13.0*sqrt(70.0))/1800.0

   if(u >= 0.0) then
      call weno7_reconstruction(wmr, qp3, qp2, qp1, q0, qm1, qm2, qm3, pcm, mu, 1.0, appm)
      call weno7_reconstruction(wpl, qm3, qm2, qm1, q0, qp1, qp2, qp3, pcm, mu, 1.0, appm)
      call PP_limiter(wq, qm1, q0, qp1, wmr, wpl, w0)
   else
      call weno7_reconstruction(wpl, qp4, qp3, qp2, qp1, q0, qm1, qm2, pcm1, mu, -1.0, appm)
      call weno7_reconstruction(wmr, qm2, qm1, q0, qp1, qp2, qp3, qp4, pcm1, mu, -1.0, appm)
      call PP_limiter(wq, qp2, qp1, q0, wmr, wpl, w0)
   endif

end subroutine weno7_reconstruction_MPP

subroutine weno7_reconstruction(wq, qm3, qm2, qm1, q0, qp1, qp2, qp3, pcm, mu, sig, appm)

   real, intent(in) :: qm3, qm2, qm1, q0, qp1, qp2, qp3, pcm, mu, sig
   real, intent(out) :: wq, appm

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
   !call apply_MP(wq, qm2, qm1, q0, qp1, qp2, pcm, mu, sig, appm)

end subroutine weno7_reconstruction

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

subroutine apply_MP(wq, qmm, qm, q0, qp, qpp, pcm, mu, sig, appm)

   real, intent(in) :: qmm, qm, q0, qp, qpp, pcm, mu, sig
   real, intent(inout) :: wq
   real, intent(out) :: appm

   real :: d0, d1, dm1, beta, ka
   real :: dm4,qlc,qmd,qul,qmin,qmax,md, dm
   real :: mmp, qmp, q0_min, q0_max, cp, r, slop
   real :: aL, aR, dA, mA, a6, s, wq_ppm, cp2, qt

   appm = 0.0
   wq_ppm = 0.0

   dm1 = qmm - 2.0*qm + q0
   d0 = qp - 2.0*q0 + qm
   d1 = qpp - 2.0*qp + q0

   call minmod3(dm4, 4.0*d0-d1, 4.0*d1-d0, d0, d1)
   call minmod3(dm, 4.0*dm1-d0, 4.0*d0-dm1, dm1, d0)
    
   beta = 4.0
   ka = 2.0
   !ka = (1.0-mu)/mu
   qul = q0 + ka*(q0-qm)
   qmd = 0.5*(q0 + qp) - 0.5*dm4
   qlc = q0 + 0.5*(q0-qm) + (beta/3.0)*dm
   !qlc = 0.5*(q0 + qul) + 0.5*ka*dm

   qmin = max(min(q0,qp,qmd),min(q0,qul,qlc))
   qmax = min(max(q0,qp,qmd),max(q0,qul,qlc))
   
   call minmod(mmp,qp-q0, ka*(q0-qm))
   qmp = q0 + mmp

   q0_min = min(q0,qmp) ; q0_max = max(q0,qmp)
   !cp = (qmax-qmin) - (q0_max-q0_min)

   cp = (wq-qmin)*(wq-qmax)
   cp2 = sign(1.0,cp)

   !if(cp2 >= 0.0) then
   !if(cp > 0.0) then
     
     !call apply_PPM(wq, qm, q0, qp, pcm, mu)

   !  aL = (5.0*q0 + (2.0*qm - qp))/6.0
   !  aL = max(min(q0,qm),aL) ; aL = min(max(q0,qm),aL)
   !  aR = (5.0*q0 + (2.0*qp - qm))/6.0
   !  aR = max(min(q0,qp),aR) ; aR = min(max(q0,qp),aR)

   !  dA = aR - aL ; mA = 0.5*(aR + aL)
   !  if(pcm <= 0.0) then
   !     aL = q0 ; aR = q0
   !  elseif(dA*(q0-mA) > (dA*dA)/6.0) then
   !     aL = (3.0*q0) - 2.0*aR
   !  elseif (dA*(q0-mA) < - (dA*dA)/6.0) then
   !     aR = (3.0*q0) - 2.0*aL
   !  endif

   !  a6 = 6.0*q0 - 3.0*(aR + aL) ! Curvature

   !  wq = aR - 0.5*mu*((aR - aL) - a6*(1.0 - (2.0/3.0)*mu))

   !  appm = 1.0

   !else 
   !  call median(md,wq,qmin,qmax)
   !  call median(md,qmin,wq,qmax)
   !  wq = md

   !endif

end subroutine apply_MP

subroutine PP_limiter(wq, qm, q0, qp, wmr, wpl, w0)

   real, intent(in) :: qm, q0, qp, wmr, wpl, w0
   real, intent(out) :: wq

   real :: qmin, qmax, theta, eps
   real :: qmin_g, qmax_g, eps2, P0, w1, wq_ppm

   P0 = (q0 - w0*wmr - w0*wpl)/(1.0 - 2.0*w0)

   qmin = min(wmr, P0, wpl)
   qmax = max(wmr, P0, wpl)

   qmin_g = 0.0; qmax_g = 1.0

   !eps = min(1.0e-13, (q0 - qmin)**6, q0)
   eps = min(1.0e-13, q0)
   theta = min(abs((qmax_g-q0)/(qmax-q0)), abs((qmin_g-q0+eps)/(qmin-q0)), 1.0)

   wq = theta*(wpl - q0) + q0

end subroutine PP_limiter

subroutine median(md,a,b,c)

   real, intent(in) :: a,b,c
   real, intent(out) :: md

   real :: mm, d1, d2

   d1 = b - a
   d2 = c - a
   
   call minmod(mm,d1,d2)
   md = a + mm

end subroutine median

subroutine minmod(mab,a,b)

   real, intent(in) :: a,b
   real, intent(out) :: mab
   
   real :: s1,s2

   mab = 0.5*(sign(1.0,a) + sign(1.0,b))*min(abs(a),abs(b))
   
end subroutine minmod 

subroutine minmod2(mab,a1,a2,a3,a4,a5,a6)

   real, intent(in) :: a1,a2,a3,a4,a5,a6
   real, intent(out) :: mab
   real :: s1,s2,s3,s4,s5,s6
   
   s1 = sign(1.0,a1); s2 = sign(1.0,a2); s3 = sign(1.0,a3)
   s4 = sign(1.0,a4); s5 = sign(1.0,a5); s6 = sign(1.0,a6)
   
   if(s1 == s2 .and. s2 == s3 .and. s3 == s4 .and. s4 == s5 .and. s5 == s6) then
     mab = s1*min(abs(a1),abs(a2),abs(a3),abs(a4),abs(a5),abs(a6))
   else
     mab = 0.0
   endif 

end subroutine minmod2 

subroutine minmod3(mab,a1,a2,a3,a4)

   real, intent(in) :: a1,a2,a3,a4
   real, intent(out) :: mab
   real :: s1,s2,s3,s4

   s1 = sign(1.0,a1); s2 = sign(1.0,a2); s3 = sign(1.0,a3)
   s4 = sign(1.0,a4)

   if(s1 == s2 .and. s2 == s3 .and. s3 == s4) then
     mab = s1*min(abs(a1),abs(a2),abs(a3),abs(a4))
   else
     mab = 0.0
   endif

end subroutine minmod3

!> \namespace mom_tracer_advect
!!
!!    This program contains the subroutines that advect tracers
!!  horizontally (i.e. along layers).
!!
!! \section section_mom_advect_intro
!!
!!  * advect_tracer advects tracer concentrations using a combination
!!  of the modified flux advection scheme from Easter (Mon. Wea. Rev.,
!!  1993) with tracer distributions given by the monotonic modified
!!  van Leer scheme proposed by Lin et al. (Mon. Wea. Rev., 1994).
!!  This scheme conserves the total amount of tracer while avoiding
!!  higher order accuracy scheme is needed, suggest monotonic
!!  piecewise parabolic method, as described in Carpenter et al.
!!  (MWR, 1990).
!!
!!  * advect_tracer has 4 arguments, described below. This
!!  subroutine determines the volume of a layer in a grid cell at the
!!  previous instance when the tracer concentration was changed, so
!!  it is essential that the volume fluxes should be correct.  It is
!!  also important that the tracer advection occurs before each
!!  calculation of the diabatic forcing.

end module MOM_tracer_advect
