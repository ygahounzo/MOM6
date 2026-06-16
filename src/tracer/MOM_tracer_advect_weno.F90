!>  This module contains the subroutines of the WENO schemes that advect tracers along coordinate surfaces.
module MOM_tracer_advect_weno

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
   use MOM_tracer_advect_schemes, only : ADVECT_WENO5, ADVECT_WENO7

   implicit none ; private

#include <MOM_memory.h>

   public weno5_reconstruction, ppmw5_reconstruction
   public weno7_reconstruction, ppmw7_reconstruction
   public PPM_reconstruction
   public rk3_substep

contains

!> One SSP-RK3 2D-unsplit sub-step for use inside advect_tracer_v0.
!! uhr/vhr must already be pre-scaled to 1/nsub of the total flux by the caller.
!! Unlike rk3_unsplit_step, this routine does NOT subtract uhh/vhh from uhr/vhr.
!! domore_k is updated from a direct flux check on uhr/vhr at the end of each layer;
!! because uhr = uhtr/nsub is constant across sub-steps, the caller does not need
!! to reset domore_k between sub-steps.

subroutine rk3_substep(G, GV, US, OBC, Reg, hprev, uhr, vhr, &
  uh_neglect, vh_neglect, domore_k, domore_j, &
  ntr, nz, isv, iev, jsv, jev, dump_cfl, &
    local_advect_scheme, Idt, CFL_subcycle)
  type(ocean_grid_type),      intent(in)    :: G
  type(verticalGrid_type),    intent(in)    :: GV
  type(unit_scale_type),      intent(in)    :: US
  type(ocean_OBC_type),       pointer       :: OBC
  type(tracer_registry_type), pointer       :: Reg
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
  integer, dimension(SZK_(GV)),                intent(inout) :: domore_k
  logical, dimension(SZJ_(G),SZK_(GV)),        intent(in)    :: domore_j !< per-row active flag
  integer,                                     intent(in)    :: ntr      !< The number of tracers
  integer,                                     intent(in)    :: nz
  integer,                                     intent(in)    :: isv   !< The starting tracer i-index to work on
  integer,                                     intent(in)    :: iev   !< The ending tracer i-index to work on
  integer,                                     intent(in)    :: jsv   !< The starting tracer j-index to work on
  integer,                                     intent(in)    :: jev   !< The ending tracer j-index to work on
  logical,                                     intent(inout)    :: dump_cfl !< flag for dumping the cfl
  integer, dimension(ntr),                     intent(in)    :: local_advect_scheme  !< list of advection schemes to use
  real,                                        intent(in)    :: Idt  !< The inverse of dt [T-1 ~> s-1]
  real,                                        intent(in)    :: CFL_subcycle  !< per-subcycle CFL limit [nondim]

  real :: uhh(SZIB_(G),SZJ_(G),SZK_(GV))
  real :: vhh(SZI_(G),SZJB_(G),SZK_(GV))
  real :: flux_x(SZIB_(G),SZJ_(G),ntr,nz)
  real :: flux_y(SZI_(G),SZJB_(G),ntr,nz)
  real :: flux_xs(SZIB_(G),SZJ_(G),ntr,nz)
  real :: flux_ys(SZI_(G),SZJB_(G),ntr,nz)
  real :: Ts1_s(SZI_(G),SZJ_(G),SZK_(GV),ntr)
  real :: Ts2_s(SZI_(G),SZJ_(G),SZK_(GV),ntr)
  real :: h_old, h_new, Ihnew_ij, dh_ij, h_neglect
  logical :: do_ij
  integer :: i, j, k, m, n
  real :: hprev_s1(SZI_(G),SZJ_(G),nz), hprev_s2(SZI_(G),SZJ_(G),nz)
  type(group_pass_type) :: pass_Ts1_hprev_s1, pass_Ts2_hprev_s2
  real :: CFL_cell_ij                     !< per-face upwind-cell outflow CFL [nondim]
  real :: scale(SZI_(G),SZJ_(G),nz)      !< outflow scale factor: min(1, CFL_subcycle/CFL_cell) [nondim]
  real :: dh(SZI_(G),SZJ_(G),nz)            !< precomputed flux divergence per cell [H]
  logical :: no_flux(SZI_(G),SZJ_(G),nz)    !< .true. if all four face fluxes are zero
  real :: tiny_h
  logical :: apply_lim_zs(ntr)  !< per-tracer Zhang-Shu limiter flag
  type(OBC_segment_type), pointer :: segment => null()
  integer :: m_zs

  h_neglect = GV%H_subroundoff
  tiny_h = GV%Angstrom_H

  do m_zs=1,ntr ; apply_lim_zs(m_zs) = Reg%Tr(m_zs)%nonneg_lim ; enddo

  flux_x = 0.0 ; flux_y = 0.0

  ! Stage 1: initialize T^n snapshot, compute mass fluxes,
  ! reconstruct stage-1 fluxes, and compute T* and h*.
  !$OMP parallel do default(shared) private(i,j,m,CFL_cell_ij,h_old,h_new,Ihnew_ij,do_ij)
  do k=1,nz
    if (domore_k(k) > 0) then

      do j=G%jsd,G%jed ; do i=G%isd,G%ied
        hprev_s1(i,j,k) = hprev(i,j,k)
        hprev_s2(i,j,k) = hprev(i,j,k)
        scale(i,j,k) = 1.0
      enddo ; enddo
      ! Default T* and T** to T^n for inactive cells.
      do m=1,ntr
        do j=G%jsd,G%jed ; do i=G%isd,G%ied
          Ts1_s(i,j,k,m) = Reg%Tr(m)%t(i,j,k)
          Ts2_s(i,j,k,m) = Reg%Tr(m)%t(i,j,k)
        enddo ; enddo
      enddo

      ! Compute per-cell outflow CFL from uhr/vhr and derive scale in one pass.
      ! Thin cells and inactive rows keep CFL=0 and scale=1.
      do j = jsv,jev ; if (domore_j(j,k)) then
        do i=isv,iev
          if ((hprev(i,j,k)*Idt > G%areaT(i,j)*tiny_h)) then
            CFL_cell_ij = (max(G%mask2dCu(I,j)*uhr(I,j,k), 0.0) &
                          -min(G%mask2dCu(I-1,j)*uhr(I-1,j,k), 0.0) &
                          +max(G%mask2dCv(i,J)*vhr(i,J,k), 0.0) &
                          -min(G%mask2dCv(i,J-1)*vhr(i,J-1,k), 0.0)) / hprev(i,j,k)
            if (CFL_cell_ij > CFL_subcycle) scale(i,j,k) = CFL_subcycle / CFL_cell_ij
          endif
        enddo
      endif ; enddo

      ! Compute mass fluxes with thin-cell zeroing and CFL scale in a single pass.
      ! uhh/vhh are written exactly once.
      do j=jsv,jev ; do I=isv-1,iev
        if (.not. domore_j(j,k)) then ; uhh(I,j,k) = 0.0 ; cycle ; endif
        if ((uhr(I,j,k) == 0.0) .or. &
            ((uhr(I,j,k) < 0.0) .and. (hprev(i+1,j,k)*Idt <= G%areaT(i+1,j)*tiny_h)) .or. &
            ((uhr(I,j,k) > 0.0) .and. (hprev(i,j,k)*Idt <= G%areaT(i,j)*tiny_h))) then
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
            ((G%mask2dCv(i,J)*vhr(i,J,k) < 0.0) .and. (hprev(i,j+1,k)*Idt <= G%areaT(i+1,j)*tiny_h)) .or. &
            ((G%mask2dCv(i,J)*vhr(i,J,k) > 0.0) .and. (hprev(i,j,k)*Idt <= G%areaT(i,j)*tiny_h))) then
          vhh(i,J,k) = 0.0
        elseif (G%mask2dCv(i,J)*vhr(i,J,k) > 0.0) then
          vhh(i,J,k) = G%mask2dCv(i,J)*vhr(i,J,k) * scale(i,J,k)
        else
          vhh(i,J,k) = G%mask2dCv(i,J)*vhr(i,J,k) * scale(i,J+1,k)
        endif
      enddo ; enddo

      call compute_flux_2d(Ts1_s, uhh, vhh, hprev(:,:,k), OBC, ntr, &
          isv, iev, jsv, jev, k, G, GV, local_advect_scheme, &
          flux_x(:,:,:,k), flux_y(:,:,:,k), domore_j, uhr, vhr, apply_lim_zs)

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
          hprev_s1(i,j,k) = h_new
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
            Ts1_s(i,j,k,m) = (Reg%Tr(m)%t(i,j,k)*h_old &
                            - (flux_x(I,j,m,k) - flux_x(I-1,j,m,k)) &
                            - (flux_y(i,J,m,k) - flux_y(i,J-1,m,k))) * Ihnew_ij
        enddo
      enddo ; enddo
    endif ! domore_k
  enddo ! Stage 1 k-loop
  dump_cfl = .false.  ! CFL diagnostics written during stage 1; suppress for stages 2 and 3

  ! T* halo exchange: pass Ts1_s(:,:,:,m) as a 3D field — ntr+1 fields total,
  ! well within MAX_DOMAIN_FIELDS regardless of nz.
  do m=1,ntr
    call create_group_pass(pass_Ts1_hprev_s1, Ts1_s(:,:,:,m), G%Domain)
  enddo
  call create_group_pass(pass_Ts1_hprev_s1, hprev_s1, G%Domain)
  call do_group_pass(pass_Ts1_hprev_s1, G%Domain)

  ! Stage 2: reconstruct from T*, compute T** and h**.
  !$OMP parallel do default(shared) private(i,j,m,h_old,h_new,Ihnew_ij,do_ij)
  do k=1,nz
    if (domore_k(k) > 0) then

      call compute_flux_2d(Ts1_s, uhh, vhh, hprev_s1(:,:,k), OBC, ntr, &
          isv, iev, jsv, jev, k, G, GV, local_advect_scheme, &
          flux_xs(:,:,:,k), flux_ys(:,:,:,k), domore_j, uhr, vhr, apply_lim_zs)

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
          hprev_s2(i,j,k) = h_new
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
            Ts2_s(i,j,k,m) = (Reg%Tr(m)%t(i,j,k)*h_old &
                            - 0.25*((flux_x(I,j,m,k) - flux_x(I-1,j,m,k)) &
                            +       (flux_y(i,J,m,k) - flux_y(i,J-1,m,k)))) * Ihnew_ij
        enddo
      enddo ; enddo
    endif ! domore_k
  enddo ! Stage 2 k-loop

  ! T** halo exchange: same pattern — ntr+1 3D fields, one do_group_pass.
  do m=1,ntr
    call create_group_pass(pass_Ts2_hprev_s2, Ts2_s(:,:,:,m), G%Domain)
  enddo
  call create_group_pass(pass_Ts2_hprev_s2, hprev_s2, G%Domain)
  call do_group_pass(pass_Ts2_hprev_s2, G%Domain)

  ! Stage 3: reconstruct from T**, combine fluxes, final update.
  !$OMP parallel do default(shared) private(i,j,m,h_old,h_new,Ihnew_ij,do_ij)
  do k=1,nz
    if (domore_k(k) > 0) then

      call compute_flux_2d(Ts2_s, uhh, vhh, hprev_s2(:,:,k), OBC, ntr, &
        isv, iev, jsv, jev, k, G, GV, local_advect_scheme, &
        flux_xs(:,:,:,k), flux_ys(:,:,:,k), domore_j, uhr, vhr, apply_lim_zs)

      ! Combine: RK3 tracers get (1/6)*(F1+F2) + (2/3)*F3.
      do m=1,ntr
        do j=jsv,jev ; do I=isv-1,iev
          flux_x(I,j,m,k) = (1.0/6.0)*flux_x(I,j,m,k) + (2.0/3.0)*flux_xs(I,j,m,k)
        enddo ; enddo
        do J=jsv-1,jev ; do i=isv,iev
          flux_y(i,J,m,k) = (1.0/6.0)*flux_y(i,J,m,k) + (2.0/3.0)*flux_ys(i,J,m,k)
        enddo ; enddo
      enddo

      ! Final tracer and thickness update — Zhang-Shu-limited WENO fluxes guarantee positivity.
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
!! Self-contained: does not delegate to compute_flux_x_2d / compute_flux_y_2d / zhang_shu_scale.
subroutine compute_flux_2d(tk, uhh_in, vhh_in, h_k, OBC, ntr, &
  is, ie, js, je, k, G, GV, advect_schemes, flux_x_out, flux_y_out, domore_j_k, &
  uhr, vhr, apply_lim)
  type(ocean_grid_type),                          intent(in)    :: G
  type(verticalGrid_type),                        intent(in)    :: GV
  integer,                                        intent(in)    :: ntr
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV),ntr), intent(in)    :: tk
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)),     intent(inout) :: uhh_in   !< zonal mass flux [H L2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)),     intent(inout) :: vhh_in   !< meridional mass flux [H L2]
  real, dimension(SZI_(G),SZJ_(G)),                intent(in)    :: h_k      !< layer thickness [H ~> m or kg m-2]
  type(ocean_OBC_type),                           pointer       :: OBC
  integer,                                        intent(in)    :: is, ie, js, je, k
  integer, dimension(ntr),                        intent(in)    :: advect_schemes
  real, dimension(SZIB_(G),SZJ_(G),ntr),          intent(out)   :: flux_x_out
  real, dimension(SZI_(G),SZJB_(G),ntr),          intent(out)   :: flux_y_out
  logical, dimension(SZJ_(G),SZK_(GV)),           intent(in)    :: domore_j_k
  real, dimension(SZIB_(G),SZJ_(G),SZK_(GV)),     intent(inout) :: uhr   !< accumulated zonal mass flux [H L2]
  real, dimension(SZI_(G),SZJB_(G),SZK_(GV)),     intent(inout) :: vhr   !< accumulated meridional mass flux [H L2]
  logical, dimension(ntr),                        intent(in)    :: apply_lim  !< per-tracer Zhang-Shu enable

  real, dimension(SZI_(G),ntr)          :: Ts2x  !< x-reconstruction stencil (per row)
  real, dimension(SZI_(G),ntr,SZJB_(G)) :: Ts2y  !< y-reconstruction stencil (all rows)
  real, dimension(SZI_(G),SZJ_(G))      :: theta  !< Zhang-Shu upwind scaling factor [nondim]
  real :: order3, order5, order7
  real :: T3(3), T7(7), wq, qext, cfl_face
  real :: hT, outgoing
  integer :: i, j, m, n, i_up, j_up, ntr_id
  type(OBC_segment_type), pointer :: segment=>NULL()

  flux_x_out(:,:,:) = 0.0
  flux_y_out(:,:,:) = 0.0

  ! Y pre-init: fill Ts2y for all j
  do j=G%jsd,G%jed ; do m=1,ntr ; do i=G%isd,G%ied
    Ts2y(i,m,j) = tk(i,j,k,m)
  enddo ; enddo ; enddo
  if (associated(OBC)) then ; if (OBC%OBC_pe) then
    do n=1,OBC%number_of_segments
      segment=>OBC%segment(n)
      if (.not. associated(segment%tr_Reg)) cycle
      do i=is,ie
        if (segment%is_N_or_S .and. i>=segment%HI%isd .and. i<=segment%HI%ied) then
          J = segment%HI%JsdB
          do m=1,segment%tr_Reg%ntseg
            ntr_id = segment%tr_reg%Tr(m)%ntr_index
            if (allocated(segment%tr_Reg%Tr(m)%tres)) then
              if (segment%direction == OBC_DIRECTION_S) then
                Ts2y(i,ntr_id,j) = segment%tr_Reg%Tr(m)%tres(i,j,k)
              else
                Ts2y(i,ntr_id,j+1) = segment%tr_Reg%Tr(m)%tres(i,j,k)
              endif
            else
              if (segment%direction == OBC_DIRECTION_S) then
                Ts2y(i,ntr_id,j) = segment%tr_Reg%Tr(m)%OBC_inflow_conc
              else
                Ts2y(i,ntr_id,j+1) = segment%tr_Reg%Tr(m)%OBC_inflow_conc
              endif
            endif
          enddo
        endif
      enddo
    enddo
  endif ; endif

  ! X reconstruction
  do j=js,je ; if (domore_j_k(j,k)) then
    do m=1,ntr ; do i=G%isd,G%ied ; Ts2x(i,m) = tk(i,j,k,m) ; enddo ; enddo
    if (associated(OBC)) then ; if (OBC%OBC_pe) then
      do n=1,OBC%number_of_segments
        segment=>OBC%segment(n)
        if (.not. associated(segment%tr_Reg)) cycle
        if (segment%is_E_or_W .and. j>=segment%HI%jsd .and. j<=segment%HI%jed) then
          I = segment%HI%IsdB
          do m=1,segment%tr_Reg%ntseg
            ntr_id = segment%tr_reg%Tr(m)%ntr_index
            if (allocated(segment%tr_Reg%Tr(m)%tres)) then
              if (segment%direction == OBC_DIRECTION_W) then
                Ts2x(i,ntr_id) = segment%tr_Reg%Tr(m)%tres(i,j,k)
              else
                Ts2x(i+1,ntr_id) = segment%tr_Reg%Tr(m)%tres(i,j,k)
              endif
            else
              if (segment%direction == OBC_DIRECTION_W) then
                Ts2x(i,ntr_id) = segment%tr_Reg%Tr(m)%OBC_inflow_conc
              else
                Ts2x(i+1,ntr_id) = segment%tr_Reg%Tr(m)%OBC_inflow_conc
              endif
            endif
          enddo
        endif
      enddo
    endif ; endif
    do m=1,ntr
      if ((advect_schemes(m) == ADVECT_WENO5) .or. (advect_schemes(m) == ADVECT_WENO7)) then
        order7 = 0.0
        do I=is-1,ie
          if (uhh_in(I,j,k) >= 0.0) then ; i_up = i ; else ; i_up = i+1 ; endif
          if (uhh_in(I,j,k) > 0.0 .and. h_k(i,j) > 0.0) then
            cfl_face = uhh_in(I,j,k) / h_k(i,j)
          elseif (uhh_in(I,j,k) < 0.0 .and. h_k(i+1,j) > 0.0) then
            cfl_face = -uhh_in(I,j,k) / h_k(i+1,j)
          else
            cfl_face = 0.0
          endif
          T3(:) = Ts2x(i_up-1:i_up+1,m) ; T7(:) = Ts2x(i_up-3:i_up+3,m)
          order3 = G%mask2dCu(I_up-2,j)*G%mask2dCu(I_up-1,j)*G%mask2dCu(I_up,j)*G%mask2dCu(I_up+1,j)
          order5 = order3*G%mask2dCu(I_up-3,j)*G%mask2dCu(I_up+2,j)
          if (advect_schemes(m) == ADVECT_WENO7) &
            order7 = order5*G%mask2dCu(I_up-4,j)*G%mask2dCu(I_up+3,j)
          if (order7 == 1.0) then
            call weno7_reconstruction(wq, T7, uhh_in(I,j,k), cfl_face, G%dxCu(i,J))
          elseif (order5 == 1.0) then
            call weno5_reconstruction(wq, T7, uhh_in(I,j,k), cfl_face, G%dxCu(i,J))
          else
            qext = G%mask2dCu(I_up,j)*G%mask2dCu(I_up-1,j)
            call PPM_reconstruction(wq, T3(1), T3(2), T3(3), uhh_in(I,j,k), cfl_face, qext)
          endif
          flux_x_out(I,j,m) = uhh_in(I,j,k)*wq
        enddo
      endif
    enddo
    if (associated(OBC)) then ; if (OBC%OBC_pe) then
      if (OBC%specified_u_BCs_exist_globally .or. OBC%open_u_BCs_exist_globally) then
        do n=1,OBC%number_of_segments
          segment=>OBC%segment(n)
          if (.not. associated(segment%tr_Reg)) cycle
          if (segment%is_E_or_W .and. j>=segment%HI%jsd .and. j<=segment%HI%jed) then
            I = segment%HI%IsdB
            if ((uhr(I,j,k) > 0.0 .and. segment%direction == OBC_DIRECTION_W) .or. &
              (uhr(I,j,k) < 0.0 .and. segment%direction == OBC_DIRECTION_E)) then
              uhh_in(I,j,k) = uhr(I,j,k)
              do m=1,segment%tr_Reg%ntseg
                ntr_id = segment%tr_reg%Tr(m)%ntr_index
                if (allocated(segment%tr_Reg%Tr(m)%tres)) then
                  flux_x_out(I,j,ntr_id) = uhh_in(I,j,k)*segment%tr_Reg%Tr(m)%tres(I,j,k)
                else
                  flux_x_out(I,j,ntr_id) = uhh_in(I,j,k)*segment%tr_Reg%Tr(m)%OBC_inflow_conc
                endif
              enddo
            endif
          endif
        enddo
      endif
      if (OBC%open_u_BCs_exist_globally) then
        do n=1,OBC%number_of_segments
          segment=>OBC%segment(n)
          I = segment%HI%IsdB
          if (segment%is_E_or_W .and. j>=segment%HI%jsd .and. j<=segment%HI%jed) then
            if (segment%specified) cycle
            if (.not. associated(segment%tr_Reg)) cycle
            if ((uhr(I,j,k) > 0.0 .and. G%mask2dT(i,j) < 0.5) .or. &
              (uhr(I,j,k) < 0.0 .and. G%mask2dT(i+1,j) < 0.5)) then
              uhh_in(I,j,k) = uhr(I,j,k)
              do m=1,segment%tr_Reg%ntseg
                ntr_id = segment%tr_reg%Tr(m)%ntr_index
                if (allocated(segment%tr_Reg%Tr(m)%tres)) then
                  flux_x_out(I,j,ntr_id) = uhh_in(I,j,k)*segment%tr_Reg%Tr(m)%tres(I,j,k)
                else
                  flux_x_out(I,j,ntr_id) = uhh_in(I,j,k)*segment%tr_Reg%Tr(m)%OBC_inflow_conc
                endif
              enddo
            endif
          endif
        enddo
      endif
    endif ; endif
  endif ; enddo ! j-loop

  ! Y reconstruction
  do J=js-1,je
    if (.not. (domore_j_k(J,k) .or. domore_j_k(J+1,k))) cycle
    do m=1,ntr
      if ((advect_schemes(m) == ADVECT_WENO5) .or. (advect_schemes(m) == ADVECT_WENO7)) then
        order7 = 0.0
        do i=is,ie
          if (vhh_in(i,J,k) >= 0.0) then ; j_up = j ; else ; j_up = j+1 ; endif
          if (vhh_in(i,J,k) > 0.0 .and. h_k(i,J) > 0.0) then
            cfl_face = vhh_in(i,J,k) / h_k(i,J)
          elseif (vhh_in(i,J,k) < 0.0 .and. h_k(i,J+1) > 0.0) then
            cfl_face = -vhh_in(i,J,k) / h_k(i,J+1)
          else
            cfl_face = 0.0
          endif
          T3(:) = Ts2y(i,m,j_up-1:j_up+1) ; T7(:) = Ts2y(i,m,j_up-3:j_up+3)
          order3 = G%mask2dCv(i,J_up-2)*G%mask2dCv(i,J_up-1)*G%mask2dCv(i,J_up)*G%mask2dCv(i,J_up+1)
          order5 = order3*G%mask2dCv(i,J_up-3)*G%mask2dCv(i,J_up+2)
          if (advect_schemes(m) == ADVECT_WENO7) &
            order7 = order5*G%mask2dCv(i,J_up-4)*G%mask2dCv(i,J_up+3)
          if (order7 == 1.0) then
            call weno7_reconstruction(wq, T7, vhh_in(i,J,k), cfl_face, G%dyCv(i,J))
          elseif (order5 == 1.0) then
            call weno5_reconstruction(wq, T7, vhh_in(i,J,k), cfl_face, G%dyCv(i,J))
          else
            qext = G%mask2dCv(i,J_up)*G%mask2dCv(i,J_up-1)
            call PPM_reconstruction(wq, T3(1), T3(2), T3(3), vhh_in(i,J,k), cfl_face, qext)
          endif
          flux_y_out(i,J,m) = vhh_in(i,J,k)*wq
        enddo
      endif
    enddo
    if (associated(OBC)) then ; if (OBC%OBC_pe) then
      if (OBC%specified_v_BCs_exist_globally .or. OBC%open_v_BCs_exist_globally) then
        do n=1,OBC%number_of_segments
          segment=>OBC%segment(n)
          if (.not. segment%specified) cycle
          if (.not. associated(segment%tr_Reg)) cycle
          if (OBC%segment(n)%is_N_or_S .and. J>=segment%HI%JsdB .and. J<=segment%HI%JedB) then
            do i=segment%HI%isd,segment%HI%ied
              if ((vhr(i,J,k) > 0.0 .and. segment%direction == OBC_DIRECTION_S) .or. &
                (vhr(i,J,k) < 0.0 .and. segment%direction == OBC_DIRECTION_N)) then
                vhh_in(i,J,k) = vhr(i,J,k)
                do m=1,segment%tr_Reg%ntseg
                  ntr_id = segment%tr_reg%Tr(m)%ntr_index
                  if (allocated(segment%tr_Reg%Tr(m)%tres)) then
                    flux_y_out(i,J,ntr_id) = vhh_in(i,J,k)*OBC%segment(n)%tr_Reg%Tr(m)%tres(i,J,k)
                  else
                    flux_y_out(i,J,ntr_id) = vhh_in(i,J,k)*OBC%segment(n)%tr_Reg%Tr(m)%OBC_inflow_conc
                  endif
                enddo
              endif
            enddo
          endif
        enddo
      endif
      if (OBC%open_v_BCs_exist_globally) then
        do n=1,OBC%number_of_segments
          segment=>OBC%segment(n)
          if (segment%specified) cycle
          if (.not. associated(segment%tr_Reg)) cycle
          if (segment%is_N_or_S .and. J>=segment%HI%JsdB .and. J<=segment%HI%JedB) then
            do i=segment%HI%isd,segment%HI%ied
              if ((vhr(i,J,k) > 0.0 .and. G%mask2dT(i,j) < 0.5) .or. &
                (vhr(i,J,k) < 0.0 .and. G%mask2dT(i,j+1) < 0.5)) then
                vhh_in(i,J,k) = vhr(i,J,k)
                do m=1,segment%tr_Reg%ntseg
                  ntr_id = segment%tr_reg%Tr(m)%ntr_index
                  if (allocated(segment%tr_Reg%Tr(m)%tres)) then
                    flux_y_out(i,J,ntr_id) = vhh_in(i,J,k)*segment%tr_Reg%Tr(m)%tres(i,J,k)
                  else
                    flux_y_out(i,J,ntr_id) = vhh_in(i,J,k)*segment%tr_Reg%Tr(m)%OBC_inflow_conc
                  endif
                enddo
              endif
            enddo
          endif
        enddo
      endif
    endif ; endif
  enddo ! J-loop

  ! Zhang-Shu positivity limiter
  do m = 1, ntr
    if (.not. apply_lim(m)) cycle
    do j = G%jsd, G%jed ; do i = G%isd, G%ied ; theta(i,j) = 1.0 ; enddo ; enddo
    do j = js, je ; if (.not. domore_j_k(j,k)) cycle
      do i = is, ie
        if (h_k(i,j) <= 0.0) cycle
        hT = h_k(i,j) * max(tk(i,j,k,m), 0.0)
        outgoing = max(flux_x_out(I,j,m), 0.0) - min(flux_x_out(I-1,j,m), 0.0) &
                  + max(flux_y_out(i,J,m), 0.0) - min(flux_y_out(i,J-1,m), 0.0)
        if (outgoing > hT) theta(i,j) = hT / outgoing
      enddo
    enddo
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
  enddo ! m

end subroutine compute_flux_2d

!> WENO5-Z + MP reconstruction at the upwind face of the donor cell.
pure subroutine weno5_reconstruction(wq, q, u, cfl, ds)
  real, intent(in)  :: q(7)      !< tracer concentration from cell i-3 to i+3 [conc]
  real, intent(in)  :: u         !< advective velocity [H L2 ~> m3 or kg]
  real, intent(in)  :: cfl    !< absolute value of the advective upwind-cell CFL number [nondim]
  real, intent(in)  :: ds        !< grid spacing [L ~> m]
  real, intent(out) :: wq        !< WENO5 reconstructed face value [conc]

  if (u >= 0.0) then
    call weno5_face(wq, q, cfl, ds, u)          ! left state at i+1/2
  else
    call weno5_face(wq, q(7:1:-1), cfl, ds, u)  ! right state at i-1/2 via mirrored stencil
  endif

end subroutine weno5_reconstruction

!> WENO5-Z + MP reconstruction at the right face of cell q(4), given a 7-point upwind-ordered stencil.
pure subroutine weno5_face(wf, q, cfl, ds, u)
  real, intent(in)  :: q(7)     !< stencil ordered upwind to downwind [conc]
  real, intent(in)  :: cfl   !< absolute value of the advective CFL number [nondim]
  real, intent(in)  :: ds       !< grid spacing [L ~> m]
  real, intent(in)  :: u         !< advective velocity [H L2 ~> m3 or kg]
  real, intent(out) :: wf       !< reconstructed value at the right face of q(4) [conc]

  real :: P0, P1, P2                   ! sub-stencil polynomial reconstructions
  real :: b0, b1, b2                   ! smoothness indicators
  real :: w0, w1, w2                   ! nonlinear weights
  real :: tau, wnorm, bmin             ! WENO-Z indicators
  real :: dm2, dm1, dd0, dd1, dd2      ! second differences
  real :: dm4p, dm4m                   ! 4th-order minmod combinations
  real :: qul, qmp, qmd, qlc           ! MP limiter reference values
  real :: qmin, qmax                   ! monotone range
  real :: q0_min, q0_max, Dqm, Dqp    ! fallback limiter values
  logical :: flag
  real, parameter :: C1_6 = 1.0/6.0
  real, parameter :: d0 = 1.0/10.0, d1 = 6.0/10.0, d2 = 3.0/10.0
  real, parameter :: alpha = 2.0
  real :: wpl, wmr, dA, mA, a6, fc, bM

  ! WENO5-Z sub-stencil reconstructions
  P0 = (2.0*q(2) - 7.0*q(3) + 11.0*q(4)) * C1_6
  b0 = q(2)*(4.0*q(2) - 19.0*q(3) + 11.0*q(4)) + q(3)*(25.0*q(3) - 31.0*q(4)) + 10.0*q(4)*q(4)

  P1 = (-q(3) + 5.0*q(4) + 2.0*q(5)) * C1_6
  b1 = q(3)*(4.0*q(3) - 13.0*q(4) + 5.0*q(5)) + q(4)*(13.0*q(4) - 13.0*q(5)) + 4.0*q(5)*q(5)

  P2 = (2.0*q(4) + 5.0*q(5) - q(6)) * C1_6
  b2 = q(4)*(10.0*q(4) - 31.0*q(5) + 11.0*q(6)) + q(5)*(25.0*q(5) - 19.0*q(6)) + 4.0*q(6)*q(6)

  bmin = min(b0, b1, b2)
  bM = (b0+b1+b2)/3.0

  ! WENO-Z nonlinear weights
  tau = abs(b2 - b0)
  w0 = d0 * weight_fac(tau, b0)
  w1 = d1 * weight_fac(tau, b1)
  w2 = d2 * weight_fac(tau, b2)
  wnorm = 1.0 / (w0 + w1 + w2)
  wf = (w0*P0 + w1*P1 + w2*P2) * wnorm

  ! MP limiter (Suresh & Huynh 1997, bounds from He et al. 2016)
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
  q0_min = min(q(4), qmp) ; q0_max = max(q(4), qmp)

  wf = min(max(qmin, qmax), wf) ; wf = max(min(qmin, qmax), wf)

  ! Near-discontinuity fallback
  fc = maxval(q(:)**2)
  flag = (((tau > fc*(ds**2)) .or. (tau > bM) .or. (bmin > fc*ds)))

  if (((qmax-qmin) > (q0_max-q0_min) .and. flag) .or. (abs(minval(q)) <= 1.0e-5) .or. &
      (wf < 1.0e-5 .and. minval(q(:)) >= 0.0)) then

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

    a6 = 6.*q(4) - 3. * (wpl + wmr) ! Curvature
    wf = (wpl - 0.5 * cfl * ((wpl - wmr) - a6 * (1. - 2./3. * cfl)))
  endif

end subroutine weno5_face

!> WENO7-Z + MP reconstruction at the upwind face of the donor cell.  
pure subroutine weno7_reconstruction(wq, q, u, cfl, ds)
  real, intent(in)  :: q(7)      !< tracer concentration from cell i-3 to i+3 [conc]
  real, intent(in)  :: u         !< advective velocity [H L2 ~> m3 or kg]
  real, intent(in)  :: cfl       !< absolute value of the advective upwind-cell CFL number [nondim]
  real, intent(in)  :: ds        !< grid spacing [L ~> m]
  real, intent(out) :: wq        !< WENO7 reconstructed face value [conc]

  if (u >= 0.0) then
    call weno7_face(wq, q, cfl, ds, u)          ! left state at i+1/2
  else
    call weno7_face(wq, q(7:1:-1), cfl, ds, u)  ! right state at i-1/2
  endif

  end subroutine weno7_reconstruction

!> WENO7-Z + MP5 reconstruction at the right face of cell q(4), given a 7-point upwind-ordered stencil.
pure subroutine weno7_face(wf, q, cfl, ds, u)
  real, intent(in)  :: q(7)     !< stencil ordered upwind to downwind [conc]
  real, intent(in)  :: cfl   !< absolute value of the advective CFL number [nondim]
  real, intent(in)  :: ds       !< grid spacing [L ~> m]
  real, intent(in)  :: u         !< advective velocity [H L2 ~> m3 or kg]
  real, intent(out) :: wf       !< reconstructed value at the right face of q(4) [conc]

  real :: P0, P1, P2, P3               ! sub-stencil polynomial reconstructions
  real :: b0, b1, b2, b3               ! smoothness indicators
  real :: w0, w1, w2, w3               ! nonlinear weights
  real :: tau, wnorm, bmin             ! WENO-Z indicators
  real :: dm2, dm1, dd0, dd1, dd2      ! second differences
  real :: dm4p, dm4m                   ! 4th-order minmod combinations
  real :: qul, qmp, qmd, qlc           ! MP limiter reference values
  real :: qmin, qmax                   ! monotone range
  real :: q0_min, q0_max, Dqm, Dqp    ! fallback limiter values
  logical :: flag
  real, parameter :: C1_12 = 1.0/12.0
  real, parameter :: d0 = 1.0/35.0, d1 = 12.0/35.0, d2 = 18.0/35.0, d3 = 4.0/35.0
  real, parameter :: alpha = 2.0
  real :: wpl, wmr, dA, mA, a6, fc, bM, r, phi

  ! WENO7-Z sub-stencil reconstructions
  P0 = (-3.0*q(1) + 13.0*q(2) - 23.0*q(3) + 25.0*q(4)) * C1_12
  b0 = q(1)*(547.0*q(1) - 3882.0*q(2) + 4642.0*q(3) - 1854.0*q(4)) + &
      q(2)*(7043.0*q(2) - 17246.0*q(3) + 7042.0*q(4)) + &
      q(3)*(11003.0*q(3) - 9402.0*q(4)) + 2107.0*q(4)*q(4)

  P1 = (q(2) - 5.0*q(3) + 13.0*q(4) + 3.0*q(5)) * C1_12
  b1 = q(2)*(267.0*q(2) - 1642.0*q(3) + 1602.0*q(4) - 494.0*q(5)) + &
      q(3)*(2843.0*q(3) - 5966.0*q(4) + 1922.0*q(5)) + &
      q(4)*(3443.0*q(4) - 2522.0*q(5)) + 547.0*q(5)*q(5)

  P2 = (-q(3) + 7.0*q(4) + 7.0*q(5) - q(6)) * C1_12
  b2 = q(3)*(547.0*q(3) - 2522.0*q(4) + 1922.0*q(5) - 494.0*q(6)) + &
      q(4)*(3443.0*q(4) - 5966.0*q(5) + 1602.0*q(6)) + &
      q(5)*(2843.0*q(5) - 1642.0*q(6)) + 267.0*q(6)*q(6)

  P3 = (3.0*q(4) + 13.0*q(5) - 5.0*q(6) + q(7)) * C1_12
  b3 = q(4)*(2107.0*q(4) - 9402.0*q(5) + 7042.0*q(6) - 1854.0*q(7)) + &
      q(5)*(11003.0*q(5) - 17246.0*q(6) + 4642.0*q(7)) + &
      q(6)*(7043.0*q(6) - 3882.0*q(7)) + 547.0*q(7)*q(7)

  bmin = min(b0, b1, b2, b3)
  bM = 0.25*(b0+b1+b2+b3)

  ! WENO7-Z nonlinear weights
  tau = abs((b0 - b3) + 3*(b1 - b2))
  w0 = d0 * weight_fac(tau, b0)
  w1 = d1 * weight_fac(tau, b1)
  w2 = d2 * weight_fac(tau, b2)
  w3 = d3 * weight_fac(tau, b3)
  wnorm = 1.0 / (w0 + w1 + w2 + w3)
  wf = (w0*P0 + w1*P1 + w2*P2 + w3*P3) * wnorm

  ! MP limiter (Suresh & Huynh 1997, bounds from He et al. 2016)
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
  q0_min = min(q(4), qmp) ; q0_max = max(q(4), qmp)

  wf = min(max(qmin, qmax), wf) ; wf = max(min(qmin, qmax), wf)

  ! Near-discontinuity fallback
  fc = maxval(q(:)**2)
  flag = (((tau > fc*(ds**2)) .or. (tau > bM) .or. (bmin > fc*ds)))

  if (((qmax-qmin) > (q0_max-q0_min) .and. flag) .or. (abs(minval(q)) <= 1.0e-5) .or. &
      (wf < 1.0e-5 .and. minval(q(:)) >= 0.0)) then

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

!> 5th-order weno z-type reconstruction flux
pure subroutine ppmw5_reconstruction(wq, q, u, cfl)
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
	real :: dm2, dd2, wpl, wmr, dA, mA, a6, disct
	logical :: lim, disc

	lim = .false.

	! Left state at i+1/2
	P0 = ((2.0*q(2) - 7.0*q(3)) + 11.0*q(4))*C1_6
	b0 = q(2)*((4.0*q(2) - 19.0*q(3)) + 11.0*q(4)) + (q(3)*(25.0*q(3) - 31.0*q(4)) + 10.0*(q(4)*q(4)))

	P1 = ((-q(3) + 5.0*q(4)) + 2.0*q(5))*C1_6
	b1 = q(3)*((4.0*q(3) - 13.0*q(4)) + 5.0*q(5)) + (q(4)*(13.0*q(4) - 13.0*q(5)) + 4.0*(q(5)*q(5)))

	P2 = ((2.0*q(4) + 5.0*q(5)) - q(6))*C1_6
	b2 = q(4)*((10.0*q(4) - 31.0*q(5)) + 11.0*q(6)) + (q(5)*(25.0*q(5) - 19.0*q(6)) + 4.0*(q(6)*q(6)))

	disc = (abs(b2-b0) >= min(b0, b1, b2) .or. &
          (abs(maxval(cfl) - minval(cfl)) > 0.0))

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
  qlc = 0.5*(3.0*q(4) - q(3)) + (4.0/3.0)*dm4m

	qmin = max(min(q(4), q(5), qmd), min(q(4), qul, qlc))
	qmax = min(max(q(4), q(5), qmd), max(q(4), qul, qlc))
	q0_min = min(q(4), qmp) ; q0_max = max(q(4), qmp)
  wpl = min(max(qmin, qmax), wpl) ; wpl = max(min(qmin, qmax), wpl)

	if ((wpl < 1.0e-5 .and. ((minval(q(:)) >= 0.0) )) .or. &
      (((qmax-qmin) > (q0_max-q0_min)) .and. disc )) then
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
				.or. (abs(maxval(cfl) - minval(cfl)) > 0.0))

	! Nonlinear weights
	tau = abs(b2 - b0)
	w0 = d0*weight_fac(tau, b0)
	w1 = d1*weight_fac(tau, b1)
	w2 = d2*weight_fac(tau, b2)

	wnorm = 1.0 / (w0 + w1 + w2)
	wmr = (w0*P0 + w1*P1 + w2*P2) * wnorm

	! MP limiter based on He et al. 2016
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
  qlc = 0.5*(3.0*q(3) - q(4)) + (4.0/3.0)*dm4m

	qmin = max(min(q(4), q(3), qmd), min(q(4), qul, qlc))
	qmax = min(max(q(4), q(3), qmd), max(q(4), qul, qlc))
	q0_min = min(q(4), qmp) ; q0_max = max(q(4), qmp)
  
	if ((wmr < 1.0e-5 .and. ((minval(q(:)) >= 0.0) )) .or. &
      (((qmax-qmin) > (q0_max-q0_min)) .and. disc )) then
		lim = .true.
	endif

	if (lim .or. (abs(minval(q)) <= 1.0e-5)) then
		wpl = min(max(q(4), q(5)), wpl) ; wpl = max(min(q(4), q(5)), wpl)
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

end subroutine ppmw5_reconstruction

!> 7th-order weno z-type reconstruction flux
pure subroutine ppmw7_reconstruction(wq, q, u, cfl)
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
  real :: dm2, dd2, wpl, wmr, dA, mA, a6
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
            .or. (abs(maxval(cfl) - minval(cfl)) > 0.0))

  wpl = (w0*P0 + w1*P1 + w2*P2 + w3*P3) !* wnorm

  ! MP limiter based on He et al. 2016
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
  qlc = 0.5*(3.0*q(4) - q(3)) + (4.0/3.0)*dm4m

  qmin = max(min(q(4), q(5), qmd), min(q(4), qul, qlc))
  qmax = min(max(q(4), q(5), qmd), max(q(4), qul, qlc))
  q0_min = min(q(4), qmp) ; q0_max = max(q(4), qmp)

  if ((wpl < 1.0e-5 .and. ((minval(q(:)) >= 0.0) )) .or. &
    (((qmax-qmin) > (q0_max-q0_min)) .and. disc )) then
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
            .or. (abs(maxval(cfl) - minval(cfl)) > 0.0))
  wmr = (w0*P0 + w1*P1 + w2*P2 + w3*P3) !* wnorm

  ! MP limiter based on He et al. 2016
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
  qlc = 0.5*(3.0*q(3) - q(4)) + (4.0/3.0)*dm4m

  qmin = max(min(q(4), q(3), qmd), min(q(4), qul, qlc))
  qmax = min(max(q(4), q(3), qmd), max(q(4), qul, qlc))
  q0_min = min(q(4), qmp) ; q0_max = max(q(4), qmp)

  if ((wmr < 1.0e-5 .and. ((minval(q(:)) >= 0.0) )) &
    .or. ((wmr-qmin)*(wmr-qmax) > 0.0) .or. &
      (((qmax-qmin) > (q0_max-q0_min)) .and. disc )) then
    lim = .true.
  endif

  if (lim .or. (abs(minval(q)) <= 1.0e-5)) then
    wpl = min(max(q(4), q(5)), wpl) ; wpl = max(min(q(4), q(5)), wpl)
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

end subroutine ppmw7_reconstruction

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

!> Compute the WENO-Z weight factor (1 + (tau/b)^2), protected against b=0.
pure function weight_fac(tau, b) result(factor)
  real, intent(in) :: tau  !< Difference of the smoothness indicator [A ~> a]
  real, intent(in) :: b    !< The smoothness indicator [A ~> a]
  real :: factor

  factor = (1.0 + (tau/(b + 1.0e-40))**2)

end function weight_fac

!> \namespace mom_tracer_advect
!!
!!  WENO5-Z and WENO7-Z tracer reconstruction schemes (Balsara et al. 2016, Borges et al. 2008)
!!  with MP5 monotonicity-preserving limiting (Suresh & Huynh 1997, He et al. 2016).
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

end module MOM_tracer_advect_weno
