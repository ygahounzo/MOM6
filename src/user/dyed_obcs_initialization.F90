! This file is part of MOM6, the Modular Ocean Model version 6.
! See the LICENSE file for licensing information.
! SPDX-License-Identifier: Apache-2.0

!> Dyed open boundary conditions; OBC_USER_CONFIG="dyed_obcs"
module dyed_obcs_initialization

use MOM_dyn_horgrid,     only : dyn_horgrid_type
use MOM_error_handler,   only : MOM_mesg, MOM_error, FATAL, WARNING, is_root_pe
use MOM_file_parser,     only : get_param, log_version, param_file_type
use MOM_get_input,       only : directories
use MOM_grid,            only : ocean_grid_type
use MOM_open_boundary,   only : ocean_OBC_type, OBC_NONE, OBC_DIRECTION_E
use MOM_open_boundary,   only : OBC_segment_type, register_segment_tracer
use MOM_tracer_registry, only : tracer_registry_type, tracer_name_lookup
use MOM_tracer_registry, only : tracer_type
use MOM_variables,       only : thermo_var_ptrs
use MOM_verticalGrid,    only : verticalGrid_type

implicit none ; private

#include <MOM_memory.h>

public dyed_obcs_set_OBC_data

integer :: ntr = 0 !< Number of dye tracers
                   !! \todo This is a module variable. Move this variable into the control structure.
real :: dye_obc_inflow = 0.0 !< Inflow value of obc dye concentration

contains

!> This subroutine sets the dye properties at open boundary conditions.
subroutine dyed_obcs_set_OBC_data(OBC, G, GV, param_file, tr_Reg)
  type(ocean_OBC_type),       pointer    :: OBC !< This open boundary condition type specifies
                                                !! whether, where, and what open boundary
                                                !! conditions are used.
  type(ocean_grid_type),      intent(in) :: G   !< The ocean's grid structure.
  type(verticalGrid_type),    intent(in) :: GV  !< The ocean's vertical grid structure.
  type(param_file_type),      intent(in) :: param_file !< A structure indicating the open file
                                                !! to parse for model parameter values.
  type(tracer_registry_type), pointer    :: tr_Reg !< Tracer registry.

  ! Local variables
  character(len=40)  :: mdl = "dyed_obcs_set_OBC_data" ! This subroutine's name.
  character(len=80)  :: name, longname
  integer :: is, ie, js, je, isd, ied, jsd, jed, m, n, nz, ntr_id
  integer :: IsdB, IedB, JsdB, JedB
  integer :: I, j, k, ntseg
  integer :: n_dye ! Number of regionsl dye tracers
  integer :: dye_obc_k_min ! First k (1=surface) at which OBC dye is applied
  real :: dye ! Inflow dye concentration [arbitrary]
  real :: dye_east_lat_max ! Paint dye_01 on eastern OBC south of this latitude [degrees_N]
  type(tracer_type), pointer      :: tr_ptr => NULL()
  type(OBC_segment_type), pointer :: segment => NULL()

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  if (.not.associated(OBC)) return

  call get_param(param_file, mdl, "NUM_DYED_TRACERS", ntr, &
                 "The number of dyed_obc tracers in this run. Each tracer "//&
                 "should have a separate boundary segment.  "//&
                 "If not present, use NUM_DYE_TRACERS.", default=-1, do_not_log=.true.)
  if (ntr == -1) then
    !for backward compatibility
    call get_param(param_file, mdl, "NUM_DYE_TRACERS", ntr, &
                   "The number of dye tracers in this run. Each tracer "//&
                   "should have a separate boundary segment.", default=0, do_not_log=.true.)
    n_dye = 0
  else
    call get_param(param_file, mdl, "NUM_DYE_TRACERS", n_dye, &
                   "The number of dye tracers in this run. Each tracer "//&
                   "should have a separate region.", default=0, do_not_log=.true.)
  endif

  call get_param(param_file, mdl, "DYE_OBC_INFLOW", dye_obc_inflow, &
                 "The OBC inflow value of dye tracers.", units="kg kg-1", &
                 default=1.0)
  call get_param(param_file, mdl, "DYE_OBC_EAST_LAT_MAX", dye_east_lat_max, &
                 "If greater than this default, dye_01 is also applied on eastern "//&
                 "OBC inflow south of this latitude (Cuba coast to the southern "//&
                 "edge).  The rest of the east wall (Florida Straits, Atlantic) stays 0.  "//&
                 "Unused when left at the default.", units="degrees_N", &
                 default=-1.0e30)
  call get_param(param_file, mdl, "DYE_OBC_K_MIN", dye_obc_k_min, &
                 "First vertical index (1 = surface) at which OBC dye is applied. "//&
                 "Layers 1:DYE_OBC_K_MIN-1 stay 0.  Default 1 paints every layer.  "//&
                 "On the 41-layer GOM hybrid grid, 17 is the first layer whose top "//&
                 "is at/below ~100 m in deep water.", default=1)
  if (dye_obc_k_min < 1) dye_obc_k_min = 1
  if (dye_obc_k_min > nz) dye_obc_k_min = nz

  if (OBC%number_of_segments < ntr) then
    call MOM_error(WARNING, "Error in dyed_obc segment setup")
    return   !!! Need a better error message here
  endif

! ! Set the inflow values of the dyes, one per segment.
! ! We know the order: north, south, east, west
  do m=1,ntr
    write(name,'("dye_",I2.2)') m+n_dye  !after regional dye tracers
    write(longname,'("Concentration of dyed_obc Tracer ",I2.2, " on segment ",I2.2)') m, m
    call tracer_name_lookup(tr_Reg, ntr_id, tr_ptr, name)

    do n=1,OBC%number_of_segments
      if (n == m) then
        dye = dye_obc_inflow
      else
        dye = 0.0
      endif
      call register_segment_tracer(tr_ptr, ntr_id, param_file, GV, &
                                   OBC%segment(n), OBC_scalar=dye)

      ! Restrict OBC dye to k >= DYE_OBC_K_MIN (default 1 = all layers).
      if (dye_obc_k_min > 1) then
        segment => OBC%segment(n)
        if (segment%on_pe) then
          ntseg = segment%tr_Reg%ntseg
          do k=1, dye_obc_k_min-1
            segment%tr_Reg%Tr(ntseg)%t(:,:,k) = 0.0
            segment%tr_Reg%Tr(ntseg)%tres(:,:,k) = 0.0
          enddo
        endif
      endif

      ! Optional: dye_01 on the eastern wall only south of Cuba.
      if ((m == 1) .and. (dye_east_lat_max > -1.0e29)) then
        segment => OBC%segment(n)
        if (segment%on_pe .and. (segment%direction == OBC_DIRECTION_E)) then
          ntseg = segment%tr_Reg%ntseg
          I = segment%Is_obc
          if ((I >= G%IsdB) .and. (I <= G%IedB)) then
            do k=dye_obc_k_min,nz
              do j=max(segment%HI%jsd, G%jsd), min(segment%HI%jed, G%jed)
                if (G%geoLatCu(I,j) <= dye_east_lat_max) then
                  segment%tr_Reg%Tr(ntseg)%t(I,j,k) = dye_obc_inflow
                  segment%tr_Reg%Tr(ntseg)%tres(I,j,k) = dye_obc_inflow
                endif
              enddo
            enddo
          endif
        endif
      endif
    enddo
  enddo

end subroutine dyed_obcs_set_OBC_data

!> \namespace dyed_obcs_initialization
!!
!! Setting dyes, one for painting the inflow on each side.
end module dyed_obcs_initialization
