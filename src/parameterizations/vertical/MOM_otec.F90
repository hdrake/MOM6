!> Implemented non-local mixing due to Ocean Thermal Energy Conversion
module MOM_otec

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_alloc
use MOM_diag_mediator, only : register_static_field, time_type, diag_ctrl
use MOM_domains,       only : pass_var
use MOM_error_handler, only : MOM_error, FATAL, WARNING, NOTE
use MOM_file_parser,   only : get_param, log_param, log_version, param_file_type
use MOM_io,            only : MOM_read_data, slasher
use MOM_grid,          only : ocean_grid_type
use MOM_unit_scaling,  only : unit_scale_type
use MOM_variables,     only : thermo_var_ptrs
use MOM_verticalGrid,  only : verticalGrid_type, get_thickness_units
use MOM_EOS,           only : calculate_density, calculate_density_derivs
use MOM_EOS,           only : EOS_type
use MOM_pipes,         only : pipe_mass_and_tracer, find_layer, var1d_type

implicit none ; private

#include <MOM_memory.h>

public otec_diabatic, otec_tracer, otec_init

!> Control structure for OTEC
type, public :: otec_CS ;
  logical :: initialized = .false. !< True if this control structure has been initialized.
  logical :: use_otec, apply_otec_thermo, apply_otec_tracer !< If true, OTEC will be applied.

  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to regulate the timing

  ! OTEC input variables
  real    :: w_cw !< Cold-water pipe velocity [Z T-1 ~> m s-1]
  real    :: w_ww !< Warm-water pipe velocity [Z T-1 ~> m s-1]

  real    :: depth_cold, depth_warm, depth_out !< Pipe depths [m]
  
  integer :: id_otec_intake_dT = -1
  integer :: id_Qcw = -1

end type otec_CS

contains

!> Represents the environmental impact of OTEC by pumping in warm near-surface water (warm intake) and cold deep water (cold intake), potentially extracting some thermal energy from them for power generation, and discharged them at a common outflow depth (effectively mixing them). The grid columns and times at which this is applied, and the parameters that prescribe the pipe parameters, are set according to various OTEC deployment scenario options.
subroutine otec_diabatic(h, tv, dt, G, GV, CS, halo)
  type(ocean_grid_type),                     intent(in)    :: G  !< The ocean's grid structure.
  type(verticalGrid_type),                   intent(in)    :: GV !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: h  !< Layer thicknesses [H ~> m or kg m-2]
  type(thermo_var_ptrs),                     intent(inout) :: tv !< A structure containing pointers
                                                                 !! to any available thermodynamic fields.
  real,                                      intent(in)    :: dt !< Time increment [T ~> s].
  type(otec_CS),                             intent(in)    :: CS !< The control structure returned by
                                                           !! a previous call to
                                                           !! otec_init.
  integer,                         optional, intent(in)    :: halo !< Halo width over which to work

  ! Local variables
  real, dimension(SZI_(G),SZJ_(G)) :: &
    intake_dT, &
    Qcw

  integer :: i, j, k, is, ie, js, je, nz, k2, k_warm, k_cold
  real :: dh_cold, dh_warm, & ! w_cw and w_ww applied over the timestep
          dMass, dSalt, dHeat, & ! Tracers being mixed and moved
          warm_layer_depth, cold_layer_depth, depth_tot ! For comparing stratification against threshold
  real, dimension(SZK_(GV)), target :: h1d, T1d, S1d
  type(var1d_type) :: v1d !< 1-dimensional copy of state variables

  k_warm = 0
  k_cold = 0
  warm_layer_depth = 0.0
  cold_layer_depth = 0.0
  depth_tot = 0.0

  v1d%h => h1d
  v1d%T => T1d
  v1d%Tr => S1d

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  if (present(halo)) then
    is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo
  endif

  if (.not. CS%initialized) call MOM_error(FATAL, "MOM_otec: "//&
         "Module must be initialized before it is used.")

  if (.not.CS%apply_otec_thermo) return

  intake_dT(:,:) = 0.0
  Qcw(:,:) = 0.0

  do j=js,je
    do i=is,ie

      ! Copy this column into a 1D array (for runtime efficiency)
      do k=1,GV%ke
        v1d%h(k) = h(i,j,k)
        v1d%T(k) = tv%T(i,j,k)
        v1d%Tr(k) = tv%S(i,j,k)
      enddo

      ! Only continue if temperature difference is larger than 20ºC (net power production)
      call find_layer(v1d%h, GV, CS%depth_warm, k_warm, warm_layer_depth)
      call find_layer(v1d%h, GV, CS%depth_cold, k_cold, cold_layer_depth)
      
      if ( (warm_layer_depth >= CS%depth_warm) .and. (cold_layer_depth >= CS%depth_cold) ) then
      
        intake_dT(i,j) = T1d(k_warm) - T1d(k_cold)

        if ( (intake_dT(i,j)   >  16.0) ) then

          call pipe_mass_and_tracer(i, j, CS%depth_cold, CS%depth_out, CS%w_cw, dt, G, GV, v1d)
          call pipe_mass_and_tracer(i, j, CS%depth_warm, CS%depth_out, CS%w_ww, dt, G, GV, v1d)

          ! Copy the 1D working arrays back into the original 3D arrays
          do k=1,GV%ke
            h(i,j,k) = h1d(k)
            tv%T(i,j,k) = T1d(k)
            tv%S(i,j,k) = S1d(k)
          enddo

          Qcw(i,j) = CS%w_cw * G%areaT(i,j)

        endif
      endif

    enddo ! i-loop
  enddo ! j-loop

  if (CS%id_otec_intake_dT > 0) then
    call post_data(CS%id_otec_intake_dT, intake_dT, CS%diag)
  endif
  if (CS%id_Qcw > 0) then
    call post_data(CS%id_Qcw, Qcw, CS%diag)
  endif

end subroutine otec_diabatic

subroutine otec_tracer(h, T, Tr, dt, G, GV, CS, halo)
  type(ocean_grid_type),                     intent(in)    :: G  !< The ocean's grid structure.
  type(verticalGrid_type),                   intent(in)    :: GV !< The ocean's vertical grid structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: h  !< Layer thicknesses [H ~> m or kg m-2]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(in)    :: T  !< Temperature [degree C]
  real, dimension(SZI_(G),SZJ_(G),SZK_(GV)), intent(inout) :: Tr !< Tracer concentration on T-cell [conc]
  real,                                      intent(in)    :: dt !< Time increment [T ~> s].
  type(otec_CS),                             intent(in)    :: CS !< The OTEC control structure
  integer,                         optional, intent(in)    :: halo !< Halo width over which to work

  ! Local variables
  integer :: i, j, k, is, ie, js, je, nz, k2, k_warm, k_cold
  real :: layer_depth, depth_tot, deltaT ! For comparing stratification against threshold
  real, dimension(SZK_(GV)), target :: h1d, T1d, Tr1d
  type(var1d_type) :: v1d !< 1-dimensional copy of state variables

  k_warm = 0
  k_cold = 0
  deltaT = 0.0
  layer_depth = 0.0
  depth_tot = 0.0

  v1d%h => h1d
  v1d%T => T1d
  v1d%Tr => Tr1d

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  ! Assume no halo for generic tracers? (Is that correct to do?)
  !if (present(halo)) then
  !  is = G%isc-halo ; ie = G%iec+halo ; js = G%jsc-halo ; je = G%jec+halo
  !endif

  if (.not. CS%initialized) call MOM_error(FATAL, "MOM_otec: "//&
         "Module must be initialized before it is used.")

  if (.not.CS%apply_otec_tracer) return

  do j=js,je
    do i=is,ie

      ! Copy this column into a 1D array (for runtime efficiency)
      do k=1,GV%ke
        v1d%h(k) = h(i,j,k)
        v1d%T(k) = T(i,j,k)
        v1d%Tr(k) = Tr(i,j,k)
      enddo

      ! Only continue if temperature difference is larger than 20ºC (net power production)
      call find_layer(v1d%h, GV, CS%depth_warm, k_warm, layer_depth)
      call find_layer(v1d%h, GV, CS%depth_cold, k_cold, layer_depth)
      deltaT = T1d(k_warm) - T1d(k_cold)

      if (deltaT > 20.0) then

        call pipe_mass_and_tracer(i, j, CS%depth_cold, CS%depth_out, CS%w_cw, dt, G, GV, v1d)
        call pipe_mass_and_tracer(i, j, CS%depth_warm, CS%depth_out, CS%w_ww, dt, G, GV, v1d)

        ! Copy the 1D working arrays back into the original 3D arrays
        do k=1,GV%ke
          Tr(i,j,k) = v1d%Tr(k)
        enddo

      endif
    enddo ! i-loop
  enddo ! j-loop
end subroutine otec_tracer

!> Initialize parameters and allocate memory associated with the OTEC module.
subroutine otec_init(Time, G, GV, param_file, diag, CS)
  type(time_type), target, intent(in)    :: Time !< Current model time.
  type(ocean_grid_type),   intent(inout) :: G    !< The ocean's grid structure.
  type(verticalGrid_type), intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(param_file_type),   intent(in)    :: param_file !< A structure to parse for run-time
                                                 !! parameters.
  type(diag_ctrl), target, intent(inout) :: diag !< Structure used to regulate diagnostic output.
  type(otec_CS),     intent(inout)       :: CS   !< OTEC heating control struct

! This include declares and sets the variable "version".
#include "version_variable.h"
  character(len=40)  :: mdl = "MOM_otec"  ! module name
  character(len=48)  :: thickness_units
  ! Local variables
  character(len=200) :: inputdir, otec_file, filename, otec_var
  logical :: use_otec, apply_otec_thermo, apply_otec_tracer
  real :: w_cw  ! A uniform pumping rate [Z T-1 ~> m s-1]
  real :: gamma ! Ratio of warm- to cold-water pumping rate
  real :: depth_cold, depth_warm, depth_out ! Depths of the pipes [m]
  integer :: i, j, isd, ied, jsd, jed, id
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  CS%initialized = .true.
  CS%diag => diag
  CS%Time => Time

  ! write parameters to the model log.
  call log_version(param_file, mdl, version, "")
  call get_param(param_file, mdl, "USE_OTEC", use_otec, &
                 "Whether or not to use the OTEC module", &
                 default=.false.)
  call get_param(param_file, mdl, "APPLY_OTEC_THERMO", apply_otec_thermo, &
                 "Whether to apply thermodynamic tendencies due to OTEC", &
                 default=.true.)
  call get_param(param_file, mdl, "APPLY_OTEC_TRACER", apply_otec_tracer, &
                 "Whether to apply passive tracer tendencies due to OTEC", &
                 default=.true.)
  call get_param(param_file, mdl, "OTEC_W_CW", w_cw, &
                 "The constant OTEC cold-water pumping rate or 0 to "//&
                 "disable OTEC.", &
                 units="m s-1", default=0.0)
  call get_param(param_file, mdl, "OTEC_GAMMA", gamma, &
                 "The ratio of warm- to cold-water OTEC pumping rates.",&
                 units="nondim", default=1.0)
  call get_param(param_file, mdl, "OTEC_COLD_DEPTH", depth_cold, &
                  "The depth of the cold-water intake for OTEC.", &
                  units="m", default=1000.0)
  call get_param(param_file, mdl, "OTEC_WARM_DEPTH", depth_warm, &
                  "The depth of the warm-water intake for OTEC.", &
                  units="m", default=20.0)
  call get_param(param_file, mdl, "OTEC_OUT_DEPTH", depth_out, &
                  "The depth of the mixed-water output for OTEC.", &
                  units="m", default=500.0)

  CS%w_cw = w_cw
  CS%w_ww = gamma*w_cw

  CS%depth_cold = depth_cold
  CS%depth_warm = depth_warm
  CS%depth_out  = depth_out

  CS%use_otec = use_otec
  CS%apply_otec_thermo = apply_otec_thermo
  CS%apply_otec_tracer = apply_otec_tracer

  ! Post 2D otec diagnostics
  CS%id_otec_intake_dT=register_diag_field('ocean_model', &
    'otec_intake_dT', diag%axesT1, Time,             &
    'Temperature difference between shallow and deep intake pipes', &
    'degC')
  CS%id_Qcw=register_diag_field('ocean_model', &
    'Qcw', diag%axesT1, Time,                  &
    'Cold water pipe transport')

end subroutine otec_init

end module MOM_otec
