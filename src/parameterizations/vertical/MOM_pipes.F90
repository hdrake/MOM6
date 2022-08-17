!> Implemented non-local mixing due to Ocean Thermal Energy Conversion
module MOM_pipes

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

implicit none ; private

#include <MOM_memory.h>

public pipe_mass_and_tracer, find_layer

!> Control structure that stores temperature and salinity as 1-dimensional arrays
!! for a single grid cell.
type, public :: var1d_type ;
  ! If allocated, the following variables have nz layers.
  real, pointer :: h(:) => Null() !< Layer thickness [m].
  real, pointer :: T(:) => NULL() !< Potential temperature [degC].
  real, pointer :: Tr(:) => NULL() !< Salinity [PSU] or [gSalt/kg], generically [ppt], or generic Tracer [conc].
end type var1d_type

contains

!> Given a depth, finds the appropriate layer that contains that depth.
subroutine find_layer(h1d, GV, target_depth, k, layer_depth)
  type(verticalGrid_type),   intent(in)  :: GV !< The ocean's vertical grid structure.
  real, dimension(SZK_(GV)), intent(in)  :: h1d !< Layer thicknesses at the grid cell [H ~> m or kg m-2].
  real,                      intent(in)  :: target_depth !< The depth we are searching for [H ~> m or kg m-2].

  integer, intent(out) :: k !< The layer corresponding to a depth of target_depth.
                            !! If k > GV%ke, then the depth does not exist at this location.
  real,    intent(out) :: layer_depth !< The depth of layer k [H ~> m or kg m-2].

  k = 0
  layer_depth = 0.0

  do while (layer_depth <= target_depth)
    k = k + 1
    if (k > GV%ke) return ! ocean is not deep enough
    layer_depth = layer_depth + h1d(k)
  enddo
end subroutine find_layer

!> Drains the layer at a cell, starting at a minimum depth of z.
!! If this depletes the layer fully, then uses layer k+1 to finish draining.
!! Upon depleting the deepest layer, stops with a warning in stdout.
subroutine mass_and_tracer_sink(v1d, GV, sink_depth, dThickness, &
                                    netMassOut, netHeatOut, netTracerOut)
  type(var1d_type),       intent(inout)    :: v1d !< 1-dimensional copies of h, T, and S
  type(verticalGrid_type),   intent(in)    :: GV !< The ocean's vertical grid structure.

  ! Mass Sink Parameters
  real,                      intent(in)    :: sink_depth   !< The depth of this mass sink [H ~> m or kg m-2].
  real,                      intent(in)    :: dThickness   !< Amount to change layer thickness [H ~> m or kg m-2]
                                                           !! Must be negative.
  real,                      intent(inout) :: netMassOut   !< The total mass being extracted [H ~> m or kg m-2].
  real,                      intent(inout) :: netHeatOut   !< The total heat being extracted [degC H ~> degC m or degC kg m-2].
  real,                      intent(inout) :: netTracerOut !< The total amount of salt being extracted
                                                           !! [ppt H ~> ppt m or ppt kg m-2].


  ! Local variables
  integer :: k
  real    :: layer_depth, dh, maximum_drainage, dh_total

  call find_layer(v1d%h, GV, sink_depth, k, layer_depth)
  if (k > GV%ke) then
    call MOM_error(WARNING, "MOM_pipes: Ocean floor reached before pipe intake depth.")
    return
  endif

  dh_total = dThickness

  do while (dh_total < 0) ! as long as there is still more to take out

    ! The maximum drainage from this layer is everything below sink_depth.
    ! Ensure the layer thickness is always at least Angstrom.
    maximum_drainage = min(v1d%h(k) - GV%Angstrom_H, layer_depth - sink_depth)
    ! Drain as much as specified by input, or until the layer vanishes.
    dh = max(dh_total, -maximum_drainage)
    v1d%h(k) = max(GV%Angstrom_H, v1d%h(k) + dh)
    dh_total = dh_total - dh
    layer_depth = layer_depth - dh ! Layer bottom has moved up

    ! Update tracers for output
    netMassOut = netMassOut - dh
    netHeatOut = netHeatOut - dh*v1d%T(k)
    netTracerOut = netTracerOut - dh*v1d%Tr(k)

    ! Increment for the next iteration
    k = k + 1
    if (k > GV%ke) then ! ocean bottom is reached
      call MOM_error(WARNING, "MOM_pipes: Ocean floor reached before pipe intake depth.")
      return
    endif

    layer_depth = layer_depth + v1d%h(k)
  enddo

end subroutine mass_and_tracer_sink

subroutine mass_and_tracer_source(v1d, GV, source_depth, netMassIn, netHeatIn, netTracerIn)
  type(var1d_type),       intent(inout)    :: v1d !< A structure containing pointers
  type(verticalGrid_type),   intent(in)    :: GV !< The ocean's vertical grid structure.
                                                   !! to any available thermodynamic fields.
  ! Source variables
  real,                      intent(in)    :: source_depth !< The depth of this mass sink [H ~> m or kg m-2].

  real, optional,            intent(in)    :: netMassIn    !< The total mass being added per unit area [H ~> m or kg m-2].
  real, optional,            intent(in)    :: netHeatIn    !< The total heat content of the water being added
                                                           !! [degC H ~> degC m or degC kg m-2].
  real, optional,            intent(in)    :: netTracerIn  !< The total amount of salt being added with the water
                                                           !! [ppt H ~> ppt m or ppt kg m-2].

  ! Local variables
  integer :: k
  real :: layer_depth, &
          oldMass, iNewMass ! Inverse of total mass before/after injection [m2 kg-1].

  ! Find the correct layer to insert.
  k = 0; layer_depth = 0.0
  call find_layer(v1d%h, GV, source_depth, k, layer_depth)
  if (k > GV%ke) then
    call MOM_error(WARNING, "MOM_otec: Ocean floor reached before pipe discharge depth.")
    return
  endif
  
  oldMass = v1d%h(k)
  ! Update mass and tracers of the layer.
  v1d%h(k) = v1d%h(k) + netMassIn
  iNewMass = 1./(v1d%h(k))

  v1d%T(k) = (oldMass*v1d%T(k) + netHeatIn) * iNewMass
  v1d%Tr(k) = (oldMass*v1d%Tr(k) + netTracerIn) * iNewMass

end subroutine mass_and_tracer_source

!> Non-local vertical transport seawater volume and tracers from point sink to point source.
subroutine pipe_mass_and_tracer(i, j, sink_depth, source_depth, pipe_velocity, dt, G, GV, v1d)
  type(var1d_type),                          intent(inout) :: v1d!< 1D column of prognostic variables
  type(ocean_grid_type),                     intent(in)    :: G  !< The ocean's grid structure.
  type(verticalGrid_type),                   intent(in)    :: GV !< The ocean's vertical grid structure.
  integer,                                   intent(in)    :: i  !< Time increment [T ~> s].
  integer,                                   intent(in)    :: j  !< Time increment [T ~> s].
  real,                                      intent(in)    :: sink_depth !< Time increment [T ~> s].
  real,                                      intent(in)    :: source_depth !< Time increment [T ~> s].
  real,                                      intent(in)    :: pipe_velocity !< Time increment [T ~> s].
  real,                                      intent(in)    :: dt !< Time increment [T ~> s].

  ! Local variables
  real :: dh, & ! pipe_velocity integrated over the timestep dt
          dMass, dHeat, dTracer, & ! Tracers being mixed and moved
          depth_tot ! To make sure ocean is deep enough to fit the pipes

  ! Skip profile if ocean depth is deeper than cold intake
  depth_tot = G%bathyT(i,j) + G%Z_ref
  if (depth_tot < max(sink_depth, source_depth)) then
    return
  endif
  
  ! Prepare tracers to be moved between layers
  ! dMass is just a shorthand--we actually are conserving volume
  dMass = 0.0; dHeat = 0.0; dTracer = 0.0
  dh = -pipe_velocity * dt
 
  call mass_and_tracer_sink(v1d, GV, sink_depth, dh, dMass, dHeat, dTracer)
  call mass_and_tracer_source(v1d, GV, source_depth, dMass, dHeat, dTracer)

end subroutine pipe_mass_and_tracer

end module MOM_pipes
