#include "os-config.h"
#include "os-preprocess.fpp"

module m_zpulse_flyfoc

#include "memory/memory.h"

use m_parameters
use m_zpulse_std
use m_zpulse_mov_wall
use m_node_conf

use m_node_conf
use m_space
use m_vdf_define
use m_grid_define

use m_fparser
use m_math
use m_fft_new

use m_emf_define, only: t_emf_bound, t_emf


implicit none


integer, parameter :: p_gaussian_new = 4 
integer, parameter :: p_laguerre_new = 5 
integer, parameter :: p_astrl = 6
integer, parameter :: p_astrl_laguerre = 7
integer, parameter :: p_astrl_discrete = 8 

integer, parameter :: p_astrl_discrete_num_pulses_max = 32 

!-----------------------------------------------------------------------------------------
type, extends( t_zpulse_mov_wall ) :: t_zpulse_flyfoc

    ! whether we should precompute the envelope and use this for injection 
    logical :: if_interp_env

    ! shared flying focus options 
    integer        :: env_type
    real(p_double) :: interp_dx2, interp_dt 
    real(p_double) :: vff     

    ! Optional laboratory geometry.  This decouples the carrier
    ! propagation direction from the flying-focus trajectory.
    ! The current input interface retains carrier_angle for the validated
    ! 2D geometry, while all stored geometry is represented by 3-vectors
    ! so the same state can be used by the forthcoming 3D implementation.
    logical :: use_lab_geometry
    real(p_double) :: carrier_angle
    real(p_double), dimension(3) :: carrier_dir
    real(p_double), dimension(3) :: lab_basis_p, lab_basis_s
    real(p_double), dimension(3) :: focus_start, focus_end, focus_dir
    real(p_double) :: focus_length, focus_t_start

    ! astrl params
    real(p_double) :: z1, z2, z_lead, zr, initial_delay, final_delay, vg
    integer :: integral_size
    integer :: focal_weight_type
    integer :: focal_delay_type
    integer :: focal_phase_type
    integer :: focal_xperp_type
    type(t_fparser) :: focal_weight_math_func
    type(t_fparser) :: focal_delay_math_func
    type(t_fparser) :: focal_phase_math_func
    type(t_fparser) :: focal_xperp_math_func
    real(p_double) :: env_norm
    real(p_double) :: dz0

    ! discrete astrl params
    integer :: astrl_discrete_num_pulses
    real(p_double), dimension(:), pointer :: astrl_discrete_focal_points => null()  
    real(p_double), dimension(:), pointer :: astrl_discrete_focal_delays => null()  
    real(p_double), dimension(:), pointer :: astrl_discrete_focal_weights => null()  
    real(p_double), dimension(:), pointer :: astrl_discrete_focal_phases => null()  
    real(p_double), dimension(:), pointer :: astrl_discrete_focal_xperps => null()  
    integer :: astrl_discrete_norm_idx 

    complex(p_double), dimension(:,:,:), pointer :: e_env_buffer_2d => null()
    complex(p_double), dimension(:,:,:,:), pointer :: e_env_buffer_3d => null()
  
    ! temporary computation buffer for astrl integral 
    complex(p_double), dimension(:), pointer :: buf_1d => null()

    ! variables used for the prepcomputed envelope for any pulse type
    real(p_double), dimension(4,2) :: e_env_grid_bounds 
    integer, dimension(4) :: e_env_grid_shape
    real(p_double), dimension(4) :: e_env_grid_dx

contains

  procedure :: read_input => read_input_zpulse_flyfoc
  procedure :: init_env_norm => init_env_norm_zpulse_flyfoc
  procedure :: init_extra_vars => init_extra_vars_zpulse_flyfoc
  procedure :: cleanup    => cleanup_zpulse_flyfoc 

  procedure :: t_duration => t_duration_zpulse_flyfoc
  procedure :: t_envelope => t_env_zpulse_flyfoc 

  procedure :: get_env_2d => get_env_2d_zpulse_flyfoc
  procedure :: get_env_3d => get_env_3d_zpulse_flyfoc

  procedure :: get_env_gaussian_2d => get_env_gaussian_2d_zpulse_flyfoc 
  procedure :: get_env_astrl_2d => get_env_astrl_2d_zpulse_flyfoc 
  procedure :: get_env_lab_2d => get_env_lab_2d_zpulse_flyfoc
  procedure :: get_field_lab_2d => get_field_lab_2d_zpulse_flyfoc
  procedure :: get_env_lab_3d => get_env_lab_3d_zpulse_flyfoc
  procedure :: get_field_lab_3d => get_field_lab_3d_zpulse_flyfoc
  procedure :: get_env_astrl_discrete_2d => get_env_astrl_discrete_2d_zpulse_flyfoc 

  procedure :: get_env_gaussian_3d => get_env_gaussian_3d_zpulse_flyfoc 
  procedure :: get_env_laguerre_3d => get_env_laguerre_3d_zpulse_flyfoc 
  procedure :: get_env_astrl_3d => get_env_astrl_3d_zpulse_flyfoc 
  procedure :: get_env_astrl_laguerre_3d => get_env_astrl_laguerre_3d_zpulse_flyfoc
  procedure :: get_env_astrl_discrete_3d => get_env_astrl_discrete_3d_zpulse_flyfoc 

  procedure :: precompute_envelope_2d => precompute_envelope_2d_zpulse_flyfoc
  procedure :: precompute_envelope_3d => precompute_envelope_3d_zpulse_flyfoc
  
  procedure :: launch_wall_2d => launch_wall_2d_zpulse_flyfoc 
  procedure :: launch_wall_3d => launch_wall_3d_zpulse_flyfoc 

  procedure :: launch_wall_t0_2d => launch_wall_t0_2d_zpulse_flyfoc 
  procedure :: launch_lab_t0_2d => launch_lab_t0_2d_zpulse_flyfoc
  procedure :: launch_lab_t0_3d => launch_lab_t0_3d_zpulse_flyfoc


end type t_zpulse_flyfoc
!-----------------------------------------------------------------------------------------



public :: t_zpulse_flyfoc

contains

!-----------------------------------------------------------------------------------------
! Integrate uniformly sampled focal-plane contributions over [z1,z2].
! The original implementation treated integral_size as both the number of
! samples and the number of sub-intervals, so z2 was never sampled.  It also
! applied Simpson 1/3 to an even number of points.  Here integral_size is the
! number of sample points, including both endpoints.  Odd point counts use
! Simpson 1/3 throughout; even point counts use Simpson 1/3 followed by a
! Simpson 3/8 panel for the final three intervals.
complex(p_double) function integrate_astrl_samples( buf, n, h ) result(val)

  complex(p_double), dimension(:), intent(in) :: buf
  integer, intent(in) :: n
  real(p_double), intent(in) :: h
  integer :: m

  val = cmplx( 0.0_p_double, 0.0_p_double, kind=p_double )

  if ( n <= 1 ) then
    return
  else if ( n == 2 ) then
    val = 0.5_p_double * h * ( buf(1) + buf(2) )
  else if ( mod(n,2) == 1 ) then
    val = h / 3.0_p_double * ( buf(1) + buf(n) + &
          4.0_p_double * sum(buf(2:n-1:2)) + &
          2.0_p_double * sum(buf(3:n-2:2)) )
  else if ( n == 4 ) then
    val = 3.0_p_double * h / 8.0_p_double * &
          ( buf(1) + 3.0_p_double*buf(2) + 3.0_p_double*buf(3) + buf(4) )
  else
    ! n is even.  Use Simpson 1/3 through point n-3 (an even number
    ! of sub-intervals), then Simpson 3/8 on points n-3:n.
    m = n - 3
    val = h / 3.0_p_double * ( buf(1) + buf(m) + &
          4.0_p_double * sum(buf(2:m-1:2)) + &
          2.0_p_double * sum(buf(3:m-2:2)) )
    val = val + 3.0_p_double * h / 8.0_p_double * &
          ( buf(n-3) + 3.0_p_double*buf(n-2) + &
            3.0_p_double*buf(n-1) + buf(n) )
  endif

end function integrate_astrl_samples
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine read_input_zpulse_flyfoc( this, input_file, g_space, bnd_con, periodic, grid, sim_options )

  use m_input_file

  implicit none

  class( t_zpulse_flyfoc ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  type( t_space ), intent(in) :: g_space
  class (t_emf_bound), intent(in) :: bnd_con
  logical, dimension(:), intent(in) :: periodic
  class( t_grid ), intent(in) :: grid
  type( t_options ), intent(in) :: sim_options

  ! local variables
  real(p_double) :: a0
  real(p_double) :: omega0
  real(p_double) :: phase
  integer        :: pol_type
  real(p_double) :: pol
  character(len=16) :: propagation
  integer         :: direction

  ! Temporal profile
  character(len=16) :: tenv_type
  real(p_double) :: tenv_fwhm
  real(p_double) :: tenv_rise, tenv_flat, tenv_fall
  real(p_double) :: tenv_duration, tenv_range
  character(len = p_max_expr_len) :: tenv_math_func
  real(p_double)  :: xi0

  ! shared flying focus params
  character(len=32) :: env_type 
  logical :: if_interp_env
  logical :: if_launch_t0
  real(p_double) :: interp_dx2, interp_dt 
  real(p_double) :: vff 
  logical :: use_lab_geometry
  real(p_double) :: carrier_angle, focus_t_start
  real(p_double), dimension(3) :: carrier_dir, focus_start, focus_end

  ! astrl params
  real(p_double) :: z1, z2, initial_delay, final_delay, vg
  integer :: integral_size
  character(len=16) :: focal_weight_type  
  character(len=16) :: focal_delay_type  
  character(len=16) :: focal_phase_type  
  character(len=16) :: focal_xperp_type  
  character(len = p_max_expr_len) :: focal_weight_math_func
  character(len = p_max_expr_len) :: focal_delay_math_func
  character(len = p_max_expr_len) :: focal_phase_math_func
  character(len = p_max_expr_len) :: focal_xperp_math_func

  ! discrete astrl params
  integer :: astrl_discrete_num_pulses
  real(p_double), dimension(p_astrl_discrete_num_pulses_max) :: astrl_discrete_focal_points
  real(p_double), dimension(p_astrl_discrete_num_pulses_max) :: astrl_discrete_focal_delays
  real(p_double), dimension(p_astrl_discrete_num_pulses_max) :: astrl_discrete_focal_weights
  real(p_double), dimension(p_astrl_discrete_num_pulses_max) :: astrl_discrete_focal_phases
  real(p_double), dimension(p_astrl_discrete_num_pulses_max) :: astrl_discrete_focal_xperps
  logical :: astrl_discrete_use_linspace
  integer :: astrl_discrete_norm_idx

  ! launch parameters
  real(p_double) :: launch_time
  logical        :: if_launch

  ! Beam profile
  real(p_double), dimension(2) :: per_center
  real(p_double), dimension(2) :: per_w0

  ! wall params
  real(p_double)  :: wall_pos  ! initial wall position
  real(p_double)  :: wall_vel  
  integer         :: ncells 
  logical         :: if_launch_from_wall
  logical         :: if_clean_fields
  logical         :: if_apply_divcorr
  logical         :: if_add_fields

  ! tmp variables
  integer :: N, i 
  real(p_k_fparse), dimension(1) :: fparser_arr 
  real(p_double) :: carrier_nmag, basis_pmag
  real(p_double), dimension(3) :: basis_ref
  logical :: carrier_dir_any_set, carrier_dir_all_set

  namelist /nl_zpulse_flyfoc/ &
        a0, omega0, phase, pol_type, pol, propagation, direction, &
        tenv_type, tenv_fwhm, tenv_rise, tenv_flat, tenv_fall, &
        tenv_duration, tenv_range, & ! tenv_launch_duration, 
        tenv_math_func, xi0, &
        per_center, per_w0, &
        env_type, &
        if_interp_env, if_launch_t0, interp_dx2, interp_dt, & 
        use_lab_geometry, carrier_angle, carrier_dir, focus_start, focus_end, focus_t_start, &
        vff, vg, z1, z2, initial_delay, final_delay, integral_size, &
        focal_weight_type, focal_weight_math_func, &
        focal_delay_type,  focal_delay_math_func, &
        focal_phase_type,  focal_phase_math_func, &
        focal_xperp_type,  focal_xperp_math_func, &
        astrl_discrete_num_pulses, &
        astrl_discrete_focal_delays, astrl_discrete_focal_points, &
        astrl_discrete_focal_weights, astrl_discrete_focal_phases, & 
        astrl_discrete_focal_xperps, &
        astrl_discrete_use_linspace, astrl_discrete_norm_idx, & 
        wall_pos, wall_vel, ncells, if_launch_from_wall, if_clean_fields, &
        if_apply_divcorr, if_add_fields, launch_time, if_launch


  integer :: ierr

  if (disp_out(input_file)) then
    SCR_ROOT(" - Reading zpulse_flyfoc configuration...")
  endif

  if_launch = .true.

  a0             = 1.0_p_double
  omega0         = 10.0_p_double
  phase          = 0.0_p_double
  pol_type       = 0
  pol            = 90.0_p_double

  propagation    = "forward"
  direction      = 1

  tenv_type       = "gaussian"
  tenv_fwhm       = -1.0_p_double
  tenv_rise       = 0.0_p_double
  tenv_flat       = 0.0_p_double
  tenv_fall       = 0.0_p_double
  xi0             = 0.0_p_double

  tenv_duration   = 0.0_p_double
  tenv_range      = -huge(1.0_p_double)
  ! tenv_launch_duration      = -huge(1.0_p_double)
  tenv_math_func  = "NO_FUNCTION_SUPPLIED!"

  per_center     = -huge(1.0_p_double)

  ! shared ff options 
  ! focal_length           = 0.0_p_double
  env_type             = "" 
  if_interp_env            = .true.
  if_launch_t0          = .false.
  interp_dx2             = -1.0_p_double
  interp_dt              = -1.0_p_double 

  ! Laboratory-geometry defaults: disabled, so legacy decks are unchanged.
  use_lab_geometry = .false.
  carrier_angle = 0.0_p_double
  carrier_dir(:) = huge(1.0_p_double)
  focus_start(:) = 0.0_p_double
  focus_end(:) = 0.0_p_double
  focus_t_start = 0.0_p_double

  ! continuous astrl
  z1 = 0.0_p_double
  z2 = 0.0_p_double
  initial_delay = 0.0_p_double
  final_delay = 0.0_p_double
  vg = 1.0_p_double  ! group velocity 
  integral_size = 1
  
  focal_weight_type  = "const"
  focal_weight_math_func  = "NO_FUNCTION_SUPPLIED!"
  focal_delay_type  = "const"
  focal_delay_math_func  = "NO_FUNCTION_SUPPLIED!"
  focal_phase_type  = "const"
  focal_phase_math_func  = "NO_FUNCTION_SUPPLIED!"
  focal_xperp_type  = "const"
  focal_xperp_math_func  = "NO_FUNCTION_SUPPLIED!"

  ! discrete astrl 
  astrl_discrete_num_pulses = 1 
  astrl_discrete_focal_delays(:) = 0.0_p_double 
  astrl_discrete_focal_points(:) = 0.0_p_double 
  astrl_discrete_focal_weights(:) = 1.0_p_double 
  astrl_discrete_focal_phases(:) = 0.0_p_double 
  astrl_discrete_focal_xperps(:) = 0.0_p_double 
  astrl_discrete_use_linspace = .true.
  astrl_discrete_norm_idx = 1 

  ! general options 
  launch_time           = 0.0_p_double
  if_launch             = .true.

  ! moving wall params
  wall_vel = 0.0_p_double
  wall_pos = -huge(1.0d0)
  ncells = 4
  if_launch_from_wall = .false.
  if_clean_fields = .false.
  if_apply_divcorr = .true.
  if_add_fields = .false.

  if( p_x_dim .eq. 1 ) then 
    SCR_ROOT(   "ERROR: zpulse_flyfoc only works for spatial dimensions > 1")
    stop 
  endif

  ! read data from file
  call get_namelist( input_file, "nl_zpulse_flyfoc", ierr )

  if ( ierr /= 0 ) then
    if ( mpi_node() == 0 ) then
      if (ierr < 0) then
        SCR_ROOT( "Error reading zpulse_flyfoc parameters" )
      else
        SCR_ROOT( "Error: zpulse_flyfoc parameters missing" )
      endif
      SCR_ROOT( "aborting...")
    endif
    stop
  endif

  read (input_file%nml_text, nml = nl_zpulse_flyfoc, iostat = ierr)
  if (ierr /= 0) then
    if ( mpi_node() == 0 ) then
      write(0,*)  ""
      write(0,*)  "   Error reading zpulse_flyfoc parameters"
      write(0,*)  "   aborting..."
    endif
    stop
  endif

  this%a0         = a0
  this%omega0     = omega0
  this%phase0     = real( phase * pi_180, p_double )

  ! Polarization parameters
  if ((pol_type <-1) .or. (pol_type >1)) then
    if ( mpi_node() == 0 ) then
      write(0,*)  ""
      write(0,*)  "   Error in zpulse_flyfoc parameters"
      write(0,*)  "   pol must be in the range [-1,0,+1]"
      write(0,*)  "   aborting..."
    endif
    stop
  endif
  this%pol_type  = pol_type
  this%pol = real( pol * pi_180, p_double )


  ! Propagation direction
  select case ( trim(propagation) )
  case ("forward")
    this%propagation = p_forward
  case ("backward")
    this%propagation = p_backward
  case default
    if ( mpi_node() == 0 ) then
      write(0,*)  ''
      write(0,*)  '   Error in zpulse_flyfoc parameters'
      write(0,*)  '   propagation must be either "forward" or "backward"'
      write(0,*)  '   aborting...'
    endif
    stop
  end select

  ! Launch direction
  if ((direction < 1) .or. (direction > p_x_dim)) then
    if ( mpi_node() == 0 ) then
      write(0,*)  ""
      write(0,*)  "   Error in zpulse_flyfoc parameters"
      write(0,*)  "   direction must be in the range [1 .. x_dim]"
      write(0,*)  "   aborting..."
    endif
    stop
  endif

  this%direction = direction

  ! Temporal envelope parameters
  select case ( trim(tenv_type) )
  case ("polynomial")
    this%tenv_type = p_polynomial

    if ( tenv_fwhm > 0.0_p_double ) then
      this%tenv_rise       = tenv_fwhm/2.0_p_double
      this%tenv_flat       = 0.0_p_double
      this%tenv_fall       = tenv_fwhm/2.0_p_double
    else
      this%tenv_rise       = tenv_rise
      this%tenv_flat       = tenv_flat
      this%tenv_fall       = tenv_fall
    endif

  case ("sin2")
    this%tenv_type = p_sin2

    if ( tenv_fwhm > 0.0_p_double ) then
      this%tenv_rise       = tenv_fwhm/2.0_p_double
      this%tenv_flat       = 0.0_p_double
      this%tenv_fall       = tenv_fwhm/2.0_p_double
    else
      this%tenv_rise       = tenv_rise
      this%tenv_flat       = tenv_flat
      this%tenv_fall       = tenv_fall
    endif

  case ("gaussian")
    this%tenv_type = p_gaussian

    if ( tenv_fwhm > 0.0_p_double  ) then
      this%tenv_duration   = tenv_fwhm/sqrt( 2.0_p_double * log(2.0_p_double))
    else
      this%tenv_duration   = tenv_duration
    endif
    this%tenv_range      = tenv_range

    if ( tenv_range == -huge(1.0_p_double) ) then
      if ( mpi_node() == 0 ) then
        write(0,*)  ''
        write(0,*)  '   Error in zpulse_flyfoc parameters'
        write(0,*)  '   When using "gaussian" temporal envelope, '
        write(0,*)  '   tenv_range must also be defined'
        write(0,*)  '   aborting...'
      endif
      stop
    endif

  case ("math")
    this%tenv_type = p_func
    this%tenv_range      = tenv_range

    if ( tenv_range == -huge(1.0_p_double) ) then
      if ( mpi_node() == 0 ) then
        write(0,*)  ''
        write(0,*)  '   Error in zpulse_flyfoc parameters'
        write(0,*)  '   When using "math" temporal envelopes, '
        write(0,*)  '   tenv_range must also be defined'
        write(0,*)  '   aborting...'
      endif
      stop
    endif

    call setup(this%tenv_math_func, trim(tenv_math_func), (/'t'/), ierr)
    if (ierr /= 0) then
      if ( mpi_node() == 0 ) then
        write(0,*)  ''
        write(0,*)  '   Error in zpulse_flyfoc parameters'
        write(0,*)  '   Supplied tenv_math_func failed to compile'
        write(0,*)  '   aborting...'
      endif
      stop
    endif

  case default
    if ( mpi_node() == 0 ) then
      write(0,*)  ''
      write(0,*)  '   Error in zpulse_flyfoc parameters'
      write(0,*)  '   lon_type must be either "polynomial", "gaussian",'
      write(0,*)  '   "sin2", "math" or "hermite"'
      write(0,*)  '   aborting...'
    endif
    stop
  end select

  this % xi0 = xi0

  this%per_w0 = per_w0

  ! Set default per_center values (this needs to happen after reading the direction )
  if ( p_x_dim > 1 ) then
   if ( per_center(1) == -huge(1.0_p_double) ) then
     if ( direction == 1 ) then
     per_center(1) = 0.5_p_double * ( xmin( g_space, 2 ) + xmax( g_space, 2 ) )
     else
     per_center(1) = 0.5_p_double * ( xmin( g_space, 1 ) + xmax( g_space, 1 ) )
     endif
   endif

   if ( p_x_dim == 3 ) then
    if ( per_center(2) == -huge(1.0_p_double) ) then
      if ( direction /= 3 ) then
      per_center(2) = 0.5_p_double * ( xmin( g_space, 3 ) + xmax( g_space, 3 ) )
      else
      per_center(2) = 0.5_p_double * ( xmin( g_space, 2 ) + xmax( g_space, 2 ) )
      endif
    endif
   endif
  endif

  this%per_center     = per_center

  this%launch_time    = launch_time
  this%if_launch = if_launch

  this%wall_pos = wall_pos
  this%wall_vel = wall_vel
  this%ncells   = ncells
  this%if_launch_from_wall = if_launch_from_wall
  this%if_clean_fields = if_clean_fields
  this%if_apply_divcorr = if_apply_divcorr
  this%if_add_fields = if_add_fields

  ! interp params
  this%if_interp_env = if_interp_env
  this%if_launch_t0     = if_launch_t0
  this%interp_dx2    = interp_dx2 
  this%interp_dt     = interp_dt 

  ! laboratory flying-focus geometry
  this%use_lab_geometry = use_lab_geometry
  this%carrier_angle = real( carrier_angle * pi_180, p_double )
  ! carrier_dir is the preferred general 3D input.  If omitted, retain the
  ! historical carrier_angle convention in the x1-x2 plane.  A partially
  ! specified vector is rejected rather than silently mixing interfaces.
  carrier_dir_any_set = any( abs(carrier_dir) < 0.5_p_double*huge(1.0_p_double) )
  carrier_dir_all_set = all( abs(carrier_dir) < 0.5_p_double*huge(1.0_p_double) )
  if( carrier_dir_any_set .and. .not. carrier_dir_all_set ) then
    SCR_ROOT( 'ERROR in zpulse_flyfoc: carrier_dir requires all three components')
    stop
  endif
  if( carrier_dir_all_set ) then
    this%carrier_dir = carrier_dir
  else
    this%carrier_dir = (/ cos(this%carrier_angle), sin(this%carrier_angle), 0.0_p_double /)
  endif
  this%lab_basis_p = 0.0_p_double
  this%lab_basis_s = 0.0_p_double
  this%focus_start = focus_start
  this%focus_end = focus_end
  this%focus_t_start = focus_t_start
  this%focus_dir = 0.0_p_double
  this%focus_length = 0.0_p_double

  ! shared ff params
  this%z1 = z1
  this%z2 = z2 
  this%initial_delay = initial_delay
  this%final_delay = final_delay
  this%vff = vff 
  this%vg = vg 

  if( this%use_lab_geometry ) then
    if( p_x_dim .lt. 2 .or. p_x_dim .gt. 3 ) then
      SCR_ROOT( 'ERROR in zpulse_flyfoc: use_lab_geometry requires 2D or 3D')
      stop
    endif
    if( .not. this%if_launch_t0 ) then
      SCR_ROOT( 'ERROR in zpulse_flyfoc: use_lab_geometry currently requires if_launch_t0 = .true.')
      stop
    endif
    if( this%if_interp_env ) then
      SCR_ROOT( 'ERROR in zpulse_flyfoc: use_lab_geometry currently requires if_interp_env = .false.')
      stop
    endif
    if( this%propagation .ne. p_forward ) then
      SCR_ROOT( 'ERROR in zpulse_flyfoc: laboratory geometry specifies carrier propagation explicitly; set propagation = "forward"')
      stop
    endif
    if( this%pol_type .ne. 0 ) then
      SCR_ROOT( 'ERROR in zpulse_flyfoc: use_lab_geometry currently supports linear polarization only (pol_type = 0)')
      stop
    endif
    if( abs(this%vff) .le. tiny(1.0_p_double) ) then
      SCR_ROOT( 'ERROR in zpulse_flyfoc: vff must be non-zero in laboratory geometry')
      stop
    endif
    if( this%vg .le. 0.0_p_double ) then
      SCR_ROOT( 'ERROR in zpulse_flyfoc: vg must be > 0 in laboratory geometry')
      stop
    endif

    ! Normalize the carrier vector and build a deterministic orthonormal
    ! transverse basis.  For carriers in the x1-x2 plane this reduces to
    ! p=(-n2,n1,0), s=(0,0,1), exactly matching the validated 2D code.
    carrier_nmag = sqrt( sum( this%carrier_dir**2 ) )
    if( carrier_nmag .le. tiny(1.0_p_double) ) then
      SCR_ROOT( 'ERROR in zpulse_flyfoc: carrier_dir must be non-zero in laboratory geometry')
      stop
    endif
    this%carrier_dir = this%carrier_dir / carrier_nmag

    if( p_x_dim .eq. 2 ) then
      if( abs(this%carrier_dir(3)) > 100.0_p_double*epsilon(1.0_p_double) ) then
        SCR_ROOT( 'ERROR in zpulse_flyfoc: 2D laboratory geometry requires carrier_dir(3)=0')
        stop
      endif
      if( abs(this%focus_start(3)) > 100.0_p_double*epsilon(1.0_p_double) .or. &
          abs(this%focus_end(3)) > 100.0_p_double*epsilon(1.0_p_double) ) then
        SCR_ROOT( 'ERROR in zpulse_flyfoc: 2D laboratory geometry requires focus_start(3)=focus_end(3)=0')
        stop
      endif
    endif

    basis_ref = (/ 0.0_p_double, 0.0_p_double, 1.0_p_double /)
    if( abs(dot_product(basis_ref,this%carrier_dir)) > 0.95_p_double ) then
      basis_ref = (/ 0.0_p_double, 1.0_p_double, 0.0_p_double /)
    endif

    this%lab_basis_p(1) = basis_ref(2)*this%carrier_dir(3) - basis_ref(3)*this%carrier_dir(2)
    this%lab_basis_p(2) = basis_ref(3)*this%carrier_dir(1) - basis_ref(1)*this%carrier_dir(3)
    this%lab_basis_p(3) = basis_ref(1)*this%carrier_dir(2) - basis_ref(2)*this%carrier_dir(1)
    basis_pmag = sqrt( sum(this%lab_basis_p**2) )
    if( basis_pmag .le. tiny(1.0_p_double) ) then
      SCR_ROOT( 'ERROR in zpulse_flyfoc: failed to construct laboratory polarization basis')
      stop
    endif
    this%lab_basis_p = this%lab_basis_p / basis_pmag

    this%lab_basis_s(1) = this%carrier_dir(2)*this%lab_basis_p(3) - this%carrier_dir(3)*this%lab_basis_p(2)
    this%lab_basis_s(2) = this%carrier_dir(3)*this%lab_basis_p(1) - this%carrier_dir(1)*this%lab_basis_p(3)
    this%lab_basis_s(3) = this%carrier_dir(1)*this%lab_basis_p(2) - this%carrier_dir(2)*this%lab_basis_p(1)

    this%focus_dir = this%focus_end - this%focus_start
    this%focus_length = sqrt( sum( this%focus_dir**2 ) )
    if( this%focus_length .le. tiny(1.0_p_double) ) then
      SCR_ROOT( 'ERROR in zpulse_flyfoc: focus_start and focus_end must define a non-zero path')
      stop
    endif
    this%focus_dir = this%focus_dir / this%focus_length
  endif

  this%zr = 0.5_p_double * omega0 * this%per_w0(1)**2 

  if( this%if_interp_env .and. &
        ( ( this%interp_dx2 .le. 0 ) .or. ( this%interp_dt .le. 0 ) ) ) then 
          SCR_ROOT( '   ERROR: when using if_interp_env = .true., interp_dx2 and interp_dt must be set')
          stop
  endif 

  select case( env_type ) 

  case( 'gaussian' )
    this%env_type = p_gaussian_new
    SCR_ROOT( 'using gaussian env' )

  case( 'laguerre' )
    this%env_type = p_laguerre_new

  case( 'astrl' )
    this%env_type = p_astrl

  case( 'astrl_laguerre' )
    this%env_type = p_astrl_laguerre

  case( 'astrl_discrete' )
    this%env_type = p_astrl_discrete

  case default  
    SCR_ROOT( 'ERROR in os-zpulse-flyfoc: invalid pulse type requested.')
    SCR_ROOT( 'Options are: gaussian, astrl, astrl_discrete ')
    stop 

  end select 

  if( this%env_type .eq. p_gaussian_new ) then 

    this%tenv_launch_duration = this%initial_delay + this%final_delay

  else if( this%env_type .eq. p_astrl ) then 

    this%integral_size = integral_size

    if( .not. this%use_lab_geometry ) then
      if( z2 <= z1 ) then
        SCR_ROOT( "ERROR in os-zpulse-flyfoc: z2 > z1 is required when using tunable_env_gaussian")
        stop
      endif
    endif

    if( integral_size .lt. 2 ) then 
      SCR_ROOT( "ERROR in os-zpulse-flyfoc: integral_size >= 2 is required")
      stop
    endif

    if( this%use_lab_geometry ) then
      this%tenv_launch_duration = abs( this%focus_length / this%vff ) + this%tenv_range &
        + this%initial_delay + this%final_delay
    else if( this%propagation .eq. p_forward ) then
      this%tenv_launch_duration = abs( ( this%vg/this%vff - 1.0_p_double ) * (z2 - z1) ) &
        + this%initial_delay + this%final_delay
    else
      this%tenv_launch_duration = abs( ( this%vg/this%vff + 1.0_p_double ) * (z2 - z1) ) &
        + this%initial_delay + this%final_delay
    endif

    select case( focal_weight_type ) 

    case( "const" ) 
      this%focal_weight_type = p_const

    case( "math" )
      this%focal_weight_type = p_func 

      call setup(this%focal_weight_math_func, &
                  trim(focal_weight_math_func), (/'z0'/), ierr)
      
      if (ierr /= 0) then
        if ( mpi_node() == 0 ) then
          write(0,*)  ''
          write(0,*)  '   Error in zpulse_flyfoc parameters'
          write(0,*)  '   Supplied focal_weight_math_func failed to compile'
          write(0,*)  '   aborting...'
        endif
        stop
      endif

    case default 
      SCR_ROOT( '   ERROR: invalid focal_weight_type: ', focal_weight_type )
      SCR_ROOT( '   Options are: const, math ')
      stop 
    
    end select 

    ! to do: set the duration correctly when using custom focal delay 
    select case( focal_delay_type ) 

    case( "const" ) 
      this%focal_delay_type = p_const

    case( "math" )
      this%focal_delay_type = p_func 

      call setup(this%focal_delay_math_func, &
                  trim(focal_delay_math_func), (/'z0'/), ierr)
      
      if (ierr /= 0) then
        if ( mpi_node() == 0 ) then
          write(0,*)  ''
          write(0,*)  '   Error in zpulse_flyfoc parameters'
          write(0,*)  '   Supplied focal_delay_math_func failed to compile'
          write(0,*)  '   aborting...'
        endif
        stop
      endif

    case default 
      SCR_ROOT( '   ERROR: invalid focal_delay_type: ', focal_delay_type )
      SCR_ROOT( '   Options are: const, math ')
      stop 
    
    end select 


    select case( focal_phase_type ) 

    case( "const" ) 
      this%focal_phase_type = p_const

    case( "math" )
      this%focal_phase_type = p_func 

      call setup(this%focal_phase_math_func, &
                  trim(focal_phase_math_func), (/'z0'/), ierr)
      
      if (ierr /= 0) then
        if ( mpi_node() == 0 ) then
          write(0,*)  ''
          write(0,*)  '   Error in zpulse_flyfoc parameters'
          write(0,*)  '   Supplied focal_phase_math_func failed to compile'
          write(0,*)  '   aborting...'
        endif
        stop
      endif

    case default 
      SCR_ROOT( '   ERROR: invalid focal_phase_type: ', focal_phase_type )
      SCR_ROOT( '   Options are: const, math ')
      stop 
    
    end select 


    select case( focal_xperp_type ) 

    case( "const" ) 
      this%focal_xperp_type = p_const

    case( "math" )
      this%focal_xperp_type = p_func 

      call setup(this%focal_xperp_math_func, &
                  trim(focal_xperp_math_func), (/'z0'/), ierr)
      
      if (ierr /= 0) then
        if ( mpi_node() == 0 ) then
          write(0,*)  ''
          write(0,*)  '   Error in zpulse_flyfoc parameters'
          write(0,*)  '   Supplied focal_xperp_math_func failed to compile'
          write(0,*)  '   aborting...'
        endif
        stop
      endif

    case default 
      SCR_ROOT( '   ERROR: invalid focal_xperp_type: ', focal_xperp_type )
      SCR_ROOT( '   Options are: const, math ')
      stop 
    
    end select 


  else if( this%env_type .eq. p_astrl_discrete ) then 

    N = astrl_discrete_num_pulses

    if( N > p_astrl_discrete_num_pulses_max ) then
      SCR_ROOT( '   ERROR: requested astrl_discrete_num_pulses = ', astrl_discrete_num_pulses )
      SCR_ROOT( '   Max num pulses = ', p_astrl_discrete_num_pulses_max )
      stop
    endif

    this%astrl_discrete_num_pulses = N
    this%astrl_discrete_norm_idx = astrl_discrete_norm_idx

    call alloc( this%astrl_discrete_focal_points, (/ N /))
    call alloc( this%astrl_discrete_focal_delays, (/ N /))
    call alloc( this%astrl_discrete_focal_weights, (/ N /))
    call alloc( this%astrl_discrete_focal_phases, (/ N /))
    call alloc( this%astrl_discrete_focal_xperps, (/ N /))

    ! astrl discrete -- focal points 
    if( astrl_discrete_use_linspace ) then 
      do i=1, N 
        this%astrl_discrete_focal_points(i) = this%z1 + (this%z2 - this%z1 ) / (N-1)* (i-1)  
      enddo
    else
      this%astrl_discrete_focal_points(1:N) = astrl_discrete_focal_points(1:N)
    endif

    ! astrl discrete --  focal delay 
    if( focal_delay_type .eq. 'math' ) then  

      call setup(this%focal_delay_math_func, &
                  trim(focal_delay_math_func), (/'z0'/), ierr)
      
      if (ierr /= 0) then
        if ( mpi_node() == 0 ) then
          write(0,*)  ''
          write(0,*)  '   Error in zpulse_flyfoc parameters'
          write(0,*)  '   Supplied focal_delay_math_func failed to compile'
          write(0,*)  '   aborting...'
        endif
        stop
      endif

      do i=1, N 
        fparser_arr(1) = this%astrl_discrete_focal_points(i)
        this%astrl_discrete_focal_delays(i) = eval( this%focal_delay_math_func, fparser_arr )
      enddo

    else
      this%astrl_discrete_focal_delays(1:N) = astrl_discrete_focal_delays(1:N) 
    endif

    ! astrl discrete -- focal weight 
    if( focal_weight_type .eq. 'math' ) then  

      call setup(this%focal_weight_math_func, &
                  trim(focal_weight_math_func), (/'z0'/), ierr)
      
      if (ierr /= 0) then
        if ( mpi_node() == 0 ) then
          write(0,*)  ''
          write(0,*)  '   Error in zpulse_flyfoc parameters'
          write(0,*)  '   Supplied focal_weight_math_func failed to compile'
          write(0,*)  '   aborting...'
        endif
        stop
      endif

      do i=1, N 
        fparser_arr(1) = this%astrl_discrete_focal_points(i)
        this%astrl_discrete_focal_weights(i) = eval( this%focal_weight_math_func, fparser_arr )
      enddo

    else
      this%astrl_discrete_focal_weights(1:N) = astrl_discrete_focal_weights(1:N) 
    endif

    ! astrl discrete --  focal phase 
    if( focal_phase_type .eq. 'math' ) then  

      call setup(this%focal_phase_math_func, &
                  trim(focal_phase_math_func), (/'z0'/), ierr)
      
      if (ierr /= 0) then
        if ( mpi_node() == 0 ) then
          write(0,*)  ''
          write(0,*)  '   Error in zpulse_flyfoc parameters'
          write(0,*)  '   Supplied focal_phase_math_func failed to compile'
          write(0,*)  '   aborting...'
        endif
        stop
      endif

      do i=1, N 
        fparser_arr(1) = this%astrl_discrete_focal_points(i)
        this%astrl_discrete_focal_phases(i) = eval( this%focal_phase_math_func, fparser_arr )
      enddo

    else
      this%astrl_discrete_focal_phases(1:N) = astrl_discrete_focal_phases(1:N) 
    endif

    ! astrl discrete --  focal xperp 
    if( focal_xperp_type .eq. 'math' ) then  

      call setup(this%focal_xperp_math_func, &
                  trim(focal_xperp_math_func), (/'z0'/), ierr)
      
      if (ierr /= 0) then
        if ( mpi_node() == 0 ) then
          write(0,*)  ''
          write(0,*)  '   Error in zpulse_flyfoc parameters'
          write(0,*)  '   Supplied focal_xperp_math_func failed to compile'
          write(0,*)  '   aborting...'
        endif
        stop
      endif

      do i=1, N 
        fparser_arr(1) = this%astrl_discrete_focal_points(i)
        this%astrl_discrete_focal_xperps(i) = eval( this%focal_xperp_math_func, fparser_arr )
      enddo

    else
      this%astrl_discrete_focal_xperps(1:N) = astrl_discrete_focal_xperps(1:N) 
    endif


    this%tenv_launch_duration = (maxval(this%astrl_discrete_focal_delays) &
                              - minval(this%astrl_discrete_focal_delays)) &
                                + this%initial_delay + this%final_delay 

  endif 


  ! Boost pulse
  this%if_boost = .false.
  
  if (sim_options%gamma > 1.0_p_double) then

    this % gamma          = sim_options%gamma
    this % gamma_times_beta = sqrt(this%gamma**2-1.0_p_double)

    if( this % propagation .eq. p_forward ) then
      this % boost_per_fac = 1.0_p_double / ( this%gamma + this%gamma_times_beta )
    else
      this % boost_per_fac = this%gamma + this%gamma_times_beta
    endif 
    SCR_ROOT('   boosting zpulse with gamma = ', this%gamma )

    this%if_boost = .true.
  else
    this % boost_per_fac = 1.0_p_double
  endif

  if( this%use_lab_geometry ) then
    if( this%env_type .ne. p_astrl ) then
      SCR_ROOT( 'ERROR in zpulse_flyfoc: laboratory geometry currently supports env_type = "astrl" only')
      stop
    endif
    if( this%focal_delay_type .ne. p_const ) then
      SCR_ROOT( 'ERROR in zpulse_flyfoc: laboratory geometry derives focal delay from focus path and vff; use focal_delay_type = "const"')
      stop
    endif
    if( this%focal_xperp_type .ne. p_const ) then
      SCR_ROOT( 'ERROR in zpulse_flyfoc: laboratory geometry currently requires focal_xperp_type = "const"')
      stop
    endif
  endif

  call this%init_extra_vars() 

  call this%init_env_norm( grid, g_space ) 

end subroutine read_input_zpulse_flyfoc
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine init_extra_vars_zpulse_flyfoc( this )

  class( t_zpulse_flyfoc ), intent(inout) :: this

  if( this%focal_delay_type .eq. p_const .and. .not. this%use_lab_geometry ) then
    ! Choose the reference focal plane so that the programmed delay is
    ! non-negative across [z1,z2].  For a carrier propagating with sign s,
    ! the moving-focus delay slope is vg/vff - s.
    if( this%propagation .eq. p_forward ) then
      if( (this%vg/this%vff - 1.0_p_double) >= 0.0_p_double ) then
        this%z_lead = this%z1
      else
        this%z_lead = this%z2
      endif
    else
      if( (this%vg/this%vff + 1.0_p_double) >= 0.0_p_double ) then
        this%z_lead = this%z1
      else
        this%z_lead = this%z2
      endif
    endif
  endif

  select case( this%env_type )

  case( p_astrl ) 

    call alloc( this%buf_1d, (/ this%integral_size /) )
    this%buf_1d(:) = 0

    if( this%use_lab_geometry ) then
      this%dz0 = this%focus_length / real( this%integral_size - 1, p_double )
    else
      this%dz0 = ( this%z2 - this%z1 ) / real( this%integral_size - 1, p_double )
    endif

  end select 

end subroutine init_extra_vars_zpulse_flyfoc
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine init_env_norm_zpulse_flyfoc( this, grid, g_space )

  class( t_zpulse_flyfoc ), intent(inout) :: this
  class( t_grid ), intent(in) :: grid
  type( t_space ),     intent(in) :: g_space

  real(p_double) :: default_env_norm 
  complex(p_double) :: tmp_env_norm 
  real(p_double) :: t, z0, x2, x3, delay, prop_sign, s0
  real(p_double), dimension(3) :: r0
  real(p_k_fparse), dimension(1) :: fparser_arr 
  real(p_double) :: gxmin2, gxmax2, dx2
  integer :: nx2 

  default_env_norm = this%boost_per_fac * this%omega0 * this%a0

  x2 = 0.0_p_double

  gxmin2 = g_space%x_bnd_initial(p_lower,2)
  gxmax2 = g_space%x_bnd_initial(p_upper,2)
  nx2 = grid%g_nx(2)
  dx2 = ( gxmax2 - gxmin2 ) / nx2
  x2 = this%per_center(1)
  x3 = this%per_center(2) 

  select case( this%env_type )

  case( p_astrl ) 
    
    ! set this so that the envelope function is unnormalized when we 
    ! get the central value 
    this%env_norm = 1.0_p_double

    if( this%use_lab_geometry ) then
      ! Normalize at the midpoint of the programmed laboratory focus path,
      ! at the instant when the moving focus reaches that point.
      s0 = 0.5_p_double * this%focus_length
      r0 = this%focus_start + s0 * this%focus_dir
      t = this%focus_t_start + s0 / this%vff
      select case( p_x_dim )
      case(2)
        call this%get_env_lab_2d( t, r0(1), r0(2), tmp_env_norm )
      case(3)
        call this%get_env_lab_3d( t, r0(1), r0(2), r0(3), tmp_env_norm )
      end select
      if( abs(tmp_env_norm) .le. tiny(1.0_p_double) ) then
        SCR_ROOT( 'ERROR in zpulse_flyfoc: zero envelope encountered while normalizing laboratory geometry')
        stop
      endif
      this%env_norm = default_env_norm / abs( tmp_env_norm )

    else

    ! normalize envelope to the central value of z 
    z0 = ( this%z1 + this%z2 ) / 2 

    ! compute the time at which we will hit this value of z 
    select case( this%focal_delay_type )

    case( p_const )

      if( this%propagation .eq. p_forward ) then
        prop_sign = 1.0_p_double
      else
        prop_sign = -1.0_p_double
      endif
      delay = ( z0 - this%z_lead ) * ( this%vg/this%vff - prop_sign )

    case( p_func )

      fparser_arr(1) = z0
      delay = eval( this%focal_delay_math_func, fparser_arr ) 

    end select

    ! Normalize at the time when the z0-focused constituent is at the
    ! centre of its temporal envelope.  t_envelope uses t-z+xi0 for a
    ! forward wave and t+z-xi0 for a backward wave.
    if( this%propagation .eq. p_forward ) then
      prop_sign = 1.0_p_double
    else
      prop_sign = -1.0_p_double
    endif
    t = prop_sign*z0 + delay - prop_sign*this%xi0 + this%lon_center()

    select case( p_x_dim )

    case(2)
      call this%get_env_astrl_2d( t, z0, x2, tmp_env_norm )
      
    case(3)
      call this%get_env_astrl_3d( t, z0, x2, x3, tmp_env_norm )
    
    end select 

    this%env_norm = default_env_norm / abs( tmp_env_norm )

    endif


  case( p_astrl_discrete )

    this%env_norm = 1.0_p_double

    z0 = this%astrl_discrete_focal_points( this%astrl_discrete_norm_idx )
    delay = this%astrl_discrete_focal_delays( this%astrl_discrete_norm_idx )
    
    if( this%propagation .eq. p_forward ) then
      prop_sign = 1.0_p_double
    else
      prop_sign = -1.0_p_double
    endif
    t = prop_sign*z0 + delay - prop_sign*this%xi0 + this%lon_center()

    select case( p_x_dim ) 

    case(2) 
      call this%get_env_astrl_discrete_2d( t, z0, x2, tmp_env_norm )

    case(3)
      call this%get_env_astrl_discrete_3d( t, z0, x2, x3, tmp_env_norm )

    end select

    this%env_norm = default_env_norm / abs( tmp_env_norm )


  ! everything else 
  case default 
    
    this%env_norm = default_env_norm

  end select 

  ! circular polarization: reduce field component amplitudes so that the complex 
  ! amplitude is the requested a0 
  if (this%pol_type /= 0) then
    this%env_norm = this%env_norm / sqrt(2.0_p_double)
  endif

end subroutine init_env_norm_zpulse_flyfoc
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
! Total duration of the laser pulse
!-----------------------------------------------------------------------------------------
function t_duration_zpulse_flyfoc( this ) result( t_duration )

  implicit none

  class( t_zpulse_flyfoc ), intent(in) :: this
  real(p_double) :: t_duration

  t_duration = this%tenv_launch_duration

  if( this % if_boost ) then
    t_duration = this % gamma * t_duration!  - this % gamma_times_beta * this % wall_pos
  endif

end function t_duration_zpulse_flyfoc
!----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
function t_env_zpulse_flyfoc( this, t_sim, z_sim ) result(t_env)

  implicit none

  class( t_zpulse_flyfoc ), intent(inout) :: this ! intent must be inout because of
                                                ! the eval( f_parser ) function
  real(p_double), intent(in) :: t_sim
  real(p_double), intent(in) :: z_sim 

  ! lab coordinates 
  real(p_double) :: z, t

  real(p_double) :: t_env

  ! local variables

  real(p_k_fparse), dimension(1) ::t_dbl

  t = t_sim
  z = z_sim 

  ! ignore pulse broadening 
  if( this%propagation .eq. p_forward ) then 
    t = t - z + this % xi0
  else
    t = t + z - this % xi0
  endif

  select case ( this%tenv_type )

  case( p_const )
    t_env = 1.0_p_double 

  case (p_polynomial)
    t_env = fenv_poly( t, this%tenv_rise, this%tenv_flat, this%tenv_fall )

  case (p_sin2)
    t_env = fenv_sin2( t, this%tenv_rise, this%tenv_flat, this%tenv_fall )

  case (p_gaussian)
    t_env = exp( - 2*((t - 0.5*this % tenv_range)/this%tenv_duration)**2 )


  case (p_func)

    if ( t > this % tenv_range ) then
      t_env = 0
    else
      t_dbl(1) = t
      t_env = eval( this%lon_math_func, t_dbl )
    endif

  case default
    t_env = 0

  end select

end function t_env_zpulse_flyfoc
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine cleanup_zpulse_flyfoc( this )

  class( t_zpulse_flyfoc ), intent(inout) :: this

  if ( associated(this%e_env_buffer_2d) ) call freemem( this%e_env_buffer_2d )
  if ( associated(this%e_env_buffer_3d) ) call freemem( this%e_env_buffer_3d )

  if ( associated(this%buf_1d)) call freemem( this%buf_1d )

  if( associated( this%astrl_discrete_focal_delays )) &
      call freemem( this%astrl_discrete_focal_delays)
  if( associated( this%astrl_discrete_focal_points )) &
      call freemem( this%astrl_discrete_focal_points)
  if( associated( this%astrl_discrete_focal_weights )) &
      call freemem( this%astrl_discrete_focal_weights)
  if( associated( this%astrl_discrete_focal_phases )) &
      call freemem( this%astrl_discrete_focal_phases)
  if( associated( this%astrl_discrete_focal_xperps )) &
      call freemem( this%astrl_discrete_focal_xperps )

end subroutine cleanup_zpulse_flyfoc
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
function gaussian_beam_2d( z, x, omega, zr, z0, xcenter ) result( ret )

  real(p_double), intent(in) :: z, x, omega, zr, z0, xcenter
  complex(p_double) :: ret
  complex(p_double) :: q 

  q = cmplx( z - z0, zr )
  ret = exp( cmplx(0,-omega)  * ( x - xcenter ) ** 2 / (2 * q) ) / sqrt( q / zr ) 

end function gaussian_beam_2d
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
function gaussian_beam_3d( z, x, y, omega, zr, z0, xcenter, ycenter ) result( ret )

  real(p_double), intent(in) :: z, x, y, omega, zr, z0, xcenter, ycenter
  complex(p_double) :: ret
  complex(p_double) :: q 

  q = cmplx( z - z0, zr )
  ret = exp( cmplx(0,-omega)  * ( (y-ycenter)**2 + (x-xcenter)**2) / (2 * q) ) / ( q / zr ) 

end function gaussian_beam_3d
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
function laguerre_beam_3d( z, x, y, omega, zr, w0, z0, xcenter, ycenter, ell, p ) result( ret )

  implicit none 

  real(p_double), intent(in) :: z, x, y, omega, zr, w0, z0, xcenter, ycenter
  integer, intent(in) ::  ell, p
  real(p_double) :: r, theta, zcent, w_of_z, Rinverse_of_z, guoy_shift
  complex(p_double) :: ret

  r = sqrt( (x - xcenter)**2 + (y - ycenter )**2 )
  theta = atan2( x - xcenter, y - ycenter )

  zcent  = z - z0 
  w_of_z = w0 * sqrt( 1 + (zcent / zr )**2 )
  Rinverse_of_z = zcent / ( zcent**2 + zr**2 )
  guoy_shift = ( 2*p + abs(ell) + 1 ) * atan2( zcent, z0 )

  ret = (w0 / w_of_z ) &
            * ( r * sqrt(2.0_p_double) / w_of_z ) ** abs(ell) & 
            * exp( - ( r / w_of_z )**2 ) &
            * laguerre( p, ell, 2 * ( r / w_of_z ) **2 ) & 
            * exp( cmplx( 0, 1) * ( - 0.5_p_double * omega * r**2 * Rinverse_of_z &
                                    + guoy_shift + ell * theta ) ) 

end function laguerre_beam_3d
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
! multilinear 3D array interpolation 
! not optimized 
!-----------------------------------------------------------------------------------------
function interp_3d( x, ix1, x1, dx, farr ) result( f )

  real(p_double), dimension(3), intent(in) :: x,  x1, dx 
  integer, dimension(3), intent(in) :: ix1
  complex(p_double), dimension(:,:,:), intent(in) :: farr 
  complex(p_double) :: f

  real(p_double), dimension(3,2) :: xdiff
  integer :: i, j, k 

  f = 0 

  xdiff(:,2) = x - x1
  xdiff(:,1) = x1 + dx - x  

  do i=1,2
    do j=1,2
      do k=1,2
        f = f + farr( ix1(1) + i, ix1(2) + j, ix1(3) + k ) * xdiff(1,i) * xdiff(2,j) * xdiff(3,k)
      enddo
    enddo
  enddo

  f = f / ( dx(1) * dx(2) * dx(3) )

end function interp_3d
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
! 4d multilinear interpolation 
! not optimized 
!-----------------------------------------------------------------------------------------
function interp_4d( x, ix1, x1, dx, farr ) result( f )

  real(p_double), dimension(4), intent(in) :: x,  x1, dx 
  integer, dimension(4), intent(in) :: ix1
  complex(p_double), dimension(:,:,:,:), intent(in) :: farr 
  complex(p_double) :: f

  real(p_double), dimension(4,2) :: xdiff
  integer :: i1, i2, i3, i4

  f = 0 

  xdiff(:,2) = x - x1
  xdiff(:,1) = x1 + dx - x  

  do i1=1,2
    do i2=1,2
      do i3=1,2
        do i4 = 1,2
          f = f + farr( ix1(1) + i1, ix1(2) + i2, ix1(3) + i3, ix1(4) + i4 ) &
                * xdiff(1,i1) * xdiff(2,i2) * xdiff(3,i3) * xdiff(4,i4)
        enddo
      enddo
    enddo
  enddo

  f = f / ( dx(1) * dx(2) * dx(3) * dx(4) )

end function interp_4d
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine precompute_envelope_2d_zpulse_flyfoc( this, e, g_space, nx_p_min, dt ) 

  integer, parameter :: rank = 2

  class( t_zpulse_flyfoc ), intent(inout) :: this
  type( t_vdf ), intent(in) :: e
  type( t_space ), intent(in) :: g_space
  integer, dimension(:), intent(in) :: nx_p_min
  real(p_double), intent(in) :: dt

  real(p_double), dimension(2,rank) :: g_x_range
  real(p_double), dimension(rank) :: ldx

  integer :: nt, nz, nx2, nx2_fft
  integer :: i0, i1, i2

  ! bounds of the z region that must be computed 
  real(p_double) :: tmin, tmax, zmin, zmax, x2min, x2max, delta_t, zwall  
  real(p_double) :: x2, z, t
  real(p_double) :: interp_dx1

  complex(p_double) :: env

  ! progress messages
  integer :: progress_counter, num_progress_intervals

  ! type(t_fft_manager) :: fft_manager 

  SCR_ROOT( "Computing pulse envelope..." )

  ! calculate spatial data 
  call get_x_bnd( g_space, g_x_range )

  ldx(1:2) = real( e%dx_(1:2), p_double )

  if( .not. this%if_launch_from_wall ) then 

    zmin = this%wall_pos - 2 * ldx(1) 
    zmax = this%wall_pos + (this%ncells + 2) * ldx(1)

  else

    if( this%propagation .eq. p_forward ) then 

      zwall = real( g_x_range(p_lower,1) )
      zmin = zwall - 2 * ldx(1)
      zmax = zwall + ( this%ncells + 2 ) * ldx(1) 

    else

      zwall = real( g_x_range(p_upper,1) )
      zmin = zwall - ( this%ncells + 2 ) * ldx(1)
      zmax = zwall + 2 * ldx(1) 

    endif

  endif

  x2min = real( g_x_range(p_lower,2) ) + ( lbound( e%f2, 3 ) + nx_p_min(2) - 3 ) * ldx(2) 
  x2max = real( g_x_range(p_lower,2) ) + ( ubound( e%f2, 3 ) + nx_p_min(2) + 3 ) * ldx(2) 

  select case( this%interpolation )
  case( p_linear, p_cubic )
    ! Nothing to change here
  case( p_quadratic, p_quartic )
    zmin = zmin + ldx(1)*0.5_p_double
    zmax = zmax + ldx(1)*0.5_p_double
    x2min = x2min + ldx(2)*0.5_p_double
    x2max = x2max + ldx(2)*0.5_p_double
  end select

  if( .not. this%if_boost ) then 
    tmin = 0 - this%interp_dt  
    tmax = this%tenv_launch_duration + this%interp_dt 
    interp_dx1 = ldx(1) 
  else 
    tmin = this%gamma_times_beta * this%wall_pos
    tmax = this%tenv_launch_duration + ldx(1) * this%gamma_times_beta
    zmin = this%gamma * zmin
    zmax = this%gamma * zmax 
    interp_dx1 = this%interp_dt
  endif

  ! temporarily disable envelope field interpolation so that we can get 
  ! the exact field values within this function call -- unset at end of this routine
  this%if_interp_env = .false.


  delta_t = this%interp_dt

  ! Number of interpolation intervals in each direction.  Multilinear
  ! interpolation needs both endpoints of every interval, so the buffer
  ! contains ninterval+1 grid points in each direction.
  nt = max( 1, ceiling( (tmax - tmin )/ this%interp_dt ) )
  nz = max( 1, ceiling( (zmax - zmin) / interp_dx1 ) )
  nx2 = max( 1, ceiling( (x2max - x2min) / this%interp_dx2 ) )

  this%e_env_grid_shape(1:3) = (/ nt + 1, nz + 1, nx2 + 1 /)

  this%e_env_grid_bounds(1,1:2) = (/ tmin, tmax /)
  this%e_env_grid_bounds(2,1:2) = (/ zmin, zmax /)
  this%e_env_grid_bounds(3,1:2) = (/ x2min, x2max /)

  this%e_env_grid_dx(1) = ( tmax - tmin ) / real(nt, p_double)
  this%e_env_grid_dx(2) = ( zmax - zmin ) / real(nz, p_double)
  this%e_env_grid_dx(3) = ( x2max - x2min ) / real(nx2, p_double)

  call alloc( this%e_env_buffer_2d, (/ nt + 1, nz + 1, nx2 + 1 /) )
  this%e_env_buffer_2d(:,:,:) = 0 

  progress_counter = 0
  num_progress_intervals = 10

  do i0 = 1, nt + 1

    if( i0 * 1.0_p_double / (nt + 1) >= ( progress_counter * 1.0_p_double / num_progress_intervals ) ) then
      SCR_ROOT( i0, " / ", nt + 1, "..." )
      progress_counter = progress_counter + 1
    endif

    t = tmin + ( i0 - 1 ) * this%e_env_grid_dx(1) 

    do i1 = 1, nz + 1
 
      z = zmin + ( i1 - 1 ) * this%e_env_grid_dx(2) 

      do i2 = 1, nx2 + 1
        x2 = x2min + (i2 - 1 ) * this%e_env_grid_dx(3)

        call this%get_env_2d( t, z, x2, env )

        this%e_env_buffer_2d(i0, i1, i2 ) =  env 

      enddo 

    enddo

  enddo

  ! print *, mpi_node(), 'max env: ' , maxval(abs(this%e_env_buffer_2d)) 

  SCR_ROOT( "Successfully computed envelope." )

  ! revert if_interp_env -- now, the fields can be interpolated. 
  this%if_interp_env = .true.

end subroutine precompute_envelope_2d_zpulse_flyfoc
!-----------------------------------------------------------------------------------------




!-----------------------------------------------------------------------------------------
subroutine precompute_envelope_3d_zpulse_flyfoc( this, e, g_space, nx_p_min, dt )

  integer, parameter :: rank = 3

  class( t_zpulse_flyfoc ), intent(inout) :: this
  type( t_vdf ), intent(in) :: e
  type( t_space ), intent(in) :: g_space
  integer, dimension(:), intent(in) :: nx_p_min
  real(p_double), intent(in) :: dt

  real(p_double), dimension(2,rank) :: g_x_range
  real(p_double), dimension(rank) :: ldx

  integer :: nt, nz, nx2, nx3
  integer :: i0, i1, i2, i3

  real(p_double) :: tmin, tmax, zmin, zmax, x2min, x2max, x3min, x3max, zwall
  real(p_double) :: x2, x3, z, t
  real(p_double) :: interp_dx1

  complex(p_double) :: env

  integer :: progress_counter, num_progress_intervals

  SCR_ROOT( "Computing 3D pulse envelope..." )

  call get_x_bnd( g_space, g_x_range )
  ldx(1:3) = real( e%dx_(1:3), p_double )

  if( .not. this%if_launch_from_wall ) then
    zmin = this%wall_pos - 2 * ldx(1)
    zmax = this%wall_pos + (this%ncells + 2) * ldx(1)
  else
    if( this%propagation .eq. p_forward ) then
      zwall = real( g_x_range(p_lower,1) )
      zmin = zwall - 2 * ldx(1)
      zmax = zwall + ( this%ncells + 2 ) * ldx(1)
    else
      zwall = real( g_x_range(p_upper,1) )
      zmin = zwall - ( this%ncells + 2 ) * ldx(1)
      zmax = zwall + 2 * ldx(1)
    endif
  endif

  x2min = real( g_x_range(p_lower,2) ) + ( lbound( e%f3, 3 ) + nx_p_min(2) - 3 ) * ldx(2)
  x2max = real( g_x_range(p_lower,2) ) + ( ubound( e%f3, 3 ) + nx_p_min(2) + 3 ) * ldx(2)
  x3min = real( g_x_range(p_lower,3) ) + ( lbound( e%f3, 4 ) + nx_p_min(3) - 3 ) * ldx(3)
  x3max = real( g_x_range(p_lower,3) ) + ( ubound( e%f3, 4 ) + nx_p_min(3) + 3 ) * ldx(3)

  select case( this%interpolation )
  case( p_linear, p_cubic )
    ! no coordinate correction
  case( p_quadratic, p_quartic )
    zmin = zmin + ldx(1)*0.5_p_double
    zmax = zmax + ldx(1)*0.5_p_double
    x2min = x2min + ldx(2)*0.5_p_double
    x2max = x2max + ldx(2)*0.5_p_double
    x3min = x3min + ldx(3)*0.5_p_double
    x3max = x3max + ldx(3)*0.5_p_double
  end select

  if( .not. this%if_boost ) then
    tmin = 0.0_p_double - this%interp_dt
    tmax = this%tenv_launch_duration + this%interp_dt
    interp_dx1 = ldx(1)
  else
    tmin = this%gamma_times_beta * this%wall_pos
    tmax = this%tenv_launch_duration + ldx(1) * this%gamma_times_beta
    zmin = this%gamma * zmin
    zmax = this%gamma * zmax
    interp_dx1 = this%interp_dt
  endif

  ! Evaluate exact envelopes while constructing the interpolation table.
  this%if_interp_env = .false.

  nt  = max( 1, ceiling( (tmax - tmin ) / this%interp_dt ) )
  nz  = max( 1, ceiling( (zmax - zmin ) / interp_dx1 ) )
  nx2 = max( 1, ceiling( (x2max - x2min) / this%interp_dx2 ) )
  nx3 = max( 1, ceiling( (x3max - x3min) / this%interp_dx2 ) )

  this%e_env_grid_shape(1:4) = (/ nt + 1, nz + 1, nx2 + 1, nx3 + 1 /)

  this%e_env_grid_bounds(1,1:2) = (/ tmin, tmax /)
  this%e_env_grid_bounds(2,1:2) = (/ zmin, zmax /)
  this%e_env_grid_bounds(3,1:2) = (/ x2min, x2max /)
  this%e_env_grid_bounds(4,1:2) = (/ x3min, x3max /)

  this%e_env_grid_dx(1) = ( tmax - tmin ) / real(nt, p_double)
  this%e_env_grid_dx(2) = ( zmax - zmin ) / real(nz, p_double)
  this%e_env_grid_dx(3) = ( x2max - x2min ) / real(nx2, p_double)
  this%e_env_grid_dx(4) = ( x3max - x3min ) / real(nx3, p_double)

  call alloc( this%e_env_buffer_3d, (/ nt + 1, nz + 1, nx2 + 1, nx3 + 1 /) )
  this%e_env_buffer_3d(:,:,:,:) = 0

  progress_counter = 0
  num_progress_intervals = 10

  do i0 = 1, nt + 1
    if( i0 * 1.0_p_double / (nt + 1) >= &
        ( progress_counter * 1.0_p_double / num_progress_intervals ) ) then
      SCR_ROOT( i0, " / ", nt + 1, "..." )
      progress_counter = progress_counter + 1
    endif

    t = tmin + ( i0 - 1 ) * this%e_env_grid_dx(1)

    do i1 = 1, nz + 1
      z = zmin + ( i1 - 1 ) * this%e_env_grid_dx(2)

      do i2 = 1, nx2 + 1
        x2 = x2min + ( i2 - 1 ) * this%e_env_grid_dx(3)

        do i3 = 1, nx3 + 1
          x3 = x3min + ( i3 - 1 ) * this%e_env_grid_dx(4)

          call this%get_env_3d( t, z, x2, x3, env )
          this%e_env_buffer_3d(i0, i1, i2, i3) = env
        enddo
      enddo
    enddo
  enddo

  SCR_ROOT( "Successfully computed 3D envelope." )
  this%if_interp_env = .true.

end subroutine precompute_envelope_3d_zpulse_flyfoc
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine get_env_2d_zpulse_flyfoc( this, t, z, x2, env )

  class( t_zpulse_flyfoc ), intent(inout) :: this
  real(p_double), intent(in) :: t, z, x2 
  complex(p_double), intent(out) :: env 

  real(p_double), dimension(3) :: x, x1, dx
  integer, dimension(3) :: ix1

  if( this%use_lab_geometry ) then
    call this%get_env_lab_2d( t, z, x2, env )
    return
  endif

  if( this%if_interp_env ) then

    dx(1:3) = this%e_env_grid_dx(1:3)

    ix1(1) = floor( (t - this%e_env_grid_bounds(1,1)) / dx(1) )
    ix1(1) = max( 0, min( ix1(1), this%e_env_grid_shape(1) - 2 ) )
    x1(1) = dx(1) * ix1(1) + this%e_env_grid_bounds(1,1)
    x(1) = t

    ix1(2) = floor( (z - this%e_env_grid_bounds(2,1)) / dx(2) )
    ix1(2) = max( 0, min( ix1(2), this%e_env_grid_shape(2) - 2 ) )
    x1(2) = dx(2) * ix1(2) + this%e_env_grid_bounds(2,1)
    x(2) = z

    ix1(3) = floor( (x2 - this%e_env_grid_bounds(3,1)) / dx(3) )
    ix1(3) = max( 0, min( ix1(3), this%e_env_grid_shape(3) - 2 ) )
    x1(3) = dx(3) * ix1(3) + this%e_env_grid_bounds(3,1)
    x(3) = x2

    env = interp_3d( x, ix1, x1, dx, this%e_env_buffer_2d )

  else

    select case( this%env_type )
  
    case( p_gaussian_new ) 
      call this%get_env_gaussian_2d( t, z, x2, env )

    case( p_astrl )
      call this%get_env_astrl_2d( t, z, x2, env ) 

    case( p_astrl_discrete ) 
      call this%get_env_astrl_discrete_2d( t, z, x2, env ) 

    end select

  endif

end subroutine get_env_2d_zpulse_flyfoc
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine get_env_3d_zpulse_flyfoc( this, t, z, x2, x3, env )

  class( t_zpulse_flyfoc ), intent(inout) :: this
  real(p_double), intent(in) :: t, z, x2, x3
  complex(p_double), intent(out) :: env 

  real(p_double), dimension(4) :: x, x1, dx
  integer, dimension(4) :: ix1

  if( this%use_lab_geometry ) then
    call this%get_env_lab_3d( t, z, x2, x3, env )
    return
  endif

  if( this%if_interp_env ) then

    dx(1:4) = this%e_env_grid_dx(1:4)

    ix1(1) = floor( (t - this%e_env_grid_bounds(1,1)) / dx(1) )
    ix1(1) = max( 0, min( ix1(1), this%e_env_grid_shape(1) - 2 ) )
    x1(1) = dx(1) * ix1(1) + this%e_env_grid_bounds(1,1)
    x(1) = t

    ix1(2) = floor( (z - this%e_env_grid_bounds(2,1)) / dx(2) )
    ix1(2) = max( 0, min( ix1(2), this%e_env_grid_shape(2) - 2 ) )
    x1(2) = dx(2) * ix1(2) + this%e_env_grid_bounds(2,1)
    x(2) = z

    ix1(3) = floor( (x2 - this%e_env_grid_bounds(3,1)) / dx(3) )
    ix1(3) = max( 0, min( ix1(3), this%e_env_grid_shape(3) - 2 ) )
    x1(3) = dx(3) * ix1(3) + this%e_env_grid_bounds(3,1)
    x(3) = x2

    ix1(4) = floor( (x3 - this%e_env_grid_bounds(4,1)) / dx(4) )
    ix1(4) = max( 0, min( ix1(4), this%e_env_grid_shape(4) - 2 ) )
    x1(4) = dx(4) * ix1(4) + this%e_env_grid_bounds(4,1)
    x(4) = x3

    env = interp_4d( x, ix1, x1, dx, this%e_env_buffer_3d )

  else

    select case( this%env_type )
  
    case( p_gaussian_new ) 
      call this%get_env_gaussian_3d( t, z, x2, x3, env )

    case( p_laguerre_new ) 
      call this%get_env_laguerre_3d( t, z, x2, x3, env )

    case( p_astrl )
      call this%get_env_astrl_3d( t, z, x2, x3, env ) 

    case( p_astrl_laguerre )
      call this%get_env_astrl_laguerre_3d( t, z, x2, x3, env ) 

    case( p_astrl_discrete ) 
      call this%get_env_astrl_discrete_3d( t, z, x2, x3, env ) 

    end select

  endif

end subroutine get_env_3d_zpulse_flyfoc
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine get_env_gaussian_2d_zpulse_flyfoc( this, t, z, x2, env )

  class( t_zpulse_flyfoc ), intent(inout) :: this
  real(p_double), intent(in) :: t, z, x2 
  complex(p_double), intent(out) :: env 

  env = gaussian_beam_2d( z, x2, this%omega0, this%zr, this%z1, this%per_center(1) ) * this%t_envelope( t, z ) 

  env = env * this%env_norm 

end subroutine get_env_gaussian_2d_zpulse_flyfoc
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine get_env_gaussian_3d_zpulse_flyfoc( this, t, z, x2, x3, env )

  class( t_zpulse_flyfoc ), intent(inout) :: this
  real(p_double), intent(in) :: t, z, x2, x3
  complex(p_double), intent(out) :: env 

  env = gaussian_beam_3d( z, x2, x3, this%omega0, this%zr, this%z1, &
                          this%per_center(1), this%per_center(2) ) &
                        * this%t_envelope( t, z ) 

  env = env * this%env_norm 

end subroutine get_env_gaussian_3d_zpulse_flyfoc
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine get_env_laguerre_3d_zpulse_flyfoc( this, t, z, x2, x3, env )

  class( t_zpulse_flyfoc ), intent(inout) :: this
  real(p_double), intent(in) :: t, z, x2, x3
  complex(p_double), intent(out) :: env 

  env = laguerre_beam_3d( z, x2, x3, this%omega0, this%zr, this%per_w0(1), this%z1, &
                          this%per_center(1), this%per_center(2), &
                          this%per_tem_mode(2), this%per_tem_mode(1) ) &
                        * this%t_envelope( t, z ) 

  env = env * this%env_norm 

end subroutine get_env_laguerre_3d_zpulse_flyfoc
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine get_env_lab_2d_zpulse_flyfoc( this, t, x1, x2, env )

  ! Evaluate the continuous ASTrL flying-focus envelope in laboratory
  ! coordinates.  The carrier direction n and the focal trajectory are
  ! independent.  Focal points lie along focus_start -> focus_end.

  class( t_zpulse_flyfoc ), intent(inout) :: this
  real(p_double), intent(in) :: t, x1, x2
  complex(p_double), intent(out) :: env

  integer :: i3, nz0
  real(p_double) :: s0, delay, t_focus, beam_coord
  real(p_double) :: n1, n2, m1, m2, u, v
  real(p_double), dimension(3) :: r0, dr
  real(p_k_fparse), dimension(1) :: fparser_arr

  nz0 = this%integral_size
  n1 = this%carrier_dir(1)
  n2 = this%carrier_dir(2)
  m1 = -n2
  m2 =  n1

  ! Absolute longitudinal coordinate along the carrier direction.  The
  ! temporal envelope propagates at vg; the optical carrier remains at c.
  beam_coord = n1*x1 + n2*x2

  do i3 = 1, nz0

    s0 = (i3 - 1) * this%dz0
    r0 = this%focus_start + s0 * this%focus_dir
    dr(:) = 0.0_p_double
    dr(1) = x1 - r0(1)
    dr(2) = x2 - r0(2)

    u = n1*dr(1) + n2*dr(2)
    v = m1*dr(1) + m2*dr(2)

    ! Required arrival time of this constituent at its waist.
    t_focus = this%focus_t_start + s0 / this%vff

    ! t_envelope(forward) evaluates tau = t-delay-z+xi0.  Passing
    ! z=beam_coord/vg makes the envelope travel along n at vg.  Choose
    ! delay so tau=lon_center at (r0,t_focus).
    delay = t_focus - ( n1*r0(1) + n2*r0(2) ) / this%vg &
            + this%xi0 - this%lon_center()

    this%buf_1d(i3) = gaussian_beam_2d( u, v, this%omega0, this%zr, &
                                         0.0_p_double, 0.0_p_double ) &
                       * this%t_envelope( t - delay, beam_coord / this%vg )

    if( this%focal_weight_type .eq. p_func ) then
      fparser_arr(1) = s0
      this%buf_1d(i3) = this%buf_1d(i3) * eval( this%focal_weight_math_func, fparser_arr )
    endif

    if( this%focal_phase_type .eq. p_func ) then
      fparser_arr(1) = s0
      this%buf_1d(i3) = this%buf_1d(i3) &
        * exp( cmplx( 0.0_p_double, eval( this%focal_phase_math_func, fparser_arr ), kind=p_double ) )
    endif

  enddo

  env = integrate_astrl_samples( this%buf_1d, nz0, this%dz0 ) * this%env_norm

end subroutine get_env_lab_2d_zpulse_flyfoc
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine get_field_lab_2d_zpulse_flyfoc( this, t, x1, x2, field )

  class( t_zpulse_flyfoc ), intent(inout) :: this
  real(p_double), intent(in) :: t, x1, x2
  complex(p_double), intent(out) :: field

  real(p_double) :: beam_coord, fast_phase
  complex(p_double) :: env

  beam_coord = this%carrier_dir(1)*x1 + this%carrier_dir(2)*x2
  fast_phase = this%omega0 * ( t - beam_coord ) - this%phase0

  call this%get_env_lab_2d( t, x1, x2, env )
  field = env * exp( cmplx( 0.0_p_double, fast_phase, kind=p_double ) )

end subroutine get_field_lab_2d_zpulse_flyfoc
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine get_env_lab_3d_zpulse_flyfoc( this, t, x1, x2, x3, env )

  ! Fully vectorial 3D laboratory flying-focus envelope.  The carrier
  ! direction and focal trajectory are independent.  The Gaussian
  ! constituent is evaluated in the orthonormal beam frame (n,p,s).

  class( t_zpulse_flyfoc ), intent(inout) :: this
  real(p_double), intent(in) :: t, x1, x2, x3
  complex(p_double), intent(out) :: env

  integer :: i3, nz0
  real(p_double) :: s0, delay, t_focus, beam_coord
  real(p_double) :: u, v1, v2
  real(p_double), dimension(3) :: r, r0, dr
  real(p_k_fparse), dimension(1) :: fparser_arr

  nz0 = this%integral_size
  r = (/ x1, x2, x3 /)
  beam_coord = dot_product( this%carrier_dir, r )

  do i3 = 1, nz0

    s0 = (i3 - 1) * this%dz0
    r0 = this%focus_start + s0 * this%focus_dir
    dr = r - r0

    u  = dot_product( this%carrier_dir, dr )
    v1 = dot_product( this%lab_basis_p, dr )
    v2 = dot_product( this%lab_basis_s, dr )

    t_focus = this%focus_t_start + s0 / this%vff
    delay = t_focus - dot_product( this%carrier_dir, r0 ) / this%vg &
            + this%xi0 - this%lon_center()

    this%buf_1d(i3) = gaussian_beam_3d( u, v1, v2, this%omega0, this%zr, &
                                         0.0_p_double, 0.0_p_double, 0.0_p_double ) &
                       * this%t_envelope( t - delay, beam_coord / this%vg )

    if( this%focal_weight_type .eq. p_func ) then
      fparser_arr(1) = s0
      this%buf_1d(i3) = this%buf_1d(i3) * eval( this%focal_weight_math_func, fparser_arr )
    endif

    if( this%focal_phase_type .eq. p_func ) then
      fparser_arr(1) = s0
      this%buf_1d(i3) = this%buf_1d(i3) &
        * exp( cmplx( 0.0_p_double, eval( this%focal_phase_math_func, fparser_arr ), kind=p_double ) )
    endif

  enddo

  env = integrate_astrl_samples( this%buf_1d, nz0, this%dz0 ) * this%env_norm

end subroutine get_env_lab_3d_zpulse_flyfoc
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine get_field_lab_3d_zpulse_flyfoc( this, t, x1, x2, x3, field )

  class( t_zpulse_flyfoc ), intent(inout) :: this
  real(p_double), intent(in) :: t, x1, x2, x3
  complex(p_double), intent(out) :: field

  real(p_double) :: beam_coord, fast_phase
  real(p_double), dimension(3) :: r
  complex(p_double) :: env

  r = (/ x1, x2, x3 /)
  beam_coord = dot_product( this%carrier_dir, r )
  fast_phase = this%omega0 * ( t - beam_coord ) - this%phase0

  call this%get_env_lab_3d( t, x1, x2, x3, env )
  field = env * exp( cmplx( 0.0_p_double, fast_phase, kind=p_double ) )

end subroutine get_field_lab_3d_zpulse_flyfoc
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine get_env_astrl_2d_zpulse_flyfoc( this, t, z, x2, env )

  class( t_zpulse_flyfoc ), intent(inout) :: this
  real(p_double), intent(in) :: t, z, x2 
  complex(p_double), intent(out) :: env 

  integer :: i3, nz0 
  real(p_double) :: z0, delay, xperp0, prop_sign, z_beam
  real(p_k_fparse), dimension(1) :: fparser_arr 

  nz0 = this%integral_size

  if( this%propagation .eq. p_forward ) then
    prop_sign = 1.0_p_double
  else
    prop_sign = -1.0_p_double
  endif

  do i3 = 1, nz0
    
    z0 = this%z1 + (i3 - 1) * this%dz0

    delay = ( z0 - this%z_lead ) * ( this%vg/this%vff - prop_sign )

    if( this%focal_xperp_type .eq. p_func ) then 
      fparser_arr(1) = z0 
      xperp0 = eval( this%focal_xperp_math_func, fparser_arr ) 
    else
      xperp0 = this%per_center(1)
    endif 

    ! Reverse the paraxial propagation coordinate for a backward beam.
    z_beam = z0 + prop_sign * ( z - z0 )
    this%buf_1d(i3) = gaussian_beam_2d( z_beam, x2, this%omega0, this%zr, z0, xperp0 ) &
                   * this%t_envelope( t - delay, z )

    if( this%focal_weight_type .eq. p_func ) then 
      fparser_arr(1) = z0 
      this%buf_1d(i3) = this%buf_1d(i3) * eval( this%focal_weight_math_func, fparser_arr )
    endif 

    if( this%focal_phase_type .eq. p_func ) then 
      fparser_arr(1) = z0 
      this%buf_1d(i3) = this%buf_1d(i3) * exp( cmplx( 0, eval( this%focal_phase_math_func, fparser_arr ) ) )
    endif 

  enddo

  env = integrate_astrl_samples( this%buf_1d, nz0, this%dz0 )

  env = env * this%env_norm 

end subroutine get_env_astrl_2d_zpulse_flyfoc
!-----------------------------------------------------------------------------------------




!-----------------------------------------------------------------------------------------
subroutine get_env_astrl_3d_zpulse_flyfoc( this, t, z, x2, x3, env )

  class( t_zpulse_flyfoc ), intent(inout) :: this
  real(p_double), intent(in) :: t, z, x2, x3
  complex(p_double), intent(out) :: env 

  integer :: i3, nz0 
  real(p_double) :: z0, delay, xperp0, yperp0, prop_sign, z_beam
  real(p_k_fparse), dimension(1) :: fparser_arr 

  nz0 = this%integral_size

  if( this%propagation .eq. p_forward ) then
    prop_sign = 1.0_p_double
  else
    prop_sign = -1.0_p_double
  endif

  do i3 = 1, nz0
    
    z0 = this%z1 + (i3 - 1) * this%dz0

    delay = ( z0 - this%z_lead ) * ( this%vg/this%vff - prop_sign ) 

    if( this%focal_xperp_type .eq. p_func ) then 
      fparser_arr(1) = z0 
      xperp0 = eval( this%focal_xperp_math_func, fparser_arr ) 
      yperp0 = this%per_center(2)
    else
      xperp0 = this%per_center(1)
      yperp0 = this%per_center(2)
    endif 

    z_beam = z0 + prop_sign * ( z - z0 )
    this%buf_1d(i3) = gaussian_beam_3d( z_beam, x2, x3, this%omega0, this%zr, z0, xperp0, yperp0 ) &
                   * this%t_envelope( t - delay, z ) 

    if( this%focal_weight_type .eq. p_func ) then 
      fparser_arr(1) = z0 
      this%buf_1d(i3) = this%buf_1d(i3) * eval( this%focal_weight_math_func, fparser_arr )
    endif 

    if( this%focal_phase_type .eq. p_func ) then 
      fparser_arr(1) = z0 
      this%buf_1d(i3) = this%buf_1d(i3) * exp( cmplx( 0, eval( this%focal_phase_math_func, fparser_arr ) ) )
    endif 

  enddo

  env = integrate_astrl_samples( this%buf_1d, nz0, this%dz0 )

  env = env * this%env_norm 

end subroutine get_env_astrl_3d_zpulse_flyfoc
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine get_env_astrl_laguerre_3d_zpulse_flyfoc( this, t, z, x2, x3, env )

  class( t_zpulse_flyfoc ), intent(inout) :: this
  real(p_double), intent(in) :: t, z, x2, x3
  complex(p_double), intent(out) :: env 

  integer :: i3, nz0 
  real(p_double) :: z0, delay, xperp0, yperp0, prop_sign, z_beam
  real(p_k_fparse), dimension(1) :: fparser_arr 

  nz0 = this%integral_size

  if( this%propagation .eq. p_forward ) then
    prop_sign = 1.0_p_double
  else
    prop_sign = -1.0_p_double
  endif

  do i3 = 1, nz0
    
    z0 = this%z1 + (i3 - 1) * this%dz0

    delay = ( z0 - this%z_lead ) * ( this%vg/this%vff - prop_sign ) 

    if( this%focal_xperp_type .eq. p_func ) then 
      fparser_arr(1) = z0 
      xperp0 = eval( this%focal_xperp_math_func, fparser_arr ) 
      yperp0 = this%per_center(2)
    else
      xperp0 = this%per_center(1)
      yperp0 = this%per_center(2)
    endif 

    z_beam = z0 + prop_sign * ( z - z0 )
    this%buf_1d(i3) = laguerre_beam_3d( z_beam, x2, x3, this%omega0, this%zr, this%per_w0(1), z0, &
                                        xperp0, yperp0, &
                                        this%per_tem_mode(2), this%per_tem_mode(1) ) &
                   * this%t_envelope( t - delay, z ) 

    if( this%focal_weight_type .eq. p_func ) then 
      fparser_arr(1) = z0 
      this%buf_1d(i3) = this%buf_1d(i3) * eval( this%focal_weight_math_func, fparser_arr )
    endif 

    if( this%focal_phase_type .eq. p_func ) then 
      fparser_arr(1) = z0 
      this%buf_1d(i3) = this%buf_1d(i3) * exp( cmplx( 0, eval( this%focal_phase_math_func, fparser_arr ) ) )
    endif 

  enddo

  env = integrate_astrl_samples( this%buf_1d, nz0, this%dz0 )

  env = env * this%env_norm 

end subroutine get_env_astrl_laguerre_3d_zpulse_flyfoc
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine get_env_astrl_discrete_2d_zpulse_flyfoc( this, t, z, x2, env )

  class( t_zpulse_flyfoc ), intent(inout) :: this
  real(p_double), intent(in) :: t, z, x2 
  complex(p_double), intent(out) :: env 

  real(p_double) :: z0, delay, weight, phase, xperp0 
  integer :: i 

  env = 0.0_p_double

  do i = 1, this%astrl_discrete_num_pulses

    z0 = this%astrl_discrete_focal_points(i)
    delay = this%astrl_discrete_focal_delays(i)
    weight = this%astrl_discrete_focal_weights(i)
    phase = this%astrl_discrete_focal_phases(i)
    xperp0 = this%astrl_discrete_focal_xperps(i)

    env = env + weight * gaussian_beam_2d( z, x2, this%omega0, this%zr, z0, this%per_center(1) ) &
                        * this%t_envelope( t - delay, z ) &
                        * exp( cmplx( 0, phase ) )

  enddo

  env = env * this%env_norm 

end subroutine get_env_astrl_discrete_2d_zpulse_flyfoc
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine get_env_astrl_discrete_3d_zpulse_flyfoc( this, t, z, x2, x3, env )

  class( t_zpulse_flyfoc ), intent(inout) :: this
  real(p_double), intent(in) :: t, z, x2, x3
  complex(p_double), intent(out) :: env 

  real(p_double) :: z0, delay, weight, phase, xperp0 
  integer :: i 

  env = 0.0_p_double

  do i = 1, this%astrl_discrete_num_pulses

    z0 = this%astrl_discrete_focal_points(i)
    delay = this%astrl_discrete_focal_delays(i)
    weight = this%astrl_discrete_focal_weights(i)
    phase = this%astrl_discrete_focal_phases(i)
    xperp0 = this%astrl_discrete_focal_xperps(i)

    env = env + weight * gaussian_beam_3d( z, x2, x3, this%omega0, this%zr, z0, &
                                           this%per_center(1), this%per_center(1) ) &
                        * this%t_envelope( t - delay, z ) &
                        * exp( cmplx( 0, phase ) )

  enddo

  env = env * this%env_norm 

end subroutine get_env_astrl_discrete_3d_zpulse_flyfoc
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine launch_wall_2d_zpulse_flyfoc( this, b, e, g_space, nx_p_min, t_sim, dt, no_co )

  implicit none

  integer, parameter :: rank = 2

  class( t_zpulse_flyfoc ), intent(inout) :: this
  type( t_vdf ), intent(inout) :: b, e
  type( t_space ), intent(in) :: g_space
  integer, dimension(:), intent(in) :: nx_p_min
  real(p_double), intent(in) :: t_sim, dt
  class( t_node_conf ), intent(in) :: no_co

  ! local variables
  real(p_double), dimension(2,rank) :: g_x_range
  real(p_double), dimension(rank) :: ldx

  real(p_double) :: zmin, rmin, zmax  
  real(p_double) :: z_sim, z_sim_2, t, t_2, r, r_2, z, z_2, fast_phase, fast_phase_2
  complex(p_double) :: field_2_0, field_0_0, field_0_2, field_2_2

  complex(p_double) :: env 

  real(p_double) :: prop_sign
  complex(p_double) :: cos_pol, sin_pol
  real(p_double) :: gpos

  integer :: gipos, ipos, i1, i2

  integer :: lbnd, rbnd, i1start, i1finish
  logical :: lbnd_in_domain, rbnd_in_domain, domain_enclosed, lbnd_right_of_domain

  ! launch at t=0 
  if( this%if_launch .and. this%if_launch_t0 ) then 
    call this % launch_wall_t0_2d( b, e, g_space, nx_p_min, t_sim, dt, no_co )
    this%if_launch = .false.
    return
  endif

  ! Turn off antenna if pulse has ended
  if ( t_sim > this%t_duration() ) then
    SCR_ROOT( "duration exceeded. " )
    SCR_ROOT( "this%t_duration(): ", this%t_duration() )
    SCR_ROOT( "t_sim: ", t_sim )
    this%if_launch = .false.
    return
  endif

    ! Global box sizes
  call get_x_bnd( g_space, g_x_range )

  ldx(1:2) = real( b%dx_(1:2), p_double )
  zmin = real( g_x_range(p_lower,1), p_double )
  zmax = real( g_x_range(p_upper,1), p_double )
  rmin = real( g_x_range(p_lower,2), p_double )
  ! rmax = real( g_x_range(p_upper,2), p_double )

  select case( this%interpolation )
  case( p_linear, p_cubic )
    ! Nothing to change here
  case( p_quadratic, p_quartic )
    zmin = zmin + ldx(1)*0.5_p_double
    zmax = zmax + ldx(1)*0.5_p_double
    rmin = rmin + ldx(2)*0.5_p_double
  end select

  cos_pol = cos( this%pol )
  sin_pol = sin( this%pol )

  if (this%pol_type /= 0) then
    sin_pol = sin_pol * cmplx( 0, this%pol ) 
  endif

  ! Current wall position
  gpos = this%wall_pos + this%wall_vel * t_sim
  gipos = int((gpos - zmin)/ldx(1)) + 1 
  ipos = gipos - nx_p_min(1) + 1

  if( (gpos - ldx(1) < zmin) .or. (gpos + (this%ncells + 1) * ldx(1) > zmax ) ) then 
    SCR_ROOT("WARNING: zpulse_mov_wall hit the boundary before pulse injection was complete. Stopping injection.")
    this%if_launch = .false.
    return
  endif

  ! left and right index bounds of the calculation in local coordinates 
  if( this%propagation == p_forward ) then 
    lbnd = ipos 
    rbnd = ipos + this%ncells
    prop_sign = 1.0_p_double
  else
    lbnd = ipos - this%ncells 
    rbnd = ipos 
    prop_sign = -1.0_p_double
  endif

  lbnd_in_domain  = (lbnd >= 1) .and. (lbnd  <= e%nx_(1))
  lbnd_right_of_domain  = (lbnd  >= e%nx_(1))
  rbnd_in_domain = (rbnd >= 1) .and. (rbnd <= e%nx_(1))
  domain_enclosed  = (lbnd < 1)  .and. (rbnd >  e%nx_(1))

  ! coordinates for local field setting in this domain 
  lbnd = max( lbnd, 1 )
  rbnd = min( rbnd, e%nx_(1) )

  if( this%if_interp_env .and. .not. associated(this%e_env_buffer_2d) ) then 
    call this % precompute_envelope_2d( e, g_space, nx_p_min, dt )
  endif

  if( lbnd_in_domain .or. rbnd_in_domain .or. domain_enclosed ) then  

    i1start = lbnd
    i1finish = rbnd

    if( .not. this%if_add_fields ) then 
      e%f2(:,i1start:i1finish,:) = 0 
      b%f2(:,i1start:i1finish,:) = 0 
    endif

    do i2 = lbound(e%f2, 3), ubound(e%f2, 3)

      r = rmin + (i2 - 1 + nx_p_min(2)) * ldx(2)
      r_2 = r + 0.5_p_double * ldx(2)
   
      do i1 = i1start, i1finish 
     
        z_sim = zmin + real( i1 - 1 + nx_p_min(1), p_double ) * ldx(1)
        z_sim_2 = z_sim + 0.5_p_double * ldx(1)

        ! convert to lab coords 
        if( this % if_boost ) then 
          call this % apply_lorentz_transformation( z_sim, t_sim, z, t )
          call this % apply_lorentz_transformation( z_sim_2, t_sim, z_2, t_2 )
        else 
          z = z_sim
          z_2 = z_sim_2
          t = t_sim 
          t_2 = t_sim
        endif 

        fast_phase = this%omega0 * ( t - prop_sign * z ) - prop_sign * this%phase0
        fast_phase_2 = this%omega0 * ( t_2 - prop_sign * z_2 ) - prop_sign * this%phase0

        ! z, r_2
        call this%get_env_2d( t, z, r_2, env )
        field_0_2 = env * exp( cmplx( 0, fast_phase ) )

        ! z, r 
        call this%get_env_2d( t, z, r, env ) 
        field_0_0 = env * exp( cmplx( 0, fast_phase ) )

        ! z_2, r 
        call this%get_env_2d( t_2, z_2, r, env )
        field_2_0 = env * exp( cmplx( 0, fast_phase_2 ) )

        ! z_2, r_2
        call this%get_env_2d( t_2, z_2, r_2, env )
        field_2_2 = env * exp( cmplx( 0, fast_phase_2 ) )

        e%f2(2,i1,i2) = e%f2(2,i1,i2) + cos_pol * real( field_0_2 ) 
        
        e%f2(3,i1,i2) = e%f2(3,i1,i2) + sin_pol * real( field_0_0 ) 
        
        b%f2(2,i1,i2) = b%f2(2,i1,i2) - sin_pol * real( field_2_0 ) * prop_sign 

        b%f2(3,i1,i2) = b%f2(3,i1,i2) + cos_pol * real( field_2_2 ) * prop_sign 

      enddo 

      ! ! set magnetic fields using the zpulse bint method 
      ! do i1 = ipos, ipos + this % ncells - 1  

      !   ! b%f2(1,i1,i2) =  0.0_p_k_fld
      !   b%f2(2,i1,i2) =  - 0.5_p_k_fld * (e%f2(3,i1,i2) + e%f2(3,i1+1,i2)) * prop_sign
      !   b%f2(3,i1,i2) =  + 0.5_p_k_fld * (e%f2(2,i1,i2) + e%f2(2,i1+1,i2)) * prop_sign

      ! enddo

    enddo

  endif

  if( this%if_clean_fields ) then 

    if( lbnd_in_domain ) then 
      e%f2(:,:i1start-1,:) = 0 
      b%f2(:,:i1start-1,:) = 0 

    else if( lbnd_right_of_domain ) then 
      e%f2(:,:,:) = 0 
      b%f2(:,:,:) = 0        

    endif
  endif

end subroutine launch_wall_2d_zpulse_flyfoc
!-----------------------------------------------------------------------------------------




!-----------------------------------------------------------------------------------------
subroutine launch_wall_t0_2d_zpulse_flyfoc( this, b, e, g_space, nx_p_min, t_sim, dt, no_co )

  implicit none

  integer, parameter :: rank = 2

  class( t_zpulse_flyfoc ), intent(inout) :: this
  type( t_vdf ), intent(inout) :: b, e
  type( t_space ), intent(in) :: g_space
  integer, dimension(:), intent(in) :: nx_p_min
  real(p_double), intent(in) :: t_sim, dt
  class( t_node_conf ), intent(in) :: no_co

  ! local variables
  real(p_double), dimension(2,rank) :: g_x_range
  real(p_double), dimension(rank) :: ldx

  real(p_double) :: zmin, rmin, zmax  
  real(p_double) :: z_sim, z_sim_2, t, t_2, r, r_2, z, z_2, fast_phase, fast_phase_2
  complex(p_double) :: field_2_0, field_0_0, field_0_2, field_2_2

  complex(p_double) :: env 
  integer :: i1, i2 
  real(p_double) :: prop_sign
  complex(p_double) :: cos_pol, sin_pol

  if( this%use_lab_geometry ) then
    call this%launch_lab_t0_2d( b, e, g_space, nx_p_min, t_sim, dt, no_co )
    return
  endif

  ! Global box sizes
  call get_x_bnd( g_space, g_x_range )

  ldx(1:2) = real( b%dx_(1:2), p_double )
  zmin = real( g_x_range(p_lower,1), p_double )
  zmax = real( g_x_range(p_upper,1), p_double )
  rmin = real( g_x_range(p_lower,2), p_double )
  ! rmax = real( g_x_range(p_upper,2), p_double )

  select case( this%interpolation )
  case( p_linear, p_cubic )
    ! Nothing to change here
  case( p_quadratic, p_quartic )
    zmin = zmin + ldx(1)*0.5_p_double
    zmax = zmax + ldx(1)*0.5_p_double
    rmin = rmin + ldx(2)*0.5_p_double
  end select

  cos_pol = cos( this%pol )
  sin_pol = sin( this%pol )

  if (this%pol_type /= 0) then
    sin_pol = sin_pol * cmplx( 0, this%pol ) 
  endif

  ! left and right index bounds of the calculation in local coordinates 
  if( this%propagation == p_forward ) then 
    prop_sign = 1.0_p_double
  else
    prop_sign = -1.0_p_double
  endif

  do i2 = lbound(e%f2, 3), ubound(e%f2, 3)

    r = rmin + (i2 - 1 + nx_p_min(2)) * ldx(2)
    r_2 = r + 0.5_p_double * ldx(2)
 
    do i1 = lbound(e%f2, 2), ubound(e%f2, 2)
   
      z_sim = zmin + real( i1 - 1 + nx_p_min(1), p_double ) * ldx(1)
      z_sim_2 = z_sim + 0.5_p_double * ldx(1)

      ! convert to lab coords 
      if( this % if_boost ) then 
        call this % apply_lorentz_transformation( z_sim, t_sim, z, t )
        call this % apply_lorentz_transformation( z_sim_2, t_sim, z_2, t_2 )
      else 
        z = z_sim
        z_2 = z_sim_2
        t = t_sim 
        t_2 = t_sim
      endif 

      fast_phase = this%omega0 * ( t - prop_sign * z ) - prop_sign * this%phase0
      fast_phase_2 = this%omega0 * ( t_2 - prop_sign * z_2 ) - prop_sign * this%phase0

      ! z, r_2
      call this%get_env_2d( t, z, r_2, env )
      field_0_2 = env * exp( cmplx( 0, fast_phase ) )

      ! z, r 
      call this%get_env_2d( t, z, r, env ) 
      field_0_0 = env * exp( cmplx( 0, fast_phase ) )

      ! z_2, r 
      call this%get_env_2d( t_2, z_2, r, env )
      field_2_0 = env * exp( cmplx( 0, fast_phase_2 ) )

      ! z_2, r_2
      call this%get_env_2d( t_2, z_2, r_2, env )
      field_2_2 = env * exp( cmplx( 0, fast_phase_2 ) )

      e%f2(2,i1,i2) = e%f2(2,i1,i2) + cos_pol * real( field_0_2 ) 
      
      e%f2(3,i1,i2) = e%f2(3,i1,i2) + sin_pol * real( field_0_0 ) 
      
      b%f2(2,i1,i2) = b%f2(2,i1,i2) - sin_pol * real( field_2_0 ) * prop_sign 

      b%f2(3,i1,i2) = b%f2(3,i1,i2) + cos_pol * real( field_2_2 ) * prop_sign 

    enddo 

  enddo

end subroutine launch_wall_t0_2d_zpulse_flyfoc
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine launch_lab_t0_2d_zpulse_flyfoc( this, b, e, g_space, nx_p_min, t_sim, dt, no_co )

  ! Initialize a 2D flying-focus field whose carrier direction is independent
  ! of its focal trajectory.  The field basis is rotated with the carrier:
  !   n = (cos(theta), sin(theta), 0)
  !   p = (-sin(theta), cos(theta), 0)
  !   s = (0,0,1)
  ! and pol mixes p and s.  This reduces to the legacy x1-propagating
  ! convention at theta=0.

  implicit none

  integer, parameter :: rank = 2
  class( t_zpulse_flyfoc ), intent(inout) :: this
  type( t_vdf ), intent(inout) :: b, e
  type( t_space ), intent(in) :: g_space
  integer, dimension(:), intent(in) :: nx_p_min
  real(p_double), intent(in) :: t_sim, dt
  class( t_node_conf ), intent(in) :: no_co

  real(p_double), dimension(2,rank) :: g_x_range
  real(p_double), dimension(rank) :: ldx
  real(p_double) :: xmin1, xmin2, x1, x2, x1h, x2h
  real(p_double) :: ct, st, cp, sp
  integer :: i1, i2
  complex(p_double) :: f_e1, f_e2, f_e3, f_b1, f_b2, f_b3

  call get_x_bnd( g_space, g_x_range )
  ldx(1:2) = real( b%dx_(1:2), p_double )
  xmin1 = real( g_x_range(p_lower,1), p_double )
  xmin2 = real( g_x_range(p_lower,2), p_double )

  select case( this%interpolation )
  case( p_linear, p_cubic )
    ! no coordinate correction
  case( p_quadratic, p_quartic )
    xmin1 = xmin1 + 0.5_p_double*ldx(1)
    xmin2 = xmin2 + 0.5_p_double*ldx(2)
  end select

  ct = this%carrier_dir(1)
  st = this%carrier_dir(2)
  cp = cos( this%pol )
  sp = sin( this%pol )

  do i2 = lbound(e%f2,3), ubound(e%f2,3)
    x2 = xmin2 + real( i2 - 1 + nx_p_min(2), p_double ) * ldx(2)
    x2h = x2 + 0.5_p_double*ldx(2)

    do i1 = lbound(e%f2,2), ubound(e%f2,2)
      x1 = xmin1 + real( i1 - 1 + nx_p_min(1), p_double ) * ldx(1)
      x1h = x1 + 0.5_p_double*ldx(1)

      ! Yee locations: E1(x+1/2,y), E2(x,y+1/2), E3(x,y);
      ! B1(x,y+1/2), B2(x+1/2,y), B3(x+1/2,y+1/2).
      call this%get_field_lab_2d( t_sim, x1h, x2,  f_e1 )
      call this%get_field_lab_2d( t_sim, x1,  x2h, f_e2 )
      call this%get_field_lab_2d( t_sim, x1,  x2,  f_e3 )
      call this%get_field_lab_2d( t_sim, x1,  x2h, f_b1 )
      call this%get_field_lab_2d( t_sim, x1h, x2,  f_b2 )
      call this%get_field_lab_2d( t_sim, x1h, x2h, f_b3 )

      ! E = cp*p + sp*s, with p=(-st,ct,0).
      e%f2(1,i1,i2) = e%f2(1,i1,i2) - st*cp*real(f_e1,p_k_fld)
      e%f2(2,i1,i2) = e%f2(2,i1,i2) + ct*cp*real(f_e2,p_k_fld)
      e%f2(3,i1,i2) = e%f2(3,i1,i2) +    sp*real(f_e3,p_k_fld)

      ! B = n x E = cp*s - sp*p.
      b%f2(1,i1,i2) = b%f2(1,i1,i2) + st*sp*real(f_b1,p_k_fld)
      b%f2(2,i1,i2) = b%f2(2,i1,i2) - ct*sp*real(f_b2,p_k_fld)
      b%f2(3,i1,i2) = b%f2(3,i1,i2) +    cp*real(f_b3,p_k_fld)

    enddo
  enddo

end subroutine launch_lab_t0_2d_zpulse_flyfoc
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine launch_lab_t0_3d_zpulse_flyfoc( this, b, e, g_space, nx_p_min, t_sim, dt, no_co )

  ! Initialize an arbitrary-orientation 3D laboratory flying-focus field.
  ! The scalar complex beam field is sampled at the native Yee location of
  ! each component; vector amplitudes are then projected using
  ! E = cos(pol)*p + sin(pol)*s and B = n x E.

  implicit none

  integer, parameter :: rank = 3
  class( t_zpulse_flyfoc ), intent(inout) :: this
  type( t_vdf ), intent(inout) :: b, e
  type( t_space ), intent(in) :: g_space
  integer, dimension(:), intent(in) :: nx_p_min
  real(p_double), intent(in) :: t_sim, dt
  class( t_node_conf ), intent(in) :: no_co

  real(p_double), dimension(2,rank) :: g_x_range
  real(p_double), dimension(rank) :: ldx, xmin, x, xh
  real(p_double), dimension(3) :: ehat, bhat
  real(p_double) :: cp, sp
  integer :: i1, i2, i3
  complex(p_double) :: f_e1, f_e2, f_e3, f_b1, f_b2, f_b3

  call get_x_bnd( g_space, g_x_range )
  ldx(1:3) = real( b%dx_(1:3), p_double )
  xmin(1:3) = real( g_x_range(p_lower,1:3), p_double )

  select case( this%interpolation )
  case( p_linear, p_cubic )
    ! no coordinate correction
  case( p_quadratic, p_quartic )
    xmin = xmin + 0.5_p_double*ldx
  end select

  cp = cos( this%pol )
  sp = sin( this%pol )
  ehat = cp*this%lab_basis_p + sp*this%lab_basis_s

  bhat(1) = this%carrier_dir(2)*ehat(3) - this%carrier_dir(3)*ehat(2)
  bhat(2) = this%carrier_dir(3)*ehat(1) - this%carrier_dir(1)*ehat(3)
  bhat(3) = this%carrier_dir(1)*ehat(2) - this%carrier_dir(2)*ehat(1)

  do i3 = lbound(e%f3,4), ubound(e%f3,4)
    x(3) = xmin(3) + real( i3 - 1 + nx_p_min(3), p_double ) * ldx(3)
    xh(3) = x(3) + 0.5_p_double*ldx(3)

    do i2 = lbound(e%f3,3), ubound(e%f3,3)
      x(2) = xmin(2) + real( i2 - 1 + nx_p_min(2), p_double ) * ldx(2)
      xh(2) = x(2) + 0.5_p_double*ldx(2)

      do i1 = lbound(e%f3,2), ubound(e%f3,2)
        x(1) = xmin(1) + real( i1 - 1 + nx_p_min(1), p_double ) * ldx(1)
        xh(1) = x(1) + 0.5_p_double*ldx(1)

        ! Yee locations:
        ! E1(x+1/2,y,z), E2(x,y+1/2,z), E3(x,y,z+1/2)
        ! B1(x,y+1/2,z+1/2), B2(x+1/2,y,z+1/2),
        ! B3(x+1/2,y+1/2,z).
        call this%get_field_lab_3d( t_sim, xh(1), x(2),  x(3),  f_e1 )
        call this%get_field_lab_3d( t_sim, x(1),  xh(2), x(3),  f_e2 )
        call this%get_field_lab_3d( t_sim, x(1),  x(2),  xh(3), f_e3 )
        call this%get_field_lab_3d( t_sim, x(1),  xh(2), xh(3), f_b1 )
        call this%get_field_lab_3d( t_sim, xh(1), x(2),  xh(3), f_b2 )
        call this%get_field_lab_3d( t_sim, xh(1), xh(2), x(3),  f_b3 )

        e%f3(1,i1,i2,i3) = e%f3(1,i1,i2,i3) + ehat(1)*real(f_e1,p_k_fld)
        e%f3(2,i1,i2,i3) = e%f3(2,i1,i2,i3) + ehat(2)*real(f_e2,p_k_fld)
        e%f3(3,i1,i2,i3) = e%f3(3,i1,i2,i3) + ehat(3)*real(f_e3,p_k_fld)

        b%f3(1,i1,i2,i3) = b%f3(1,i1,i2,i3) + bhat(1)*real(f_b1,p_k_fld)
        b%f3(2,i1,i2,i3) = b%f3(2,i1,i2,i3) + bhat(2)*real(f_b2,p_k_fld)
        b%f3(3,i1,i2,i3) = b%f3(3,i1,i2,i3) + bhat(3)*real(f_b3,p_k_fld)

      enddo
    enddo
  enddo

end subroutine launch_lab_t0_3d_zpulse_flyfoc
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine launch_wall_3d_zpulse_flyfoc( this, b, e, g_space, nx_p_min, t_sim, dt, no_co )


  implicit none

  integer, parameter :: rank = 3

  class( t_zpulse_flyfoc ), intent(inout) :: this
  type( t_vdf ), intent(inout) :: b, e
  type( t_space ), intent(in) :: g_space
  integer, dimension(:), intent(in) :: nx_p_min
  real(p_double), intent(in) :: t_sim, dt
  class( t_node_conf ), intent(in) :: no_co

  ! local variables
  real(p_double), dimension(2,rank) :: g_x_range
  real(p_double), dimension(rank) :: ldx

  real(p_double) :: zmin, r1min, r2min, zmax  
  real(p_double) :: z_sim, z_sim_2, t, t_2, r1, r1_2, r2, r2_2, z, z_2, fast_phase, fast_phase_2
  complex(p_double) :: field_2_0_2, field_0_0_2, field_0_2_0, field_2_2_0

  complex(p_double) :: env 

  real(p_double) :: prop_sign
  complex(p_double) :: cos_pol, sin_pol
  real(p_double) :: gpos

  integer :: gipos, ipos, i1, i2, i3 

  integer :: lbnd, rbnd, i1start, i1finish
  logical :: lbnd_in_domain, rbnd_in_domain, domain_enclosed, lbnd_right_of_domain

  if( this%use_lab_geometry ) then
    call this%launch_lab_t0_3d( b, e, g_space, nx_p_min, t_sim, dt, no_co )
    return
  endif
    

  ! Turn off antenna if pulse has ended
  if ( t_sim > this%t_duration() ) then
    SCR_ROOT( "duration exceeded. " )
    SCR_ROOT( "this%t_duration(): ", this%t_duration() )
    SCR_ROOT( "t_sim: ", t_sim )
    this%if_launch = .false.
    return
  endif

    ! Global box sizes
  call get_x_bnd( g_space, g_x_range )

  ldx(1:3) = real( b%dx_(1:3), p_double )
  zmin = real( g_x_range(p_lower,1), p_double )
  zmax = real( g_x_range(p_upper,1), p_double )
  r1min = real( g_x_range(p_lower,2), p_double )
  r2min = real( g_x_range(p_lower,3), p_double )
  ! rmax = real( g_x_range(p_upper,2), p_double )

  select case( this%interpolation )
  case( p_linear, p_cubic )
    ! Nothing to change here
  case( p_quadratic, p_quartic )
    zmin = zmin + ldx(1)*0.5_p_double
    zmax = zmax + ldx(1)*0.5_p_double
    r1min = r1min + ldx(2)*0.5_p_double
    r2min = r2min + ldx(3)*0.5_p_double
  end select

  cos_pol = cos( this%pol )
  sin_pol = sin( this%pol )

  if (this%pol_type /= 0) then
    sin_pol = sin_pol * cmplx( 0, this%pol ) 
  endif

  ! Current wall position
  gpos = this%wall_pos + this%wall_vel * t_sim
  gipos = int((gpos - zmin)/ldx(1)) + 1 
  ipos = gipos - nx_p_min(1) + 1

  if( (gpos - ldx(1) < zmin) .or. (gpos + (this%ncells + 1) * ldx(1) > zmax ) ) then 
    SCR_ROOT("WARNING: zpulse_mov_wall hit the boundary before pulse injection was complete. Stopping injection.")
    this%if_launch = .false.
    return
  endif

  ! left and right index bounds of the calculation in local coordinates 
  if( this%propagation == p_forward ) then 
    lbnd = ipos 
    rbnd = ipos + this%ncells
    prop_sign = 1.0_p_double
  else
    lbnd = ipos - this%ncells 
    rbnd = ipos 
    prop_sign = -1.0_p_double
  endif

  lbnd_in_domain  = (lbnd >= 1) .and. (lbnd  <= e%nx_(1))
  lbnd_right_of_domain  = (lbnd  >= e%nx_(1))
  rbnd_in_domain = (rbnd >= 1) .and. (rbnd <= e%nx_(1))
  domain_enclosed  = (lbnd < 1)  .and. (rbnd >  e%nx_(1))

  ! coordinates for local field setting in this domain 
  lbnd = max( lbnd, 1 )
  rbnd = min( rbnd, e%nx_(1) )

  if( this%if_interp_env .and. .not. associated(this%e_env_buffer_3d) ) then 
    call this % precompute_envelope_3d( e, g_space, nx_p_min, dt )
  endif

  if( lbnd_in_domain .or. rbnd_in_domain .or. domain_enclosed ) then  

    i1start = lbnd
    i1finish = rbnd

    if( .not. this%if_add_fields ) then 
      e%f3(:,i1start:i1finish,:,:) = 0 
      b%f3(:,i1start:i1finish,:,:) = 0 
    endif

    do i2 = lbound(e%f3, 3), ubound(e%f3, 3)

      r1 = r1min + (i2 - 1 + nx_p_min(2)) * ldx(2)
      r1_2 = r1 + 0.5_p_double * ldx(2)
   
      do i3 = lbound(e%f3, 4), ubound(e%f3, 4)

        r2 = r2min + (i3 - 1 + nx_p_min(3)) * ldx(3)
        r2_2 = r2 + 0.5_p_double * ldx(3)
     
        do i1 = i1start, i1finish 
       
          z_sim = zmin + real( i1 - 1 + nx_p_min(1), p_double ) * ldx(1)
          z_sim_2 = z_sim + 0.5_p_double * ldx(1)

          ! convert to lab coords 
          if( this % if_boost ) then 
            call this % apply_lorentz_transformation( z_sim, t_sim, z, t )
            call this % apply_lorentz_transformation( z_sim_2, t_sim, z_2, t_2 )
          else 
            z = z_sim
            z_2 = z_sim_2
            t = t_sim 
            t_2 = t_sim
          endif 

          fast_phase = this%omega0 * ( t - prop_sign * z ) - prop_sign * this%phase0
          fast_phase_2 = this%omega0 * ( t_2 - prop_sign * z_2 ) - prop_sign * this%phase0

          call this%get_env_3d( t, z, r1, r2, env )
          field_0_2_0 = env * exp( cmplx( 0, fast_phase ) )

          ! z, r 
          call this%get_env_3d( t, z, r1, r2_2, env ) 
          field_0_0_2 = env * exp( cmplx( 0, fast_phase ) )

          ! z_2, r 
          call this%get_env_3d( t_2, z_2, r1, r2_2, env )
          field_2_0_2 = env * exp( cmplx( 0, fast_phase_2 ) )

          ! z_2, r_2
          call this%get_env_3d( t_2, z_2, r1_2, r2, env )
          field_2_2_0 = env * exp( cmplx( 0, fast_phase_2 ) )

          e%f3(2,i1,i2,i3) = e%f3(2,i1,i2,i3) + cos_pol * real( field_0_2_0 ) 
          
          e%f3(3,i1,i2,i3) = e%f3(3,i1,i2,i3) + sin_pol * real( field_0_0_2 ) 
          
          b%f3(2,i1,i2,i3) = b%f3(2,i1,i2,i3) - sin_pol * real( field_2_0_2 ) * prop_sign 

          b%f3(3,i1,i2,i3) = b%f3(3,i1,i2,i3) + cos_pol * real( field_2_2_0 ) * prop_sign 

        enddo

      enddo 

    enddo

  endif


  if( this%if_clean_fields ) then 

    if( lbnd_in_domain ) then 
      e%f3(:,:i1start-1,:,:) = 0 
      b%f3(:,:i1start-1,:,:) = 0 

    else if( lbnd_right_of_domain ) then 
      e%f3(:,:,:,:) = 0 
      b%f3(:,:,:,:) = 0        

    endif
  endif

end subroutine launch_wall_3d_zpulse_flyfoc
!-----------------------------------------------------------------------------------------



end module m_zpulse_flyfoc


