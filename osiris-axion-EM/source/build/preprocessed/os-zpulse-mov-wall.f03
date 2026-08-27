# 1 "zpulse/os-zpulse-mov-wall.f03"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "zpulse/os-zpulse-mov-wall.f03"
# 1 "./os-config.h" 1
! Configuration file for osiris

! ----------------------------------------------------------------------------------------
! Algorithm options
! ----------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------------
! System options
! ----------------------------------------------------------------------------------------

! MPI supports MPI_IN_PLACE in global operations
!#define __HAS_MPI_IN_PLACE__

! Use OpenMP
!#define _OPENMP

! Use 1 optimized code
!#define 1
!#define 1
!#define SIMD_AVX
!#define SIMD_BGQ
!#define SIMD_MIC

! Use SION for checkpointing (not available on all systems)
!#define __RST_IO__ = __RST_SION__

! Use log files
!#define __USE_LOG__

! Use MPE for logging and profiling
!#define __USE_MPE__

! Compiler does not fully support the sizeof intrinsic which is a Fortran 2003 feature
!#define __NO_SIZEOF__

! Use the PAPI library for profiling
!#define __USE_PAPI__


! ----------------------------------------------------------------------------------------
! Distribution Options
! ----------------------------------------------------------------------------------------
!
! Optional modules to be removed for distribution
! These are all turned on by default unless the __DISTRO__ preprocessor macro is defined
!
! ----------------------------------------------------------------------------------------



! Include Ionization module


! Include binary collisions module


! Include particle tracking module


! Include perfectly matched layers boundary conditions for EMF


! Include spin advance module which is disabled by default
! #define __HAS_SPIN__

! Include a debug flag for cylindrical modes simulations
! #define __CYL_MODES_DEBUG__

! Compile functions for reporting B fields in the RaDiO module
! #define __HAS_RAD_BFLD__
# 2 "zpulse/os-zpulse-mov-wall.f03" 2
# 1 "./os-preprocess.fpp" 1
!
! File: os-preprocess.fpp
!
! A set of preprocessing macros for osiris, and the fpp/cpp preprocessor
!
!



! Macros for the IBM CPP/GNU preprocessor




! Assertion macro



! Debug functions
# 33 "./os-preprocess.fpp"
! SCR hostname functions







! LOG functions




! ERROR functions






! WARNING functions




! functions for restart io
# 3 "zpulse/os-zpulse-mov-wall.f03" 2

module m_zpulse_mov_wall

use m_system
use m_parameters
use m_math

use m_zpulse_std
use m_zpulse_wall

use m_node_conf
use m_space
use m_vdf_define
use m_grid_define

use m_fparser

use m_emf_define, only: t_emf_bound, t_emf
use m_emf_bound, only: is_open

implicit none

private


!-----------------------------------------------------------------------------------------
type, extends(t_zpulse) :: t_zpulse_mov_wall

  ! whether to initialize fields fully at t=0
  logical :: if_launch_t0

  integer :: tenv_type

  real(p_double) :: inj_time

  real(p_double) :: tenv_rise, tenv_flat, tenv_fall
  real(p_double) :: tenv_duration, tenv_range, tenv_launch_duration
  type(t_fparser) :: tenv_math_func
  real(p_double) :: xi0
  integer :: ncells

  logical :: if_tenv_tilt
  real(p_double), dimension(2) :: tenv_tilt

  real(p_double) :: wall_pos ! initial wall position
  real(p_double) :: wall_vel ! velocity (must be >= 0 and <= 1)
  logical :: if_launch_from_wall

  logical :: if_clean_fields
  logical :: if_clean_x2
  real(p_double) :: min_x2_clean, max_x2_clean

  real(p_double) :: clean_fields_delta
  real(p_double) :: clean_fields_time

  logical :: if_apply_divcorr
  logical :: if_add_fields

contains

  procedure :: init => init_zpulse_mov_wall
  procedure :: read_input => read_input_zpulse_mov_wall
  procedure :: launch => launch_zpulse_mov_wall
  procedure :: t_duration => t_duration_zpulse_mov_wall
  procedure :: t_envelope => t_env_zpulse_mov_wall
  ! procedure :: clean_fields => clean_fields_zpulse_mov_wall
  procedure :: lon_center => lon_center_zpulse_mov_wall

  ! procedure :: launch_wall_1d => launch_wall_1d_zpulse_mov_wall
  procedure :: launch_wall_2d => launch_wall_2d_zpulse_mov_wall
  procedure :: launch_wall_3d => launch_wall_3d_zpulse_mov_wall

  procedure :: apply_divcorr_2d => apply_divcorr_2d_zpulse_mov_wall
  procedure :: apply_divcorr_3d => apply_divcorr_3d_zpulse_mov_wall

  ! procedure :: clean_fields_2d => clean_fields_2d_zpulse_mov_wall
  ! procedure :: clean_fields_3d => clean_fields_3d_zpulse_mov_wall


end type t_zpulse_mov_wall
!-----------------------------------------------------------------------------------------



public :: t_zpulse_mov_wall

contains

!-----------------------------------------------------------------------------------------
subroutine init_zpulse_mov_wall( this, restart, t, dt, b, g_space, nx_p_min, no_co, interpolation )

  implicit none

  class( t_zpulse_mov_wall ), intent(inout) :: this
  logical, intent(in) :: restart
  real( p_double ), intent(in) :: t, dt
  type( t_vdf ) , intent(in) :: b
  type( t_space ), intent(in) :: g_space
  integer, intent(in), dimension(:) :: nx_p_min
  class( t_node_conf ), intent(in) :: no_co
  integer, intent(in) :: interpolation

  logical :: if_launch

  ! Need to preserve if_launch value for this zpulse
  if_launch = this%if_launch

  call this % t_zpulse % init( restart, t, dt, b, g_space, nx_p_min, no_co, interpolation )

  this%if_launch = if_launch

  if ( this%if_launch .and. ( t > this%launch_time + this % t_duration() ) ) then
    this%if_launch = .false.
  endif

end subroutine init_zpulse_mov_wall
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine read_input_zpulse_mov_wall( this, input_file, g_space, bnd_con, periodic, grid, sim_options )

  use m_input_file

  implicit none

  class( t_zpulse_mov_wall ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  type( t_space ), intent(in) :: g_space
  class (t_emf_bound), intent(in) :: bnd_con
  logical, dimension(:), intent(in) :: periodic
  class( t_grid ), intent(in) :: grid
  type( t_options ), intent(in) :: sim_options

  real(p_double) :: a0
  real(p_double) :: omega0
  real(p_double) :: phase
  integer :: pol_type
  real(p_double) :: pol
  character(len=16) :: propagation

  integer :: chirp_order
  real(p_double), dimension(p_max_chirp_order) :: chirp_coefs

  logical :: if_launch_t0

  ! Temporal profile
  character(len=16) :: tenv_type

  real(p_double) :: tenv_fwhm

  real(p_double) :: tenv_rise, tenv_flat, tenv_fall
  real(p_double) :: tenv_duration, tenv_range, tenv_launch_duration
  real(p_double) :: xi0

  character(len = p_max_expr_len) :: tenv_math_func
  real(p_double), dimension(2) :: tenv_tilt

  ! Beam profile
  character(len=16) :: per_type
  real(p_double), dimension(2) :: per_center

  integer :: per_chirp_order
  real(p_double), dimension(p_max_chirp_order) :: per_chirp_coefs

  real(p_double), dimension(2) :: per_w0, per_fwhm, per_focus
  integer, dimension(2) :: per_tem_mode

  real(p_double), dimension(2,2) :: per_w0_asym, per_fwhm_asym, per_focus_asym
  real(p_double), dimension(2) :: per_asym_trans

  integer :: per_n, per_0clip
  real(p_double) :: per_kt

  real(p_double) :: launch_time
  logical :: if_launch

  real(p_double) :: wall_pos, wall_vel, clean_fields_delta
  integer :: ncells
  logical :: if_clean_fields, if_launch_from_wall, if_apply_divcorr, if_add_fields

  namelist /nl_zpulse_mov_wall/ &
    a0, omega0,phase, pol_type, pol, propagation, &
    chirp_order, chirp_coefs, xi0, tenv_fwhm, &
    if_launch_t0, tenv_type, tenv_rise, tenv_flat, tenv_fall, &
    tenv_duration, tenv_range, tenv_launch_duration, tenv_math_func, tenv_tilt, &
    per_type, per_center, per_w0, per_fwhm, per_focus, &
    per_chirp_order, per_chirp_coefs, &
    per_w0_asym, per_fwhm_asym, per_focus_asym, per_asym_trans, &
    per_n, per_0clip, per_kt, per_tem_mode, &
    launch_time, if_launch, wall_pos, wall_vel, if_launch_from_wall, ncells, &
    if_clean_fields, clean_fields_delta, if_apply_divcorr, if_add_fields

  integer :: i, ierr

  if (disp_out(input_file)) then
    if (mpi_node()==0) print *," - Reading zpulse_mov_wall configuration..."
  endif

  if_launch = .true.

  a0 = 1.0_p_double
  omega0 = 10.0_p_double
  phase = 0.0_p_double
  pol_type = 0
  pol = 90.0_p_double

  chirp_order = 0
  chirp_coefs = 0.0_p_double

  propagation = "forward"

  if_launch_t0 = .false.
  tenv_type = "polynomial"

  tenv_fwhm = -1.0_p_double

  tenv_rise = 0.0_p_double
  tenv_flat = 0.0_p_double
  tenv_fall = 0.0_p_double

  xi0 = 0.0_p_double


  tenv_duration = 0.0_p_double
  tenv_range = -huge(1.0_p_double)
  tenv_launch_duration = -huge(1.0_p_double)
  tenv_math_func = "NO_FUNCTION_SUPPLIED!"

  tenv_tilt = 0.0_p_double

  per_type = "plane"
  per_center = -huge(1.0_p_double)

  per_chirp_order = 0
  per_chirp_coefs = 0.0_p_double

  per_w0 = 0.0_p_double
  per_fwhm = 0.0_p_double
  per_focus = -huge(1.0_p_double)
  per_tem_mode = 0 ! tem mode, integer, default = 0,0

  per_w0_asym = 0.0_p_double
  per_fwhm_asym = 0.0_p_double
  per_focus_asym = -huge(1.0_p_double)
  per_asym_trans = 0.0_p_double

  per_n = 0
  per_0clip = 1
  per_kt = 0.0_p_double

  wall_vel = 0.0_p_double
  wall_pos = -huge(1.0d0)
  ncells = 4
  if_clean_fields = .false.
  clean_fields_delta = 0.0_p_double
  if_launch_from_wall = .false.

  if_apply_divcorr = .true.
  if_add_fields = .false.


  launch_time = 0.0_p_double


  ! read data from file
  call get_namelist( input_file, "nl_zpulse_mov_wall", ierr )

  if ( ierr /= 0 ) then
    if ( mpi_node() == 0 ) then
      if (ierr < 0) then
        print *, "Error reading zpulse_mov_wall parameters"
      else
        print *, "Error: zpulse_mov_wall parameters missing"
      endif
      print *, "aborting..."
    endif
    stop
  endif

  read (input_file%nml_text, nml = nl_zpulse_mov_wall, iostat = ierr)
  if (ierr /= 0) then
    if ( mpi_node() == 0 ) then
      write(0,*) ""
      write(0,*) "   Error reading zpulse_mov_wall parameters"
      write(0,*) "   aborting..."
    endif
    stop
  endif

  this%a0 = a0
  this%omega0 = omega0
  this%phase0 = real( phase * pi_180, p_double )

  ! Polarization parameters
  if ((pol_type <-1) .or. (pol_type >1)) then
    if ( mpi_node() == 0 ) then
      write(0,*) ""
      write(0,*) "   Error in zpulse_mov_wall parameters"
      write(0,*) "   pol must be in the range [-1,0,+1]"
      write(0,*) "   aborting..."
    endif
    stop
  endif
  this%pol_type = pol_type
  this%pol = real( pol * pi_180, p_double )

  this % xi0 = xi0


  ! Frequency chirp parameters
  if ((chirp_order < 0) .or. (chirp_order > p_max_chirp_order)) then
    if ( mpi_node() == 0 ) then
      write(0,*) ""
      write(0,*) "   Error in zpulse_mov_wall parameters"
      write(0,*) "   chirp_order must be in the range [0,",p_max_chirp_order,"]"
      write(0,*) "   aborting..."
    endif
    stop
  endif
  this%chirp_order = chirp_order
  this%chirp_coefs = chirp_coefs

  if ((per_chirp_order < 0) .or. (per_chirp_order > p_max_chirp_order)) then
    if ( mpi_node() == 0 ) then
      write(0,*) ""
      write(0,*) "   Error in zpulse_mov_wall parameters"
      write(0,*) "   per_chirp_order must be in the range [0,",p_max_chirp_order,"]"
      write(0,*) "   aborting..."
    endif
    stop
  endif
  this%per_chirp_order = per_chirp_order
  this%per_chirp_coefs = per_chirp_coefs

  select case ( trim(propagation) )
  case ("forward")
    this%propagation = p_forward
  case ("backward")
    this%propagation = p_backward
  case default
    if ( mpi_node() == 0 ) then
      write(0,*) ''
      write(0,*) '   Error in zpulse parameters'
      write(0,*) '   propagation must be either "forward" or "backward"'
      write(0,*) '   aborting...'
    endif
    stop
  end select

  ! Temporal envelope tilt
  this%if_tenv_tilt = .false.
  if ( p_x_dim > 1 ) then
    do i = 1, p_x_dim-1
      if ( tenv_tilt(i) /= 0.0_p_double ) then
        if (( tenv_tilt(i) < -45.0_p_double ) .or. (tenv_tilt(i) > 45.0_p_double )) then
          if ( mpi_node() == 0 ) then
            write(0,*) ""
            write(0,*) "   Error in zpulse_mov_wall parameters"
            write(0,*) "   tenv_tilt must be in the range [ -45 , 45 ] deg"
            write(0,*) "   aborting..."
          endif
          stop
        endif
        this%tenv_tilt(i) = tan( real( tenv_tilt(i) * pi_180, p_double ))
        if ( this%propagation == p_forward ) this%tenv_tilt(i) = -this%tenv_tilt(i)
        this%if_tenv_tilt = .true.
      else
        this%tenv_tilt(i) = 0.0_p_double
      endif
    enddo
  endif


  ! Temporal envelope parameters
  select case ( trim(tenv_type) )
  case ("polynomial")
    this%tenv_type = p_polynomial

    if ( tenv_fwhm > 0.0_p_double ) then
      this%tenv_rise = tenv_fwhm/2.0_p_double
      this%tenv_flat = 0.0_p_double
      this%tenv_fall = tenv_fwhm/2.0_p_double
    else
      this%tenv_rise = tenv_rise
      this%tenv_flat = tenv_flat
      this%tenv_fall = tenv_fall
    endif

  case ("sin2")
    this%tenv_type = p_sin2

    if ( tenv_fwhm > 0.0_p_double ) then
      this%tenv_rise = tenv_fwhm/2.0_p_double
      this%tenv_flat = 0.0_p_double
      this%tenv_fall = tenv_fwhm/2.0_p_double
    else
      this%tenv_rise = tenv_rise
      this%tenv_flat = tenv_flat
      this%tenv_fall = tenv_fall
    endif

  case ("gaussian")
    this%tenv_type = p_gaussian

    if ( tenv_fwhm > 0.0_p_double ) then
      this%tenv_duration = tenv_fwhm/sqrt( 2.0_p_double * log(2.0_p_double))
    else
      this%tenv_duration = tenv_duration
    endif
    this%tenv_range = tenv_range
    this%tenv_launch_duration = tenv_launch_duration

    if ( tenv_launch_duration == -huge(1.0_p_double) ) then
      if ( mpi_node() == 0 ) then
        write(0,*) ''
        write(0,*) '   Error in zpulse_mov_wall parameters'
        write(0,*) '   When using "gaussian" temporal envelope, '
        write(0,*) '   tenv_launch_duration must also be defined'
        write(0,*) '   aborting...'
      endif
      stop
    endif

  case ("math")
    this%tenv_type = p_func
    this%tenv_range = tenv_range

    if ( tenv_range == -huge(1.0_p_double) ) then
      if ( mpi_node() == 0 ) then
        write(0,*) ''
        write(0,*) '   Error in zpulse_mov_wall parameters'
        write(0,*) '   When using "math" temporal envelopes, '
        write(0,*) '   tenv_range must also be defined'
        write(0,*) '   aborting...'
      endif
      stop
    endif

    call setup(this%tenv_math_func, trim(tenv_math_func), (/'t'/), ierr)
    if (ierr /= 0) then
      if ( mpi_node() == 0 ) then
        write(0,*) ''
        write(0,*) '   Error in zpulse_mov_wall parameters'
        write(0,*) '   Supplied tenv_math_func failed to compile'
        write(0,*) '   aborting...'
      endif
      stop
    endif

  case default
    if ( mpi_node() == 0 ) then
      write(0,*) ''
      write(0,*) '   Error in zpulse_mov_wall parameters'
      write(0,*) '   lon_type must be either "polynomial", "gaussian",'
      write(0,*) '   "sin2", "math" or "hermite"'
      write(0,*) '   aborting...'
    endif
    stop
  end select



  ! perpendiular profile
  select case ( trim(per_type) )
  case ("plane")
    this%per_type = p_plane

  case ("gaussian","hermite")

    this%per_type = p_hermite_gaussian
    if ( p_x_dim > 1 ) then
      this%per_tem_mode = 0
      do i = 1, p_x_dim - 1
        this%per_tem_mode(i) = per_tem_mode(i)
        if ( per_tem_mode(i) > 7 ) then
          if (mpi_node() == 0) then
            write(0,*) ''
            write(0,*) '   Error in zpulse_mov_wall parameters'
            write(0,*) '   Only TEM modes up to 7 have been implemented. Please contact the'
            write(0,*) '   development team if you need higher order modes.'
            write(0,*) '   aborting...'
          endif
          stop
        endif
      enddo
    endif

    ! get main spot size and focal plane
    this%per_w0 = 0.0_p_double
    if (per_fwhm(1) > 0.0_p_double) then
      this%per_w0 = per_fwhm(1)/sqrt(4.0_p_double * log(2.0_p_double))
    else
      this%per_w0 = per_w0(1)
    endif

    if ( this%per_w0(1) <= 0.0_p_double ) then
      if (mpi_node() == 0) then
        write(0,*) ''
        write(0,*) '   Error in zpulse_mov_wall parameters'
        write(0,*) '   for gaussian pulses per_w0 or per_fwhm must be > 0'
        write(0,*) '   aborting...'
      endif
      stop
    endif

    if ( per_focus(1) == -huge(1.0_p_double) ) then
      if (mpi_node() == 0) then
        write(0,*) ''
        write(0,*) '   Error in zpulse_mov_wall parameters'
        write(0,*) '   for gaussian pulses per_focus must be defined'
        write(0,*) '   aborting...'
      endif
      stop
    endif
    this%per_focus = per_focus(1)

    ! check for astigmatic beam
    if (p_x_dim == 3) then
      if (per_fwhm(2) > 0.0_p_double) then
        this%per_w0(2) = per_fwhm(2)/sqrt(4.0_p_double * log(2.0_p_double))
      else if ( per_w0(2) > 0.0_p_double ) then
        this%per_w0(2) = per_w0(2)
      endif

      if ( per_focus(2) /= -huge(1.0_p_double) ) then
        this%per_focus(2) = per_focus(2)
      else
        this%per_focus(2) = this%per_focus(1)
      endif

      if ((this%per_w0(2) /= this%per_w0(1)) .or. (this%per_focus(2) /= this%per_focus(1))) then
        this%per_type = p_hermite_gaussian_astigmatic
      endif

    endif

  case ("laguerre")

    this%per_type = p_laguerre_gaussian
    if ( p_x_dim > 1 ) then
      this%per_tem_mode = 0
      do i = 1, 2
        this%per_tem_mode(i) = per_tem_mode(i)
        if ( per_tem_mode(i) > 7 ) then
          if (mpi_node() == 0) then
            print *, ''
            print *, '   Error in zpulse_mov_wall parameters'
            print *, '   Only Laguerre modes up to 7 have been implemented. Please contact the'
            print *, '   development team if you need higher order modes.'
            print *, '   aborting...'
          endif
          stop
        endif
      enddo
    endif


    ! get main spot size and focal plane
    this%per_w0 = 0.0_p_double
    if (per_fwhm(1) > 0.0_p_double) then
      this%per_w0 = per_fwhm(1)/sqrt(4.0_p_double * log(2.0_p_double))
    else
      this%per_w0 = per_w0(1)
    endif

    if ( this%per_w0(1) <= 0.0_p_double ) then
      if (mpi_node() == 0) then
        print *, ''
        print *, '   Error in zpulse_mov_wall parameters'
        print *, '   for Laguerre-Gaussian pulses per_w0 or per_fwhm must be > 0'
        print *, '   aborting...'
      endif
      stop
    endif

    if ( per_focus(1) == -huge(1.0_p_double) ) then
      if (mpi_node() == 0) then
        print *, ''
        print *, '   Error in zpulse_mov_wall parameters'
        print *, '   for Laguerre-Gaussian pulses per_focus must be defined'
        print *, '   aborting...'
      endif
      stop
    endif
    this%per_focus = per_focus(1)

    ! default parameters
    if (per_fwhm(1) > 0.0_p_double) then
      this%per_w0(1) = per_fwhm(1)/sqrt(4.0_p_double * log(2.0_p_double))
    else
      this%per_w0(1) = per_w0(1)
    endif
    this%per_focus(1) = per_focus(1)

  case default
    if (mpi_node() == 0) then
      write(0,*) ''
      write(0,*) '   Error in zpulse_mov_wall parameters'
      write(0,*) '   per_type must be one of the following:'
      write(0,*) '   "plane", "gaussian"/"hermite"'
      write(0,*) '   aborting...'
    endif
    stop
  end select

  ! Set default per_center values
  if ( p_x_dim > 1 ) then
    if ( per_center(1) == -huge(1.0_p_double) ) then
        per_center(1) = 0.5_p_double * ( xmin( g_space, 2 ) + xmax( g_space, 2 ) )
    endif

    if ( p_x_dim == 3 ) then
      if ( per_center(2) == -huge(1.0_p_double) ) then
          per_center(2) = 0.5_p_double * ( xmin( g_space, 3 ) + xmax( g_space, 3 ) )
      endif
    endif
  endif

  this%per_center = per_center

  ! if ( wall_vel < 0.0 .or. wall_vel > 1.0 ) then
  ! if ( mpi_node() == 0 ) then
  ! write(0,*) ""
  ! write(0,*) "   Error in zpulse_mov_wall parameters"
  ! write(0,*) "   When using a moving wall, the velocity (wall_vel) must be 0 <= v <= 1"
  ! write(0,*) "   aborting..."
  ! endif
  ! stop
  ! endif

  if ( wall_pos == -huge(1.0d0) ) then
    if ( mpi_node() == 0 ) then
      write(0,*) ""
      write(0,*) "   Error in zpulse_mov_wall parameters"
      write(0,*) "   When using a moving wall, the initial position (wall_pos) must be set"
      write(0,*) "   aborting..."
    endif
    stop
  endif

  this%if_launch_t0 = if_launch_t0

  this%wall_pos = wall_pos
  this%wall_vel = wall_vel
  this%ncells = ncells
  this%if_clean_fields = if_clean_fields
  this%clean_fields_delta = clean_fields_delta
  this%if_launch_from_wall = if_launch_from_wall
  this%if_apply_divcorr = if_apply_divcorr
  this%if_add_fields = if_add_fields

  this%launch_time = launch_time
  this%if_launch = if_launch

  ! Boost pulse
  this%if_boost = .false.

  if (sim_options%gamma > 1.0_p_double) then

    this % gamma = sim_options%gamma
    this % gamma_times_beta = sqrt(this%gamma**2-1.0_p_double)

    if( this % propagation .eq. p_forward ) then
      this % boost_per_fac = this % gamma - this % gamma_times_beta
    else
      this % boost_per_fac = this % gamma + this % gamma_times_beta
    endif
    if (mpi_node()==0) print *,'   boosting zpulse with gamma = ', this%gamma

    this%if_boost = .true.
  else
    this % boost_per_fac = 1.0_p_double
  endif

  this % if_clean_fields = if_clean_fields

  this % clean_fields_time = this % t_duration() + this % clean_fields_delta

end subroutine
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
function t_env_zpulse_mov_wall( this, t_sim, z_sim ) result(t_env)

  implicit none

  class( t_zpulse_mov_wall ), intent(inout) :: this ! intent must be inout because of
                                                ! the eval( f_parser ) function
  real(p_double), intent(in) :: t_sim
  real(p_double), intent(in) :: z_sim

  ! lab coordinates
  real(p_double) :: z, t

  real(p_double) :: t_env

  real(p_k_fparse), dimension(1) ::t_dbl

  ! convert to lab frame coords
  if( this % if_boost ) then
    call this % apply_lorentz_transformation( z_sim, t_sim, z, t )
  else
    z = z_sim
    t = t_sim
  endif

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

end function t_env_zpulse_mov_wall
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
! Total duration of the laser pulse
!-----------------------------------------------------------------------------------------
function t_duration_zpulse_mov_wall( this ) result( t_duration )

  implicit none

  class( t_zpulse_mov_wall ), intent(in) :: this
  real(p_double) :: t_duration

  select case ( this%tenv_type )
  case (p_polynomial, p_sin2)
    t_duration = this%tenv_rise + this%tenv_flat + this%tenv_fall

  case default
    t_duration = this%tenv_launch_duration

  end select

  if( this%propagation .eq. p_forward ) then
    t_duration = t_duration - this % xi0
  else
    t_duration = t_duration + this % xi0
  endif

  if( this % if_boost ) then
    t_duration = this % gamma * t_duration! - this % gamma_times_beta * this % wall_pos
  endif

end function t_duration_zpulse_mov_wall
!----------------------------------------------------------------------------------------



!----------------------------------------------------------------------------------------
function lon_center_zpulse_mov_wall( this ) result( lon_center )

  implicit none

  class( t_zpulse_mov_wall ), intent(in) :: this
  real(p_double) :: lon_center

  select case ( this%tenv_type )

  case (p_polynomial, p_sin2)
    ! the center of the pulse is define as being
    ! half way in the flat region
    lon_center = this % tenv_rise + this % tenv_flat / 2

  case (p_gaussian, p_func)
    lon_center = 0.5*this % tenv_range

  case default
    lon_center = 0

  end select

  lon_center = lon_center

end function lon_center_zpulse_mov_wall
!----------------------------------------------------------------------------------------


! to do: fix the below code to use the algorithms in the higher
! dimensional versions
! !-----------------------------------------------------------------------------------------
! subroutine launch_wall_1d_zpulse_mov_wall( this, b, e, g_space, nx_p_min, t, dt, no_co )

! implicit none

! integer, parameter :: rank = 1

! class( t_zpulse_mov_wall ), intent(inout) :: this
! type( t_vdf ), intent(inout) :: b, e
! type( t_space ), intent(in) :: g_space
! integer, dimension(:), intent(in) :: nx_p_min
! real(p_double), intent(in) :: t, dt
! class( t_node_conf ), intent(in) :: no_co


! ! local variables
! real(p_double), dimension(2,rank) :: g_x_range
! real(p_double), dimension(2), parameter :: t_range = 0
! real(p_double), dimension(rank) :: ldx

! real(p_double) :: zmin

! real(p_double) :: cos_pol, sin_pol
! real(p_double) :: amp, pha, gpos, old_gpos

! integer :: gipos, ipos, i1, di, old_gipos
! real(p_double) :: t_0, z

! ! Turn off antenna if pulse has ended
! if ( t > this%t_duration() ) then
! this%if_launch = .false.
! return
! endif

! ! Global box sizes
! call get_x_bnd( g_space, g_x_range )
! ldx(1) = real( b%dx_(1), p_double )
! zmin = real( g_x_range(p_lower,1), p_double )

! cos_pol = cos( this%pol )
! sin_pol = sin( this%pol )
! amp = this%omega0 * this%a0

! ! Initial wall position
! old_gpos = this%wall_pos
! old_gipos = int((old_gpos - zmin)/ldx(1)) + 1
! if ( old_gipos < 1 ) old_gipos = 1

! ! Current wall position
! gpos = this%wall_pos + this%wall_vel*t
! gipos = int((gpos - zmin)/ldx(1)) + 1
! if ( gipos < 1 ) gipos = 1

! ! Number of cells moved
! di = old_gipos - gipos

! if ( gipos >= nx_p_min( 1 ) .and. gipos < nx_p_min( 1 ) + b%nx_(1) - 1 ) then

! ipos = gipos - nx_p_min(1) + 1

! do i1 = ipos, ipos+2
! ! t_0 = t + (di - (i1-ipos))*ldx(1)
! z = gipos + ldx(1) * ( i1 - ipos )

! pha = amp * this % t_envelope( t, z ) * cos( this%omega0*(z-t) + this%phase0 )

! e%f1(1, i1) = 0
! e%f1(2, i1) = real( pha * cos_pol, p_k_fld )
! e%f1(3, i1) = real( pha * sin_pol, p_k_fld )

! enddo

! do i1 = ipos, ipos+1
! b%f1(1,i1) = 0
! b%f1(2,i1) = - 0.5_p_k_fld * (e%f1(3,i1) + e%f1(3,i1+1))
! b%f1(3,i1) = + 0.5_p_k_fld * (e%f1(2,i1) + e%f1(2,i1+1))
! enddo


! do i1 = 0, ipos-1
! b%f1(:,i1) = 0
! e%f1(:,i1) = 0
! enddo

! endif

! end subroutine launch_wall_1d_zpulse_mov_wall
! !-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine launch_wall_2d_zpulse_mov_wall( this, b, e, g_space, nx_p_min, t_sim, dt, no_co )

  implicit none

  ! dummy variables

  integer, parameter :: rank = 2

  class( t_zpulse_mov_wall ), intent(inout) :: this
  type( t_vdf ), intent(inout) :: b, e
  type( t_space ), intent(in) :: g_space
  integer, dimension(:), intent(in) :: nx_p_min
  real(p_double), intent(in) :: t_sim, dt
  class( t_node_conf ), intent(in) :: no_co


  ! local variables

  real(p_double), dimension(2,rank) :: g_x_range
  real(p_double), dimension(rank) :: ldx

  real(p_double) :: zmin, zmax, rmin, r, r_2, z_2

  real(p_double) :: cos_pol, sin_pol, prop_sign
  real(p_double) :: amp, lenv, lenv_2, gpos

  integer :: gipos, ipos, i1, i2
  real(p_double) :: z

  integer :: lbnd, rbnd, i1start, i1finish
  logical :: lbnd_in_domain, rbnd_in_domain, domain_enclosed, lbnd_right_of_domain

  real(p_double) :: phase0_temp, pol_temp

  ! Turn off antenna if pulse has ended
  if ( t_sim > this % t_duration() ) then
    if (mpi_node()==0) print *,"duration exceeded. "
    if (mpi_node()==0) print *,"this%t_duration(): ", this%t_duration()
    this%if_launch = .false.
    return
  ! else
    ! print *, "launching mov wall"
  endif

  if ( this%propagation == p_forward ) then
    prop_sign = 1.0_p_double
  else
    prop_sign = -1.0_p_double
  endif

  ! Global box sizes
  call get_x_bnd( g_space, g_x_range )

  ldx(1:2) = real( b%dx_(1:2), p_double )
  zmin = real( g_x_range(p_lower,1), p_double )
  zmax = real( g_x_range(p_upper,1), p_double )
  rmin = real( g_x_range(p_lower,2), p_double )

  select case( this%interpolation )
  case( p_linear, p_cubic )
    ! Nothing to change here
  case( p_quadratic, p_quartic )
    zmin = zmin + ldx(1)*0.5_p_double
    rmin = rmin + ldx(2)*0.5_p_double
  end select

  cos_pol = cos( this%pol )
  sin_pol = sin( this%pol )
  amp = this%omega0 * this%a0

  if (this%pol_type /= 0) then
    amp = amp / sqrt(2.0_p_double)
  endif

  ! start near domain boundary
  if( this%if_launch_from_wall ) then

    if( this%propagation == p_forward ) then
      gipos = 1
    else
      gipos = floor( ( zmax - zmin ) / ldx(1) ) - 1
    endif

  ! else compute the wall position in global coordinates, stored as gpos
  else
    gpos = this%wall_pos + this%wall_vel*t_sim
    gipos = int((gpos - zmin)/ldx(1)) + 1
  endif

  ipos = gipos - nx_p_min(1) + 1

  ! left and right index bounds of the calculation in local coordinates
  if( this%propagation == p_forward ) then
    lbnd = ipos
    rbnd = ipos + this%ncells
  else
    lbnd = ipos - this%ncells
    rbnd = ipos
  endif

  lbnd_right_of_domain = ( lbnd >= e%nx_(1) )
  lbnd_in_domain = (lbnd >= 1) .and. (lbnd <= e%nx_(1))
  rbnd_in_domain = (rbnd >= 1) .and. (rbnd <= e%nx_(1))
  domain_enclosed = (lbnd < 1) .and. (rbnd > e%nx_(1))

  ! if (mpi_node()==0) print *,gipos, gipos + this%ncells
  ! print *, e%nx_(1), nx_p_min(1), lbnd, rbnd, e%lbound(2), e%ubound(2), lbnd_in_domain, rbnd_in_domain, domain_enclosed

  ! coordinates for local field setting in this domain
  lbnd = max( lbnd, 1 )
  rbnd = min( rbnd, e%nx_(1) )

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

        z = zmin + real( i1 - 1 + nx_p_min(1), p_double ) * ldx(1)
        z_2 = z + 0.5_p_double * ldx(1)

        lenv = amp * this % t_envelope( t_sim, z )
        lenv_2 = amp * this % t_envelope( t_sim, z_2 )

        ! e%f2(1, i1, i2) = 0.0_p_double
        e%f2(2, i1, i2) = e%f2(2, i1, i2) &
                + real( this%boost_per_fac * lenv * this % per_envelope_2d( z, r_2, t_sim ) * cos_pol, p_k_fld )

        e%f2(3, i1, i2) = e%f2(3, i1, i2) &
                + real( this%boost_per_fac * lenv * this % per_envelope_2d( z, r, t_sim ) * sin_pol, p_k_fld )

        ! b%f2(1, i1, i2) = 0.0_p_double
        b%f2(2, i1, i2) = b%f2(2, i1, i2) &
                + real( - this%boost_per_fac * lenv_2 * this % per_envelope_2d( z_2, r, t_sim ) * sin_pol * prop_sign, p_k_fld )

        b%f2(3, i1, i2) = b%f2(3, i1, i2) &
                + real( this%boost_per_fac * lenv_2 * this % per_envelope_2d( z_2, r_2, t_sim ) * cos_pol * prop_sign, p_k_fld )

      enddo
    enddo

    if (this%pol_type /= 0) then

      phase0_temp = this%phase0
      pol_temp = this%pol

      this%phase0 = this%phase0 + real( sign( pi_2, real( this%pol_type, p_double )), p_double )
      this%pol = this%pol + real( pi_2, p_double )

      cos_pol = cos( this%pol )
      sin_pol = sin( this%pol )

      do i2 = lbound(e%f2, 3), ubound(e%f2, 3)

        r = rmin + (i2 - 1 + nx_p_min(2)) * ldx(2)
        r_2 = r + 0.5_p_double * ldx(2)

        do i1 = i1start, i1finish

          z = zmin + real( i1 - 1 + nx_p_min(1), p_double ) * ldx(1)
          lenv = amp * this % t_envelope( t_sim, z )
          lenv_2 = amp * this % t_envelope( t_sim, z_2 )

          e%f2(2, i1, i2) = e%f2(2, i1, i2) + real( this%boost_per_fac &
                  * lenv * this % per_envelope_2d( z, r_2, t_sim ) * cos_pol, p_k_fld )

          e%f2(3, i1, i2) = e%f2(3, i1, i2) + real( this%boost_per_fac &
                  * lenv * this % per_envelope_2d( z, r, t_sim ) * sin_pol, p_k_fld )

          b%f2(2, i1, i2) = b%f2(2, i1, i2) - real( this%boost_per_fac &
                  * lenv_2 * this % per_envelope_2d( z_2, r, t_sim ) * sin_pol * prop_sign, p_k_fld )

          b%f2(3, i1, i2) = b%f2(3, i1, i2) + real( this%boost_per_fac &
                  * lenv_2 * this % per_envelope_2d( z_2, r_2, t_sim ) * cos_pol * prop_sign, p_k_fld )

        enddo
      enddo

      this%phase0 = phase0_temp
      this%pol = pol_temp

    endif

    ! ! set magnetic fields using the zpulse bint method
    ! do i2 = lbound(e%f2, 3), ubound(e%f2, 3)
    ! do i1 = i1start, i1finish - 1

    ! b%f2(1,i1,i2) = 0.0_p_k_fld
    ! b%f2(2,i1,i2) = - 0.5_p_k_fld * (e%f2(3,i1,i2) + e%f2(3,i1+1,i2)) * prop_sign
    ! b%f2(3,i1,i2) = + 0.5_p_k_fld * (e%f2(2,i1,i2) + e%f2(2,i1+1,i2)) * prop_sign

    ! enddo
    ! enddo

  endif


  if( this%if_clean_fields ) then

    if( lbnd_right_of_domain ) then
      e%f2(:,:,:) = 0
      b%f2(:,:,:) = 0
    else if( lbnd_in_domain ) then
      e%f2(:,:i1start-1,:) = 0
      b%f2(:,:i1start-1,:) = 0
    endif

  endif

end subroutine launch_wall_2d_zpulse_mov_wall
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine launch_wall_3d_zpulse_mov_wall( this, b, e, g_space, nx_p_min, t_sim, dt, no_co )

  implicit none

  ! dummy variables

  integer, parameter :: rank = 2

  class( t_zpulse_mov_wall ), intent(inout) :: this
  type( t_vdf ), intent(inout) :: b, e
  type( t_space ), intent(in) :: g_space
  integer, dimension(:), intent(in) :: nx_p_min
  real(p_double), intent(in) :: t_sim, dt
  class( t_node_conf ), intent(in) :: no_co


  ! local variables

  real(p_double), dimension(2,rank) :: g_x_range
  real(p_double), dimension(rank) :: ldx

  real(p_double) :: zmin, zmax, r1min, r2min, r1, r2, r1_2, r2_2, z_2

  real(p_double) :: cos_pol, sin_pol, prop_sign
  real(p_double) :: amp, lenv, lenv_2, gpos

  integer :: gipos, ipos, i1, i2, i3
  real(p_double) :: z

  integer :: lbnd, rbnd, i1start, i1finish
  logical :: lbnd_in_domain, rbnd_in_domain, domain_enclosed, lbnd_right_of_domain

  real(p_double) :: phase0_temp, pol_temp

  ! Turn off antenna if pulse has ended
  if ( t_sim > this % t_duration() ) then
    if (mpi_node()==0) print *,"duration exceeded. "
    if (mpi_node()==0) print *,"this%t_duration(): ", this%t_duration()
    this%if_launch = .false.
    return
  ! else
    ! print *, "launching mov wall"
  endif

  if ( this%propagation == p_forward ) then
    prop_sign = 1.0_p_double
  else
    prop_sign = -1.0_p_double
  endif

  ! Global box sizes
  call get_x_bnd( g_space, g_x_range )

  ldx(1:3) = real( b%dx_(1:3), p_double )
  zmin = real( g_x_range(p_lower,1), p_double )
  zmax = real( g_x_range(p_upper,1), p_double )
  r1min = real( g_x_range(p_lower,2), p_double )
  r2min = real( g_x_range(p_lower,3), p_double )

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
  amp = this%omega0 * this%a0

  if (this%pol_type /= 0) then
    amp = amp / sqrt(2.0_p_double)
  endif

  ! start near domain boundary
  if( this%if_launch_from_wall ) then

    if( this%propagation == p_forward ) then
      gipos = 1
    else
      gipos = floor( ( zmax - zmin ) / ldx(1) ) - 1
    endif

  ! else compute the wall position in global coordinates, stored as gpos
  else
    gpos = this%wall_pos + this%wall_vel*t_sim
    gipos = int((gpos - zmin)/ldx(1)) + 1
  endif

  ipos = gipos - nx_p_min(1) + 1

  ! left and right index bounds of the calculation in local coordinates
  if( this%propagation == p_forward ) then
    lbnd = ipos
    rbnd = ipos + this%ncells
  else
    lbnd = ipos - this%ncells
    rbnd = ipos
  endif

  lbnd_in_domain = (lbnd >= 1) .and. (lbnd <= e%nx_(1))
  rbnd_in_domain = (rbnd >= 1) .and. (rbnd <= e%nx_(1))
  domain_enclosed = (lbnd < 1) .and. (rbnd > e%nx_(1))

  ! coordinates for local field setting in this domain
  lbnd = max( lbnd, 1 )
  rbnd = min( rbnd, e%nx_(1) )

  if( lbnd_in_domain .or. rbnd_in_domain .or. domain_enclosed ) then

    i1start = lbnd
    i1finish = rbnd

    if( .not. this%if_add_fields ) then
      e%f2(:,i1start:i1finish,:) = 0
      b%f2(:,i1start:i1finish,:) = 0
    endif

    do i2 = lbound(e%f3, 3), ubound(e%f3, 3)

      r1 = r1min + (i2 - 1 + nx_p_min(2)) * ldx(2)
      r1_2 = r1 + 0.5_p_double * ldx(2)

      do i3 = lbound(e%f3, 4), ubound(e%f3, 4)

        r2 = r2min + (i2 - 1 + nx_p_min(3)) * ldx(3)
        r2_2 = r2 + 0.5_p_double * ldx(3)

        do i1 = i1start, i1finish

          z = zmin + real( i1 - 1 + nx_p_min(1), p_double ) * ldx(1)
          z_2 = z + 0.5_p_double * ldx(1)

          lenv = amp * this % t_envelope( t_sim, z )
          lenv_2 = amp * this % t_envelope( t_sim, z_2 )

          ! e%f2(1, i1, i2) = 0.0_p_double
          e%f3(2, i1, i2, i3) = e%f3(2, i1, i2, i3) + real( this%boost_per_fac &
                        * lenv * this % per_envelope_3d( z, r1_2, r2, t_sim ) * cos_pol, p_k_fld )

          e%f3(3, i1, i2, i3) = e%f3(3, i1, i2, i3) + real( this%boost_per_fac &
                        * lenv * this % per_envelope_3d( z, r1, r2_2, t_sim ) * sin_pol, p_k_fld )

          ! b%f2(1, i1, i2) = 0.0_p_double
          b%f3(2, i1, i2, i3) = b%f3(2, i1, i2, i3) + real( - this%boost_per_fac &
                        *lenv_2 * this % per_envelope_3d( z_2, r1, r2_2, t_sim ) &
                        * sin_pol * prop_sign, p_k_fld )

          b%f3(3, i1, i2, i3) = b%f3(3, i1, i2, i3) + real( this%boost_per_fac &
                        * lenv_2 * this % per_envelope_3d( z_2, r1_2, r2, t_sim ) &
                        * cos_pol * prop_sign, p_k_fld )

        enddo
      enddo
    enddo

   if (this%pol_type /= 0) then

      phase0_temp = this%phase0
      pol_temp = this%pol

      this%phase0 = this%phase0 + real( sign( pi_2, real( this%pol_type, p_double )), p_double )
      this%pol = this%pol + real( pi_2, p_double )

      cos_pol = cos( this%pol )
      sin_pol = sin( this%pol )

      do i2 = lbound(e%f3, 3), ubound(e%f3, 3)

        r1 = r1min + (i2 - 1 + nx_p_min(2)) * ldx(2)
        r1_2 = r1 + 0.5_p_double * ldx(2)

        do i3 = lbound(e%f3, 3), ubound(e%f3, 3)

          r2 = r2min + (i2 - 1 + nx_p_min(3)) * ldx(3)
          r2_2 = r2 + 0.5_p_double * ldx(3)

          do i1 = i1start, i1finish

            z = zmin + real( i1 - 1 + nx_p_min(1), p_double ) * ldx(1)
            lenv = amp * this % t_envelope( t_sim, z )
            lenv_2 = amp * this % t_envelope( t_sim, z_2 )

            e%f3(2, i1, i2, i3) = e%f3(2, i1, i2, i3) + real( this%boost_per_fac &
                    * lenv * this % per_envelope_3d( z, r1_2, r2, t_sim ) * cos_pol, p_k_fld )

            e%f3(3, i1, i2, i3) = e%f3(3, i1, i2, i3) + real( this%boost_per_fac &
                    * lenv * this % per_envelope_3d( z, r1, r2_2, t_sim ) * sin_pol, p_k_fld )

            b%f3(2, i1, i2, i3) = b%f3(2, i1, i2, i3) - real( this%boost_per_fac &
                    * lenv_2 * this % per_envelope_3d( z_2, r1, r2_2, t_sim ) &
                    * sin_pol * prop_sign, p_k_fld )

            b%f3(3, i1, i2, i3) = b%f3(3, i1, i2, i3) + real( this%boost_per_fac &
                    * lenv_2 * this % per_envelope_3d( z_2, r1_2, r2, t_sim ) &
                    * cos_pol * prop_sign, p_k_fld )

          enddo
        enddo
      enddo

      this%phase0 = phase0_temp
      this%pol = pol_temp

    endif

  endif

end subroutine launch_wall_3d_zpulse_mov_wall
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine apply_divcorr_2d_zpulse_mov_wall( this, b, e, g_space, nx_p_min, t, dt, no_co )

  implicit none
  integer, parameter :: rank = 2

  class( t_zpulse_mov_wall ), intent(inout) :: this
  type( t_vdf ), intent(inout) :: b, e
  type( t_space ), intent(in) :: g_space
  integer, dimension(:), intent(in) :: nx_p_min
  real(p_double), intent(in) :: t, dt
  class( t_node_conf ), intent(in) :: no_co

  real(p_double) :: e2p, e2m, e1, b2p, b2m, b1
  integer :: i1, i2

  integer :: lbnd, rbnd
  logical :: lbnd_in_domain, rbnd_in_domain, domain_enclosed

  integer :: i1start, i1finish, direction, i1_initial_data
  integer :: recv_neighbor, send_neighbor
  logical :: start_in_domain, finish_in_domain

  real(p_double), dimension(2,rank) :: g_x_range
  real(p_double), dimension(rank) :: ldx
  real(p_double) :: zmin, zmax, rmin, z_2
  real(p_double) :: gpos, gipos, ipos

  integer :: size
  integer :: ierr
  integer, dimension( MPI_STATUS_SIZE ) :: status
  real(p_double), dimension(:,:), allocatable :: divcorr_buf

  ! start and finish are the start and finish values of the calculation
  ! for forward pulse, start at the end of the wall in the forward direction
  ! and finish at the start of the pulse
  ! the field outside the fields we just set (up to ipos + ncells ) is assumed to be correct


  ! ! Turn off antenna if pulse has ended
  ! if ( t > this % t_duration() ) then
  ! print *, "duration exceeded. "
  ! print *, "this%t_duration(): ", this%t_duration()
  ! print *, "t: ", t
  ! this%if_launch = .false.
  ! return
  ! ! else
  ! ! print *, "launching mov wall"
  ! endif

  ! Global box sizes
  call get_x_bnd( g_space, g_x_range )

  ldx(1:2) = real( b%dx_(1:2), p_double )
  zmin = real( g_x_range(p_lower,1), p_double )
  zmax = real( g_x_range(p_upper,1), p_double )
  rmin = real( g_x_range(p_lower,2), p_double )

  select case( this%interpolation )
  case( p_linear, p_cubic )
    ! Nothing to change here
  case( p_quadratic, p_quartic )
    zmin = zmin + ldx(1)*0.5_p_double
    rmin = rmin + ldx(2)*0.5_p_double
  end select

  ! left and right bounds of the calculation
  if( this%if_launch_t0 ) then

      lbnd = 1 - nx_p_min(1)
      rbnd = int( (zmax - zmin)/ldx(1) ) - nx_p_min(1)

  else

    ! Current wall position
    gpos = this%wall_pos + this%wall_vel*t
    gipos = int((gpos - zmin)/ldx(1)) + 1
    ipos = gipos - + 1

    if( (gpos - ldx(1) < zmin) .or. (gpos + (this%ncells + 1) * ldx(1) > zmax ) ) then
      if (mpi_node()==0) print *,"WARNING: zpulse_mov_wall hit the boundary before pulse injection was complete. Stopping injection."
      this%if_launch = .false.
      return
    endif

    if( this%propagation == p_forward ) then
      lbnd = ipos
      rbnd = ipos + this%ncells + 1
    else
      lbnd = ipos - this%ncells - 1
      rbnd = ipos
    endif

  endif

  lbnd_in_domain = ( lbnd >= 1) .and. (lbnd <= e%nx_(1))
  rbnd_in_domain = (rbnd >= 1) .and. (rbnd <= e%nx_(1))
  domain_enclosed = (lbnd < 1) .and. (rbnd > e%nx_(1))

  ! wall outside is completely outside this domain: do nothing
  if( .not. this%if_launch_t0 .and. .not.( lbnd_in_domain .or. rbnd_in_domain .or. domain_enclosed ) ) return

  if( this%propagation == p_forward ) then
    direction = -1
    recv_neighbor = p_upper
    send_neighbor = p_lower
    i1start = min( rbnd - 1, e%nx_(1) )
    i1finish = max( lbnd, 1 )
    i1_initial_data = rbnd
    start_in_domain = rbnd_in_domain
    finish_in_domain = lbnd_in_domain
  else
    direction = 1
    recv_neighbor = p_lower
    send_neighbor = p_upper
    i1start = max( 1, lbnd + 1 )
    i1finish = min( rbnd, e%nx_(1) )
    i1_initial_data = lbnd
    start_in_domain = lbnd_in_domain
    finish_in_domain = rbnd_in_domain
  endif

  ! print *, i1start, i1finish, e%nx_(1), lbnd_in_domain, rbnd_in_domain, domain_enclosed

  size = e%nx_(2)
  allocate( divcorr_buf( 2, size ) )

  ! check if integral start is in domain; else we need calculated values
  ! from nearest neighbor
  if( .not. start_in_domain ) then
    ! print *, "calling recv."
    call mpi_recv( divcorr_buf, 2 * size, MPI_DOUBLE, &
                   no_co%neighbor( recv_neighbor, 1 ) - 1, &
                   0, no_co%comm, status, ierr )

    print *, "recv"

    call check_error(ierr,"zpulse wall: recv fail",p_err_mpi,"zpulse/os-zpulse-mov-wall.f03",1487)

  endif

  do i2 = 1, e%nx_(2)

    ! starting value: either use value from previous node if the
    ! divergence integral is not starting locally
    if( start_in_domain ) then
      e1 = e%f2( 1, i1_initial_data, i2 )
      b1 = b%f2( 1, i1_initial_data, i2 )
    else
      e1 = divcorr_buf(1,i2)
      b1 = divcorr_buf(2,i2)
    endif

    do i1 = i1start, i1finish, direction

      ! old method:
      ! e2p = e%f2(2,i1+1,i2)
      ! e2m = e%f2(2,i1+1,i2-1)

      e2p = e%f2(2,i1,i2)
      e2m = e%f2(2,i1,i2-1)

      e1 = e1 - direction * ldx(1)*(e2p - e2m)/ldx(2)
      e%f2(1, i1, i2) = real( e1, p_k_fld )

      b2p = b%f2(2,i1,i2+1)
      b2m = b%f2(2,i1,i2)

      b1 = b1 - direction * ldx(1) * (b2p - b2m) / ldx(2)
      b%f2(1, i1, i2) = real(b1, p_k_fld )

    enddo ! i1

  enddo ! i2

  ! stop

  ! if integral finish is not in domain, we need to send the results of the integral
  ! to be finished in the next domain
  if( .not. finish_in_domain ) then

    divcorr_buf(1,1:size) = e%f2(1, i1finish, 1:size)
    divcorr_buf(2,1:size) = b%f2(1, i1finish, 1:size)

    call mpi_send( divcorr_buf, 2 * size, MPI_DOUBLE, &
                   no_co%neighbor( send_neighbor, 1 ) - 1, &
                   0, no_co%comm, ierr )

    call check_error(ierr,"zpulse wall: send fail",p_err_mpi,"zpulse/os-zpulse-mov-wall.f03",1538)

  endif

 deallocate( divcorr_buf )

end subroutine apply_divcorr_2d_zpulse_mov_wall
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
subroutine apply_divcorr_3d_zpulse_mov_wall( this, b, e, g_space, nx_p_min, t, dt, no_co )

  implicit none
  integer, parameter :: rank = 2

  class( t_zpulse_mov_wall ), intent(inout) :: this
  type( t_vdf ), intent(inout) :: b, e
  type( t_space ), intent(in) :: g_space
  integer, dimension(:), intent(in) :: nx_p_min
  real(p_double), intent(in) :: t, dt
  class( t_node_conf ), intent(in) :: no_co

 ! real(p_double) :: e2p, e2m, e1, b2p, b2m, b1
 ! integer :: i1, i2

 ! integer :: lbnd, rbnd
 ! logical :: lbnd_in_domain, rbnd_in_domain, domain_enclosed

 ! integer :: i1start, i1finish, direction, i1_initial_data
 ! integer :: recv_neighbor, send_neighbor
 ! logical :: start_in_domain, finish_in_domain

 ! real(p_double), dimension(2,rank) :: g_x_range
 ! real(p_double), dimension(rank) :: ldx
 ! real(p_double) :: zmin, zmax, rmin
 ! real(p_double) :: gpos, gipos, ipos

 ! integer :: size
 ! integer :: ierr
 ! integer, dimension( MPI_STATUS_SIZE ) :: status
 ! real(p_double), dimension(:,:), allocatable :: divcorr_buf

 ! ! start and finish are the start and finish values of the calculation
 ! ! for forward pulse, start at the end of the wall in the forward direction
 ! ! and finish at the start of the pulse
 ! ! the field outside the fields we just set (up to ipos + ncells ) is assumed to be correct

 ! ! Global box sizes
 ! call get_x_bnd( g_space, g_x_range )

 ! ldx(1:2) = real( b%dx_(1:2), p_double )
 ! zmin = real( g_x_range(p_lower,1), p_double )
 ! zmax = real( g_x_range(p_upper,1), p_double )
 ! rmin = real( g_x_range(p_lower,2), p_double )

 ! select case( this%interpolation )
 ! case( p_linear, p_cubic )
 ! ! Nothing to change here
 ! case( p_quadratic, p_quartic )
 ! zmin = zmin + ldx(1)*0.5_p_double
 ! rmin = rmin + ldx(2)*0.5_p_double
 ! end select

 ! ! start near domain boundary
 ! if( this%if_launch_from_wall ) then

 ! if( this%propagation == p_forward ) then
 ! gipos = 1
 ! else
 ! gipos = floor( ( zmax - zmin ) / ldx(1) ) - 1
 ! endif

 ! ! else compute the wall position in global coordinates, stored as gpos
 ! else
 ! gpos = this%wall_pos + this%wall_vel*t
 ! gipos = int((gpos - zmin)/ldx(1)) + 1

 ! if( (gpos - ldx(1) < zmin) .or. (gpos + (this%ncells + 1) * ldx(1) > zmax ) ) then
 ! if (mpi_node()==0) print *,"WARNING: zpulse_mov_wall hit the boundary before pulse injection was complete. Stopping injection."
 ! this%if_launch = .false.
 ! return
 ! endif

 ! endif

 ! ipos = gipos - nx_p_min(1) + 1

 ! ! left and right bounds of the calculation
 ! if( this%propagation == p_forward ) then
 ! lbnd = ipos
 ! rbnd = ipos + this%ncells + 1
 ! else
 ! lbnd = ipos - this%ncells - 1
 ! rbnd = ipos
 ! endif

 ! lbnd_in_domain = ( lbnd >= 1) .and. (lbnd <= e%nx_(1))
 ! rbnd_in_domain = (rbnd >= 1) .and. (rbnd <= e%nx_(1))
 ! domain_enclosed = (lbnd < 1) .and. (rbnd > e%nx_(1))


 ! ! wall outside is completely outside this domain: do nothing
 ! if( .not.( lbnd_in_domain .or. rbnd_in_domain .or. domain_enclosed ) ) return

 ! if( this%propagation == p_forward ) then
 ! direction = -1
 ! recv_neighbor = p_upper
 ! send_neighbor = p_lower
 ! i1start = min( rbnd - 1, e%nx_(1) )
 ! i1finish = max( lbnd, 1 )
 ! i1_initial_data = rbnd
 ! start_in_domain = rbnd_in_domain
 ! finish_in_domain = lbnd_in_domain
 ! else
 ! direction = 1
 ! recv_neighbor = p_lower
 ! send_neighbor = p_upper
 ! i1start = max( 1, lbnd + 1 )
 ! i1finish = min( rbnd, e%nx_(1) )
 ! i1_initial_data = lbnd
 ! start_in_domain = lbnd_in_domain
 ! finish_in_domain = rbnd_in_domain
 ! endif

  ! print *, i1start, i1finish, e%nx_(1), lbnd_in_domain, rbnd_in_domain, domain_enclosed

  ! size = e%nx_(2)
  ! allocate( divcorr_buf( 2, size ) )

  ! ! check if integral start is in domain; else we need calculated values
  ! ! from nearest neighbor
  ! if( .not. start_in_domain ) then
  ! ! print *, "calling recv."
  ! call mpi_recv( divcorr_buf, 2 * size, MPI_DOUBLE, &
  ! no_co%neighbor( recv_neighbor, 1 ) - 1, &
  ! 0, no_co%comm, status, ierr )

  ! call check_error(ierr,"zpulse wall: recv fail",p_err_mpi,"zpulse/os-zpulse-mov-wall.f03",1677)

  ! endif

  ! do i3 = 1, e%nx_(3)

  ! do i2 = 1, e%nx_(2)

  ! ! starting value: either use value from previous node if the
  ! ! divergence integral is not starting locally
  ! if( start_in_domain ) then
  ! e1 = e%f3( 1, i1_initial_data, i2 )
  ! b1 = b%f3( 1, i1_initial_data, i2 )
  ! else
  ! e1 = divcorr_buf(1,i2,i3)
  ! b1 = divcorr_buf(2,i2,i3)
  ! endif

  ! do i1 = i1start, i1finish, direction

  ! ! old method:
  ! ! e2p = e%f2(2,i1+1,i2)
  ! ! e2m = e%f2(2,i1+1,i2-1)

  ! e2p = e%f2(2,i1,i2)
  ! e2m = e%f2(2,i1,i2-1)

  ! e1 = e1 - direction * ldx(1)*(e2p - e2m)/ldx(2)
  ! e%f2(1, i1, i2) = real( e1, p_k_fld )

  ! b2p = b%f2(2,i1,i2+1)
  ! b2m = b%f2(2,i1,i2)

  ! b1 = b1 - direction * ldx(1) * (b2p - b2m) / ldx(2)
  ! b%f2(1, i1, i2) = real(b1, p_k_fld )

  ! enddo ! i1

  ! enddo ! i2
  ! enddo ! i3

  ! ! stop

  ! ! if integral finish is not in domain, we need to send the results of the integral
  ! ! to be finished in the next domain
  ! if( .not. finish_in_domain ) then

  ! divcorr_buf(1,:) = e%f2(1, i1finish, 1:size)
  ! divcorr_buf(2,:) = b%f2(1, i1finish, 1:size)

  ! call mpi_send( divcorr_buf, 2 * size, MPI_DOUBLE, &
  ! no_co%neighbor( send_neighbor, 1 ) - 1, &
  ! 0, no_co%comm, ierr )

  ! call check_error(ierr,"zpulse wall: send fail",p_err_mpi,"zpulse/os-zpulse-mov-wall.f03",1731)

  ! endif

 ! deallocate( divcorr_buf )

end subroutine apply_divcorr_3d_zpulse_mov_wall
!-----------------------------------------------------------------------------------------



! !-----------------------------------------------------------------------------------------
! subroutine clean_fields_2d_zpulse_mov_wall( this, b, e, g_space, nx_p_min, t, dt, no_co )

! implicit none
! integer, parameter :: rank = 2

! class( t_zpulse_mov_wall ), intent(inout) :: this
! type( t_vdf ), intent(inout) :: b, e
! type( t_space ), intent(in) :: g_space
! integer, dimension(:), intent(in) :: nx_p_min
! real(p_double), intent(in) :: t, dt
! class( t_node_conf ), intent(in) :: no_co

! call this % get_wall_vars_2d( g_space, lbnd, rbnd )

! lbnd_in_domain = (lbnd >= 1) .and. (lbnd <= e%nx_(1))
! rbnd_in_domain = (rbnd >= 1) .and. (rbnd <= e%nx_(1))
! domain_enclosed = (lbnd < 1) .and. (rbnd > e%nx_(1))

! ! coordinates for local field setting in this domain
! lbnd = max( lbnd, 1 )
! rbnd = min( rbnd, e%nx_(1) )


! end subroutine clean_fields_2d_zpulse_mov_wall
! !-----------------------------------------------------------------------------------------


! !-----------------------------------------------------------------------------------------
! subroutine get_wall_vars_2d_zpulse_mov_wall( this, g_space, nx_p_min, t, &
! zmin, rmin )

! implicit none
! integer, parameter :: rank = 2

! class( t_zpulse_mov_wall ), intent(inout) :: this
! type( t_vdf ), intent(inout) :: b, e
! type( t_space ), intent(in) :: g_space
! integer, dimension(:), intent(in) :: nx_p_min
! real(p_double), intent(in) :: t

! ! Global box sizes
! call get_x_bnd( g_space, g_x_range )

! ldx(1:2) = real( b%dx_(1:2), p_double )
! zmin = real( g_x_range(p_lower,1), p_double )
! zmax = real( g_x_range(p_upper,1), p_double )
! rmin = real( g_x_range(p_lower,2), p_double )

! select case( this%interpolation )
! case( p_linear, p_cubic )
! ! Nothing to change here
! case( p_quadratic, p_quartic )
! zmin = zmin + ldx(1)*0.5_p_double
! rmin = rmin + ldx(2)*0.5_p_double
! end select

! ! start near domain boundary
! if( this%if_launch_from_wall ) then

! if( this%propagation == p_forward ) then
! gipos = 1
! else
! gipos = floor( ( zmax - zmin ) / ldx(1) ) - 1
! endif

! ! else compute the wall position in global coordinates, stored as gpos
! else
! gpos = this%wall_pos + this%wall_vel*t
! gipos = int((gpos - zmin)/ldx(1)) + 1
! endif

! ipos = gipos - nx_p_min(1) + 1

! ! left and right index bounds of the calculation in local coordinates
! if( this%propagation == p_forward ) then
! lbnd = ipos
! rbnd = ipos + this%ncells
! else
! lbnd = ipos - this%ncells
! rbnd = ipos
! endif

! end subroutine get_wall_vars_2d_zpulse_mov_wall
! !-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
subroutine launch_zpulse_mov_wall( this, emf, g_space, &
                                   nx_p_min, g_nx, t, dt, no_co )

  implicit none

  class( t_zpulse_mov_wall ), intent(inout) :: this
  class( t_emf ) , intent(inout) :: emf

  type( t_space ), intent(in) :: g_space
  integer, intent(in), dimension(:) :: nx_p_min, g_nx
  real(p_double), intent(in) :: t
  real(p_double), intent(in) :: dt
  class( t_node_conf ), intent(in) :: no_co

  real(p_double) :: phase0_temp, pol_temp

  if ( this%if_launch .and. ( t > this%launch_time + this % t_duration() ) ) then
    this%if_launch = .false.
  endif

  if ( this%if_launch .and. ( t >= this%launch_time ) ) then

    select case (emf%b%x_dim_)
    ! case(1)
    ! call this % launch_wall_1d( emf%b, emf%e, g_space, nx_p_min, t, dt, no_co )

    case(2)
      call this % launch_wall_2d( emf%b, emf%e, g_space, nx_p_min, t, dt, no_co )

      if( this%if_apply_divcorr ) then
        call this % apply_divcorr_2d( emf%b, emf%e, g_space, nx_p_min, t, dt, no_co )
      endif

    case(3)
      call this % launch_wall_3d( emf%b, emf%e, g_space, nx_p_min, t, dt, no_co )

      if( this%if_apply_divcorr ) then
        call this % apply_divcorr_3d( emf%b, emf%e, g_space, nx_p_min, t, dt, no_co )
      endif

    case default
      write(err_buf__,*) 'Not implemented yet';call err__("zpulse/os-zpulse-mov-wall.f03",1871)
      call abort_program( p_err_notimplemented )

    end select



  endif

  ! if ( this%if_clean_fields .and. ( t >= this%clean_fields_time ) ) then
  ! call this % clean_fields( emf%b, emf%e, g_space, nx_p_min, t, dt )
  ! this % if_clean_fields = .false.
  ! endif

end subroutine launch_zpulse_mov_wall
!-----------------------------------------------------------------------------------------



! !-----------------------------------------------------------------------------------------
! subroutine clean_fields_zpulse_mov_wall( this, b, e, g_space, nx_p_min, t, dt )

! implicit none

! integer, parameter :: rank = 2

! class( t_zpulse_mov_wall ), intent(inout) :: this
! type( t_vdf ), intent(inout) :: b, e
! type( t_space ), intent(in) :: g_space
! integer, dimension(:), intent(in) :: nx_p_min
! real(p_double), intent(in) :: t, dt


! ! local variables
! real(p_double), dimension(2,rank) :: g_x_range
! real(p_double), dimension(rank) :: ldx

! integer :: gipos, ipos, i1
! real(p_double) :: zmin, gpos

! print *, "INFO: clearing fields."

! call get_x_bnd( g_space, g_x_range )

! ldx(1:2) = real( b%dx_(1:2), p_double )
! zmin = real( g_x_range(p_lower,1), p_double )

! select case( this%interpolation )
! case( p_linear, p_cubic )
! ! Nothing to change here
! case( p_quadratic, p_quartic )
! zmin = zmin + ldx(1)*0.5_p_double
! end select

! ! original wall position
! gpos = this%wall_pos ! + this%wall_vel*t
! gipos = int((gpos - zmin)/ldx(1)) + 1
! ipos = gipos - nx_p_min(1) + 1

! do i1 = 1, ipos + this % ncells

! e%f2(1, i1, :) = 0.0_p_double
! e%f2(2, i1, :) = 0.0_p_double
! e%f2(3, i1, :) = 0.0_p_double

! b%f2(1, i1, :) = 0.0_p_double
! b%f2(2, i1, :) = 0.0_p_double
! b%f2(3, i1, :) = 0.0_p_double

! enddo

! end subroutine clean_fields_zpulse_mov_wall
! !-----------------------------------------------------------------------------------------





end module m_zpulse_mov_wall
