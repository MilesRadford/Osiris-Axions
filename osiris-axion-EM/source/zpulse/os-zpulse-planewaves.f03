#include "os-config.h"
#include "os-preprocess.fpp"
 
module m_zpulse_planewaves

#include "memory/memory.h"

use m_zpulse_std

use m_node_conf
use m_restart
use m_parameters
use m_math

use m_vdf_define
use m_vdf_math
use m_vdf 

use m_emf_define, only: t_emf_bound, t_emf

use m_space
use m_grid_define
use m_fparser
use m_utilities

implicit none

private 

! testing, debug and control
logical, parameter :: control = .true.

! longitudinal profiles
integer, parameter :: k_lon_plane 	 = 0
integer, parameter :: k_lon_gaussian = 1
integer, parameter :: k_lon_sin2     = 2

! perpendicular profiles
integer, parameter :: k_per_plane 		      = 0
integer, parameter :: k_per_gaussian          = 1
integer, parameter :: k_per_laguerre_gaussian = 2

! just the imaginaty unity
complex(p_double), parameter :: im_unit = (0.0_p_double, 1.0_p_double)

! parameters for integration ranges
integer, parameter :: n_sig_gaussian = 3
integer, parameter :: n_sig_sin2     = 5

type, extends( t_zpulse ) :: t_zpulse_planewaves

	! new variables go here
	real(p_double), dimension(3) :: focus
	real(p_double), dimension(2) :: angles

	logical                      :: paraxial
	real(p_double)               :: abouraddy_angle
	integer, dimension(3)        :: num_waves
	logical  					 :: normalize
	real(p_double)               :: normalization

	real(p_double)               :: lon_position
	real(p_double), dimension(2) :: lon_int_range

	real(p_double), dimension(2) :: per1_int_range
	real(p_double), dimension(2) :: per2_int_range

	real(p_double) 				 :: tolerance

	logical 					 :: wall_injection ! not yet functioning 

	contains

	procedure :: init     					=> init_zpulse_planewaves
	procedure :: launch     				=> launch_zpulse_planewaves
	procedure :: read_input 				=> read_input_zpulse_planewaves 
	procedure :: check_dimensionality  	 	=> check_dimensionality_zpulse_planewaves

	procedure :: k_function 				=> k_function_zpulse_planewaves
	procedure :: k_lon_function 			=> k_lon_function_zpulse_planewaves
	procedure :: k_per_function 			=> k_per_function_zpulse_planewaves

end type t_zpulse_planewaves

public :: t_zpulse_planewaves

public :: k_per_gaussian, k_per_laguerre_gaussian
public :: im_unit

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Initialize zpulse_planewaves
subroutine init_zpulse_planewaves( this, restart, t , dt, b, g_space, nx_p_min, no_co, interpolation )

	implicit none

	class( t_zpulse_planewaves ), intent(inout) :: this
	logical, intent(in) :: restart
	real( p_double ), intent(in) :: t, dt
	type( t_vdf ) ,  intent(in)  :: b
	type( t_space ), intent(in) :: g_space
	integer, intent(in), dimension(:) :: nx_p_min
	class( t_node_conf ), intent(in) :: no_co
	integer, intent(in) :: interpolation

	logical :: if_launch

	! Need to preserve if_launch value for this zpulse
	if_launch = this%if_launch

	!call this % t_zpulse % init( restart, t, dt, b, g_space, nx_p_min, no_co, interpolation )

	this%if_launch = if_launch

	if ( this%if_launch .and. t > this%launch_time  ) then
		this%if_launch = .false.
	endif

end subroutine init_zpulse_planewaves
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Reading the input file
subroutine read_input_zpulse_planewaves( this, input_file, g_space, bnd_con, periodic, grid, sim_options )

	use m_input_file

	implicit none

	class( t_zpulse_planewaves ), intent(inout) :: this
	class( t_input_file ), intent(inout)     :: input_file
	type( t_space ), intent(in)             :: g_space
	class (t_emf_bound), intent(in)         :: bnd_con
	logical, dimension(:), intent(in)       :: periodic
	class( t_grid ), intent(in)             :: grid
	type( t_options ), intent(in)           :: sim_options

	! ----------------------------------------------------------------------
	! local variables
	real(p_double)               :: a0
	real(p_double)               :: omega0
	integer                      :: pol_type
	real(p_double)               :: pol
	real(p_double), dimension(2) :: angles
	real(p_double), dimension(3) :: focus
	real(p_double)               :: phase
	real(p_double)               :: tolerance
	logical 	               	 :: wall_injection

	character(len=8)             :: preset

	logical                      :: paraxial
	character(len=16)            :: pulse_type
	integer, dimension(3)        :: num_waves
	real(p_double)     			 :: abouraddy_angle

	character(len=16)            :: lon_type
	real(p_double)               :: lon_position
	real(p_double)               :: lon_duration, lon_fwhm, lon_range 			! for gaussian profile
	real(p_double) 			     :: lon_rise, lon_flat, lon_fall 				! for sin2 profile
	real(p_double), dimension(2) :: lon_int_range

	character(len=16)            :: per_type
	real(p_double), dimension(2) :: per_w0, per_fwhm
	integer, dimension(2)        :: per_tem_mode
	real(p_double), dimension(2) :: per1_int_range, per2_int_range

	logical                      :: if_launch
	real(p_double)               :: launch_time

	namelist /nl_zpulse_planewaves/ &
		a0, omega0, pol_type, pol, angles, focus, phase, preset, &	
		lon_type, lon_position, lon_duration, lon_fwhm, lon_range, lon_rise, lon_flat, lon_fall, lon_int_range, &
		per_type, per_w0, per_fwhm, per_tem_mode, per1_int_range, per2_int_range, &
		paraxial, num_waves, pulse_type, abouraddy_angle, tolerance, wall_injection, if_launch

	! ----------------------------------------------------------------------
	! auxiliar variables
	integer 			:: ierr
	real(p_double) 		:: sig ! auxiliar to define integration ranges
	real(p_double)      :: dper, dlon

	! ----------------------------------------------------------------------
	! defining default values
	if_launch          = .true. 
	launch_time        =  0.0_p_double

	a0                 =  1.0_p_double
	omega0             = 6.28_p_double
	pol_type           =  0
	pol                =  0.0_p_double
	angles             = (/ 0.0_p_double, 0.0_p_double /)
	focus              = (/ 0.0_p_double, 0.0_p_double, 0.0_p_double /)
	phase              =  0.0_p_double

	preset             = "0"

	paraxial           = .false.
	pulse_type         = "time"
	num_waves          = (/ 30, 30, 30 /)
	abouraddy_angle    = 0.0_p_double
	tolerance          = 0.0_p_double
	wall_injection     = .false.

	lon_type           = "sin2"
	lon_position       =  0.0_p_double
	lon_duration       = 10.0_p_double
	lon_range          = -huge(1.0_p_double)
	lon_fwhm           =  0.0_p_double
	lon_rise    	   =  pi
	lon_flat           =  0.0_p_double
	lon_fall           =  pi
	lon_int_range      = (/ 0.0_p_double, 0.0_p_double /)

	per_type           = "gaussian"
	per_w0             =  2.0_p_double
	per_fwhm           =  0.0_p_double
	per_tem_mode       = 0 
	per1_int_range     = (/ 0.0_p_double, 0.0_p_double /)
	per2_int_range     = (/ 0.0_p_double, 0.0_p_double /)

	! ----------------------------------------------------------------------
	! Get namelist text from input file

	SCR_ROOT(" - Reading zpulse_planewaves configuration...")

	call get_namelist( input_file, "nl_zpulse_planewaves", ierr )

	if (ierr /= 0) then
	if ( mpi_node() == 0 ) then
		if ( ierr < 0 ) then
		print *, "Error reading zpulse_planewaves parameters"
		else
		print *, "Error: zpulse_planewaves parameters missing"
		endif
		print *, "aborting..."
	endif
	stop
	endif

	! ----------------------------------------------------------------------
	! read input file
	read (input_file%nml_text, nml = nl_zpulse_planewaves, iostat = ierr)
	if (ierr /= 0) then
	if ( mpi_node() == 0 ) then
		write(0,*)  ""
		write(0,*)  "   Error reading zpulse_planewaves parameters"
		write(0,*)  "   aborting..."
	endif
	stop
	endif

	! ----------------------------------------------------------------------
	! save fundamental parameters
	this%a0             = a0
	this%omega0         = omega0
	this%focus          = focus
	this%phase0         = real( phase * pi_180, p_double )
	this%tolerance      = tolerance
	this%wall_injection = wall_injection

	if ((pol_type <-1) .or. (pol_type >1)) then
		if ( mpi_node() == 0 ) then
			write(0,*)  ""
			write(0,*)  "   Error in zpulse parameters"
			write(0,*)  "   pol must be in the range [-1,0,+1]"
			write(0,*)  "   aborting..."
		endif
		stop
	endif
	this%pol_type  = pol_type

	select case ( trim(preset) )
	case ("+x1x2"); angles(1) =    0; angles(2) =   0; pol =   0
	case ("-x1x2"); angles(1) =  180; angles(2) =   0; pol = 180
	case ("+x1x3"); angles(1) =    0; angles(2) =   0; pol =  90
	case ("-x1x3"); angles(1) =  180; angles(2) =   0; pol =  90
	case ("+x2x1"); angles(1) =   90; angles(2) =   0; pol = 180
	case ("-x2x1"); angles(1) =  270; angles(2) =   0; pol =   0
	case ("+x2x3"); angles(1) =   90; angles(2) =   0; pol =  90
	case ("-x2x3"); angles(1) =  270; angles(2) =   0; pol =  90
	case ("+x3x1"); angles(1) =    0; angles(2) =  90; pol = 270
	case ("-x3x1"); angles(1) =    0; angles(2) = -90; pol =  90
	case ("+x3x2"); angles(1) =    0; angles(2) =  90; pol =   0
	case ("-x3x2"); angles(1) =    0; angles(2) = -90; pol =   0
	case ("0"); continue
	case default
		if ( mpi_node() == 0 ) then
			write(0,*)  ''
			write(0,*)  '   Error in zpulse_planewaves parameters'
			write(0,*)  '   preset input does not exist '
			write(0,*)  '   aborting...'
		endif
		stop
	end select
	
	this%angles(1)    	= real( angles(1) * pi_180, p_double )
	this%angles(2)    	= real( angles(2) * pi_180, p_double )
	this%pol          	= real( pol       * pi_180, p_double )	

	this%paraxial     	= paraxial
	this%num_waves    	= num_waves
	this%normalize    	= .false.
	this%normalization	= 1.0_p_double

	! ----------------------------------------------------------------------
	! save longitudinal profile parameters
	this%lon_position = lon_position

	select case ( trim(lon_type) )
	case ("plane")
		this%lon_type = k_lon_plane
		this%lon_int_range(1) = this%omega0
		this%lon_int_range(2) = this%omega0 + 1
		this%num_waves(1) = 1
	case ("gaussian")
		this%lon_type = k_lon_gaussian

		if ( lon_fwhm > 0.0_p_double  ) then
			this%lon_duration   = lon_fwhm/sqrt( 2.0_p_double * log(2.0_p_double))
		else
			this%lon_duration   = lon_duration
		endif
		this%lon_range          = lon_range

		if ( lon_range == -huge(1.0_p_double) ) then
			if ( mpi_node() == 0 ) then
			write(0,*)  ''
			write(0,*)  '   Error in zpulse_planewaves parameters'
			write(0,*)  '   When using "gaussian" longitudinal envelopes, '
			write(0,*)  '   lon_range must also be defined'
			write(0,*)  '   aborting...'
			endif
			stop
		endif

		if ( lon_int_range(1) == 0 .and. lon_int_range(2) == 0 ) then
			sig = 2.0_p_double*sqrt( 2.0_p_double * log(2.0_p_double)) / this%lon_duration
			this%lon_int_range = (/ this%omega0 - n_sig_gaussian*sig, this%omega0 + n_sig_gaussian*sig /)
		else 
			this%lon_int_range(1) = lon_int_range(1)
			this%lon_int_range(2) = lon_int_range(2)
		endif

		! check periodicity
		dlon = ( this%lon_int_range(2) - this%lon_int_range(1) ) / ( this%num_waves(1) - 1 )
		if ( 2*pi/dlon < 2*this%lon_duration ) then
			if ( mpi_node() == 0 ) then
				write(0,*)  ''
				write(0,*)  '   Error in zpulse_planewaves parameters'
				write(0,*)  '   not be able to correctly resolve pulse in longitudinal direction'
				write(0,*)  '   consider more waves than the minimum'
				write(0,*)  '   minimum num_waves(1) =', int( (2*this%lon_duration) * ( this%lon_int_range(2) - this%lon_int_range(1) ) / 2 / pi ) + 2
				write(0,*)  '   aborting...'
			endif
			stop
		endif

	case ("sin2")
		this%lon_type = k_lon_sin2

		if ( lon_fwhm > 0.0_p_double ) then
			this%lon_rise       = lon_fwhm/2.0_p_double
			this%lon_flat       = 0.0_p_double
			this%lon_fall       = lon_fwhm/2.0_p_double
		else
			this%lon_rise       = lon_rise
			this%lon_flat       = lon_flat
			this%lon_fall       = lon_fall
		endif

		if ( lon_int_range(1) == 0 .and. lon_int_range(2) == 0 ) then
			sig = 2.0_p_double*pi / (this%lon_rise + 0.1*this%lon_flat + this%lon_fall)
			this%lon_int_range = (/ this%omega0 - n_sig_sin2*sig, this%omega0 + n_sig_sin2*sig /)
		else 
			this%lon_int_range(1) = lon_int_range(1)
			this%lon_int_range(2) = lon_int_range(2)
		endif

		! check periodicity
		dlon = ( this%lon_int_range(2) - this%lon_int_range(1) ) / ( this%num_waves(1) - 1 )
		if ( 2*pi/dlon < this%lon_rise + this%lon_flat + this%lon_fall ) then
			if ( mpi_node() == 0 ) then
				write(0,*)  ''
				write(0,*)  '   Error in zpulse_planewaves parameters'
				write(0,*)  '   not be able to correctly resolve pulse in longitudinal direction'
				write(0,*)  '   consider more waves than the minimum'
				write(0,*)  '   minimum num_waves(1) =', int( (this%lon_rise + this%lon_flat + this%lon_fall) * ( this%lon_int_range(2) - this%lon_int_range(1) ) / 2 / pi ) + 2
				write(0,*)  '   aborting...'
			endif
			stop
		endif

	case default
		if ( mpi_node() == 0 ) then
			write(0,*)  ''
			write(0,*)  '   Error in zpulse_planewaves parameters'
			write(0,*)  '   lon_type must be "gaussian" or "sin2" '
			write(0,*)  '   no other type has been implemented yet '
			write(0,*)  '   aborting...'
		endif
		stop
	end select

	! ----------------------------------------------------------------------
	! save perpendicular profile parameters
	select case ( trim(per_type) )
	case ("plane")
		this%per_type = k_per_plane
		this%per1_int_range(1) = 0.0_p_double
		this%per1_int_range(2) = 1.0_p_double
		this%per2_int_range(1) = 0.0_p_double
		this%per2_int_range(2) = 1.0_p_double
		this%num_waves(2) = 1
		this%num_waves(3) = 1
	case ("gaussian")
		this%per_type  = k_per_gaussian

		! get main spot size and focal plane
		this%per_w0 = 0.0_p_double
		if (per_fwhm(1) > 0.0_p_double) then
			this%per_w0 = per_fwhm(1)/sqrt(4.0_p_double * log(2.0_p_double))
		else
			this%per_w0 = per_w0(1)
		endif

		! set integration ranges
		if ( per1_int_range(1) == 0 .and. per1_int_range(2) == 0 ) then
			sig = sqrt( 2.0_p_double ) / this%per_w0(1)
			this%per1_int_range = (/ - n_sig_gaussian*sig, n_sig_gaussian*sig /)
		else 
			this%per1_int_range = per1_int_range
		endif

		if ( per2_int_range(1) == 0 .and. per2_int_range(2) == 0 ) then
			sig = sqrt( 2.0_p_double ) / this%per_w0(1)
			this%per2_int_range = (/ - n_sig_gaussian*sig, n_sig_gaussian*sig /)
		else 
			this%per2_int_range = per2_int_range
		endif
		

	case ("laguerre")
		this%per_type = k_per_laguerre_gaussian
		this%per_tem_mode = 0
		this%per_tem_mode = per_tem_mode
		if ( per_tem_mode(1) > 7 .or. per_tem_mode(2) > 7 ) then
			if ( mpi_node() == 0 ) then
				print *, ''
				print *, '   Error in zpulse parameters'
				print *, '   Only Laguerre modes up to 7 have been implemented. Please contact the'
				print *, '   development team if you need higher order modes.'
				print *, '   aborting...'
			endif
			stop
		endif

		! get main spot size and focal plane
		this%per_w0 = 0.0_p_double
		if (per_fwhm(1) > 0.0_p_double) then
			this%per_w0 = per_fwhm(1)/sqrt(4.0_p_double * log(2.0_p_double))
		else
			this%per_w0 = per_w0(1)
		endif

		! set integration ranges
		if ( per1_int_range(1) == 0 .and. per1_int_range(2) == 0 ) then
			sig = sqrt( 2.0_p_double ) / this%per_w0(1)
			this%per1_int_range = (/ - n_sig_gaussian*sig, n_sig_gaussian*sig /)
		else 
			this%per1_int_range = per1_int_range
		endif

		if ( per2_int_range(1) == 0 .and. per2_int_range(2) == 0 ) then
			sig = sqrt( 2.0_p_double ) / this%per_w0(1)
			this%per2_int_range = (/ - n_sig_gaussian*sig, n_sig_gaussian*sig /)
		else 
			this%per2_int_range = per2_int_range
		endif

		! take care of normalization
		if ( this%per_tem_mode(1) == 0 .and. this%per_tem_mode(2) == 0 ) then
			this%normalization = 4*pi**2
		else if ( this%per_tem_mode(1) == 0 .and. abs(this%per_tem_mode(2)) == 1 ) then
			this%normalization = 2*sqrt(2.0_p_double)*pi**(5/2)
		else
			this%normalize = .true.
		endif

	case default
		if ( mpi_node() == 0 ) then
			write(0,*)  ''
			write(0,*)  '   Error in zpulse_planewaves parameters'
			write(0,*)  '   per_type must be "gaussian'
			write(0,*)  '   no other type has been implemented yet '
			write(0,*)  '   aborting...'
		endif
		stop
	end select

	! extra parameters
	this%if_launch      = if_launch
	this%launch_time    = launch_time

	select case ( trim(pulse_type) )
	case ("time")
		this%abouraddy_angle = 0.0_p_double
	case ("space")
		this%abouraddy_angle = pi / 2
	case ("abouraddy")
		if ( abouraddy_angle == 45 ) then
			if ( mpi_node() == 0 ) then
				write(0,*)  ''
				write(0,*)  '   Error in zpulse_planewaves parameters'
				write(0,*)  '   abouraddy_angle cannot be 45 '
				write(0,*)  '   aborting...'
			endif
			stop
		endif
		if ( abouraddy_angle < 0 .or. abouraddy_angle > 180 ) then
			if ( mpi_node() == 0 ) then
				write(0,*)  ''
				write(0,*)  '   Error in zpulse_planewaves parameters'
				write(0,*)  '   abouraddy_angle has to be in range [0; 180] '
				write(0,*)  '   aborting...'
			endif
			stop
		endif
		this%abouraddy_angle = real( abouraddy_angle * pi_180, p_double )	
	case default
		if ( mpi_node() == 0 ) then
			write(0,*)  ''
			write(0,*)  '   Error in zpulse_planewaves parameters'
			write(0,*)  '   pulse_type must be "time", "space" or "abouraddy" '
			write(0,*)  '   aborting...'
		endif
		stop
	end select

	! adapt integration range according to pulse type
	if ( this%per_type == k_per_plane .and. this%lon_type == k_lon_plane ) then
		continue
	else if ( per1_int_range(1) == 0 .and. per1_int_range(2) == 0 .and. per2_int_range(1) == 0 .and. per2_int_range(2) == 0 .and.lon_int_range(1) == 0 .and. lon_int_range(2) == 0 ) then
			
		dper = min( this%per1_int_range(2), this%omega0 )
		dlon = this%lon_int_range(2) - this%omega0

		if ( this%abouraddy_angle < pi/4 .or. this%abouraddy_angle > pi/2 ) then	
			if (  dper**2*(tan(this%abouraddy_angle)**2-1) + (1-tan(this%abouraddy_angle))**2*this%omega0**2 > 0 ) then
				this%lon_int_range(1) = ( (1-tan(this%abouraddy_angle))*this%omega0*tan(this%abouraddy_angle) + sqrt( dper**2*(tan(this%abouraddy_angle)**2-1) + (1-tan(this%abouraddy_angle))**2*this%omega0**2 ) ) / (1-tan(this%abouraddy_angle)**2) - dlon
			else 
				this%lon_int_range(1) = 0
			endif
			this%lon_int_range(2) = this%omega0 + dlon
		else 
			this%lon_int_range(1) = this%omega0 - dlon
			this%lon_int_range(2) = ( (1-tan(this%abouraddy_angle))*this%omega0*tan(this%abouraddy_angle) - sqrt( dper**2*(tan(this%abouraddy_angle)**2-1) + (1-tan(this%abouraddy_angle))**2*this%omega0**2 ) ) / (1-tan(this%abouraddy_angle)**2) + dlon
		endif

		if ( this%lon_int_range(1) < 0 ) then
			this%lon_int_range(1) = 1e-5 * this%omega0
		endif

	endif

	! normalize if needed 
	if ( this%normalize ) then 
		if ( mpi_node() == 0 ) then
			write(0,*)  ''
			write(0,*)  '   Error in zpulse_planewaves parameters'
			write(0,*)  '   normalization not yet implemented '
			write(0,*)  '   aborting...'
		endif
		!stop
	endif

	! Boost pulse
	this%if_boost = .false.
	if (sim_options%gamma > 1.0_p_double) then

		this%gamma = sim_options%gamma
		SCR_ROOT('   boosting zpulse_planewaves with gamma = ', this%gamma )
		this%if_boost = .true.

	endif

	! Write periodicity
	SCR_ROOT(" - x1 periodicity:", 2*pi/( this%lon_int_range(2) - this%lon_int_range(1) ) * ( this%num_waves(1) - 1 ))
	SCR_ROOT(" - x2 periodicity:", 2*pi/( this%per1_int_range(2) - this%per1_int_range(1) ) * ( this%num_waves(2) - 1 ))
	SCR_ROOT(" - x3 periodicity:", 2*pi/( this%per2_int_range(2) - this%per2_int_range(1) ) * ( this%num_waves(3) - 1 ))


end subroutine read_input_zpulse_planewaves
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Check if parameters are good to use with current dimensionality
subroutine check_dimensionality_zpulse_planewaves( this ) 

	implicit none 

	class( t_zpulse_planewaves ), intent(inout) :: this

	if ( p_x_dim == 1) then
		if ( mpi_node() == 0 ) then
			write(0,*)  ''
			write(0,*)  '   Error in zpulse_planewaves dimensionality'
			write(0,*)  '   not usable in 1d'
			write(0,*)  '   aborting...'
		endif
		stop
	endif

	if ( p_x_dim < 3 .and. this%per_type == k_per_laguerre_gaussian ) then
		if ( mpi_node() == 0 ) then
			write(0,*)  ''
			write(0,*)  '   Error in zpulse_planewaves dimensionality'
			write(0,*)  '   laguerre-gaussian only usable in 3d'
			write(0,*)  '   aborting...'
		endif
		stop
	endif

end subroutine check_dimensionality_zpulse_planewaves
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Launch 2D pulse routine non paraxial
subroutine launch_zpulse_planewaves_2d_normal( this, e, b, t, dt, nx_p_min, g_space )

	implicit none

	integer, parameter                   :: rank = 2

	class( t_zpulse_planewaves ), intent(inout) :: this
	type( t_vdf ), intent(inout)         :: e, b
	real(p_double), intent(in)           :: t, dt
	integer, dimension(:), intent(in)    :: nx_p_min
	type( t_space ), intent(in)          :: g_space

	! local variables
	real(p_double), dimension(2,rank)    :: g_x_range
	real(p_double), dimension(rank)      :: ldx

	real(p_double)                       :: x1min, x1, x1_2, dx1_2
	real(p_double)                       :: x2min, x2, x2_2, dx2_2
	real(p_double) 						 :: max_amp, sqrt2

	! other variables
	complex(p_double)                    :: jones_vec1, jones_vec2
	real(p_double), dimension(rank)      :: x_per, xper
	real(p_double), dimension(rank,rank) :: rot
	real(p_double)                       :: k_lon, k_per, dk_lon, dk_per, k1, k2, omega
	complex(p_double) 					 :: e1_fac, e2_fac, e3_fac, b1_fac, b2_fac, b3_fac
	complex(p_double)                    :: f_k
	real(p_double) 						 :: amp
	integer                              :: ik1, ik2, i1, i2
	complex(p_double)  					 :: exp_x1, exp_x1_2, exp_x2, exp_x2_2
	real(p_double) 						 :: beta, beta1, beta2, beta3
	complex(p_double) 					 :: w_b, k1_b, k2_b, e1_b, e2_b, e3_b, b1_b, b2_b, b3_b

	! executable statements
	call get_x_bnd( g_space, g_x_range )
	ldx(1) = b%dx_(1)
	ldx(2) = b%dx_(2)

	! spatial coordinates
	! XMIN ALREADY INCLUDES FOCUS POSITION
	x1min = g_x_range(p_lower,1) + real(nx_p_min(1)-1, p_double)*ldx(1) - this%focus(1)
	x2min = g_x_range(p_lower,2) + real(nx_p_min(2)-1, p_double)*ldx(2) - this%focus(2)
	dx1_2 = ldx(1)/2.0_p_double
	dx2_2 = ldx(2)/2.0_p_double

	! trignometry
	rot(1,:) = (/ cos( this%angles(1) ), -sin( this%angles(1) ) /)
	rot(2,:) = (/ sin( this%angles(1) ),  cos( this%angles(1) ) /)
	sqrt2 = sqrt( 2.0_p_double )

	! testing, debug and control
	if ( control .and. mpi_node() == 0 ) then
		write(0,*) "<control> ", "lon_int_range", this%lon_int_range
		write(0,*) "<control> ", "per1_int_range", this%per1_int_range
	endif

	! plane wave coordinates
	if ( this%num_waves(1) > 1 ) then
		dk_lon = ( this%lon_int_range(2)  - this%lon_int_range(1)  ) / (this%num_waves(1) - 1)
	else; dk_lon = 1.0_p_double; endif
	if ( this%num_waves(2) > 1 ) then
		dk_per = ( this%per1_int_range(2) - this%per1_int_range(1) ) / (this%num_waves(2) - 1)
	else; dk_per = 1.0_p_double; endif
	! calculate maximum amplitude for tolerance
	! assuming the largest contribution is the principal wave ...
	max_amp = abs( this%a0 * this%k_function( this%omega0, 0.0_p_double, 0.0_p_double ) * dk_lon * dk_per )

	! iterate for every plane wave
	do ik1 = 0, this%num_waves(1) - 1 
		! calculate direction of each wave - longitudinal
		k_lon = this%lon_int_range(1)  + real(ik1) * dk_lon
	do ik2 = 0, this%num_waves(2) - 1
		! calculate direction of each wave - perpendicular
		k_per = this%per1_int_range(1) + real(ik2) * dk_per
	! --------------------------------------------------

		! calculate direction of wave in rotated space
		k1    = k_lon * rot(1,1) + k_per * rot(1,2)
		k2    = k_lon * rot(2,1) + k_per * rot(2,2)

		! calculate omega corresponding to this wave
		omega = sqrt( k_lon**2 + k_per**2 )
		!omega = 2/dt * asin( dt * sqrt( sin(k1*ldx(1)/2)**2/ldx(1)**2 + sin(k2*ldx(2)/2)**2/ldx(2)**2 ) )

		! calculate intensity of each wave
		f_k = this%a0 * this%k_function( k_lon, k_per, 0.0_p_double ) * dk_lon * dk_per
			
		! tolerance for speep up
		if ( abs(f_k) < this%tolerance * max_amp ) then
			! testing, debug and control
			if ( control .and. mpi_node() == 0 ) then
				write(0,*) "<control> ", "cycled (tolerance) >", k_lon, k_per
			endif
			cycle 
		endif

		! get perpendicular unitary direction
		x_per = (/ -k_per/sqrt(k_lon**2+k_per**2), k_lon/sqrt(k_lon**2+k_per**2) /)

		! calculate direction of perpendicular unitary direction in rotated space
		xper = (/ x_per(1)*rot(1,1)+x_per(2)*rot(1,2), x_per(1)*rot(2,1)+x_per(2)*rot(2,2) /)

		! get jones vector for polarization
		call get_jones_vector( this, k_lon, k_per, 0.0_p_double, jones_vec1, jones_vec2 )

		! calculate component for each field
		e1_fac = ( ( xper(1)                    ) * jones_vec1 + ( xper(1)                    ) * jones_vec2 ) / sqrt2
		e2_fac = ( ( xper(2)                    ) * jones_vec1 + ( xper(2)                    ) * jones_vec2 ) / sqrt2
		e3_fac = ( (          + im_unit         ) * jones_vec1 + (          - im_unit         ) * jones_vec2 ) / sqrt2
		b1_fac = ( (          - im_unit*xper(1) ) * jones_vec1 + (          + im_unit*xper(1) ) * jones_vec2 ) / sqrt2
		b2_fac = ( (          - im_unit*xper(2) ) * jones_vec1 + (          + im_unit*xper(2) ) * jones_vec2 ) / sqrt2
		b3_fac = ( ( 1.0_p_double               ) * jones_vec1 + ( 1.0_p_double               ) * jones_vec2 ) / sqrt2

		! boost pulse
		if (this%if_boost) then
			! save betas
			beta1 = sqrt( 1 - 1/this%gamma**2 )
			beta2 = 0.0_p_double
			beta3 = 0.0_p_double
			beta  = sqrt( beta1**2 + beta2**2 + beta3**2 )
			! boost k-four-vector
			w_b = omega; k1_b = k1; k2_b = k2
			omega = this%gamma * (w_b - beta1*k1_b - beta2*k2_b)
			k1    = k1 + (this%gamma-1) * (beta1*k1_b + beta2*k2_b) * beta1 / beta**2 - this%gamma * w_b * beta1
			k2    = k2 + (this%gamma-1) * (beta1*k1_b + beta2*k2_b) * beta2 / beta**2 - this%gamma * w_b * beta2
			! boost field factors
			e1_b = e1_fac; e2_b = e2_fac; e3_b = e3_fac; b1_b = b1_fac; b2_b = b2_fac; b3_b = b3_fac;
			e1_fac = this%gamma * ( e1_b + beta2*b3_b - beta3*b2_b ) - (this%gamma-1) * (e1_b*beta1 + e2_b*beta2 + e3_b*beta3) * beta1 / beta**2
			e2_fac = this%gamma * ( e2_b - beta1*b3_b + beta3*b1_b ) - (this%gamma-1) * (e1_b*beta1 + e2_b*beta2 + e3_b*beta3) * beta2 / beta**2
			e3_fac = this%gamma * ( e3_b + beta1*b2_b - beta2*b1_b ) - (this%gamma-1) * (e1_b*beta1 + e2_b*beta2 + e3_b*beta3) * beta3 / beta**2
			b1_fac = this%gamma * ( b1_b - beta2*e3_b + beta3*e2_b ) - (this%gamma-1) * (b1_b*beta1 + b2_b*beta2 + b3_b*beta3) * beta1 / beta**2
			b2_fac = this%gamma * ( b2_b + beta1*e3_b - beta3*e1_b ) - (this%gamma-1) * (b1_b*beta1 + b2_b*beta2 + b3_b*beta3) * beta2 / beta**2
			b3_fac = this%gamma * ( b3_b - beta1*e2_b + beta2*e1_b ) - (this%gamma-1) * (b1_b*beta1 + b2_b*beta2 + b3_b*beta3) * beta3 / beta**2
		endif

		! start integration scheme
		f_k = f_k * omega * exp( - im_unit * omega * (t + this%lon_position) + im_unit * this%phase0 )

		! iterate for every space point
		do i1 = b%lbound(2), b%ubound(2)
			! calculate grid positions and exponentials for x1
			x1 = x1min + real(i1-1, p_double) * ldx(1)
			x1_2 = x1 + dx1_2
			exp_x1   = exp( im_unit * k1*x1   )
			exp_x1_2 = exp( im_unit * k1*x1_2 )
		do i2 = b%lbound(3), b%ubound(3)
			! calculate grid positions and exponentials for x2
			x2 = x2min + real(i2-1, p_double) * ldx(2)
			x2_2 = x2 + dx2_2
			exp_x2   = exp( im_unit * k2*x2   )
			exp_x2_2 = exp( im_unit * k2*x2_2 )
		! --------------------------------------------------

			! e3 at (x1, x2)
			e%f2(3,i1,i2) = e%f2(3,i1,i2) + real( f_k * exp_x1   * exp_x2   * e3_fac, p_double)

			! b3 at (x1_2, x2_2)
			b%f2(3,i1,i2) = b%f2(3,i1,i2) + real( f_k * exp_x1_2 * exp_x2_2 * b3_fac, p_double)

			! e1, b2 at (x1_2, x2)
			e%f2(1,i1,i2) = e%f2(1,i1,i2) + real( f_k * exp_x1_2 * exp_x2   * e1_fac, p_double)
			b%f2(2,i1,i2) = b%f2(2,i1,i2) + real( f_k * exp_x1_2 * exp_x2   * b2_fac, p_double)

			! e2, b1 at (x1, x2_2)
			e%f2(2,i1,i2) = e%f2(2,i1,i2) + real( f_k * exp_x1   * exp_x2_2 * e2_fac, p_double)
			b%f2(1,i1,i2) = b%f2(1,i1,i2) + real( f_k * exp_x1   * exp_x2_2 * b1_fac, p_double)

		enddo
		enddo

	enddo
	! testing, debug and control
	if ( control .and. mpi_node() == 0 ) then
		write(0,*) "<control> ", "injected long ", ik1, "of", this%num_waves(1)
	endif
	enddo

end subroutine launch_zpulse_planewaves_2d_normal
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Launch 3D pulse routine
subroutine launch_zpulse_planewaves_3d_normal( this, e, b, t, dt, nx_p_min, g_space )

	implicit none

	integer, parameter                   :: rank = 3

	class( t_zpulse_planewaves ), intent(inout) :: this
	type( t_vdf ), intent(inout)         :: e, b
	real(p_double), intent(in)           :: t, dt
	integer, dimension(:), intent(in)    :: nx_p_min
	type( t_space ), intent(in)          :: g_space

	! local variables
	real(p_double), dimension(2,rank)    :: g_x_range
	real(p_double), dimension(rank)      :: ldx

	real(p_double)                       :: x1min, x1, x1_2, dx1_2
	real(p_double)                       :: x2min, x2, x2_2, dx2_2
	real(p_double)                       :: x3min, x3, x3_2, dx3_2
	real(p_double) 						 :: max_amp, sqrt2

	! other variables
	complex(p_double)                    :: jones_vec1, jones_vec2
	real(p_double), dimension(rank,rank) :: rot
	real(p_double), dimension(rank)		 :: x_per1, x_per2, xper1, xper2
	real(p_double)                       :: k_long, k_per1, k_per2, dk_long, dk_per1, dk_per2, omega, k1, k2, k3
	complex(p_double)                    :: f_k
	complex(p_double) 					 :: e1_fac, e2_fac, e3_fac, b1_fac, b2_fac, b3_fac
	integer                              :: ik1, ik2, ik3, i1, i2, i3
	complex(p_double)  					 :: exp_x1, exp_x1_2, exp_x2, exp_x2_2, exp_x3, exp_x3_2
	real(p_double) 						 :: aux_boost, beta, beta1, beta2, beta3
	complex(p_double) 					 :: w_b, k1_b, k2_b, k3_b, e1_b, e2_b, e3_b, b1_b, b2_b, b3_b

	! executable statements
	call get_x_bnd( g_space, g_x_range )
	ldx(1) = b%dx_(1)
	ldx(2) = b%dx_(2)
	ldx(3) = b%dx_(3)

	! spatial coordinates
	! XMIN ALREADY INCLUDES FOCUS POSITION
	x1min = g_x_range(p_lower,1) + real(nx_p_min(1)-1, p_double)*ldx(1) - this%focus(1)
	x2min = g_x_range(p_lower,2) + real(nx_p_min(2)-1, p_double)*ldx(2) - this%focus(2)
	x3min = g_x_range(p_lower,3) + real(nx_p_min(3)-1, p_double)*ldx(3) - this%focus(3)
	dx1_2 = ldx(1)/2.0_p_double
	dx2_2 = ldx(2)/2.0_p_double
	dx3_2 = ldx(3)/2.0_p_double

	! trignometry 
	rot(1,:) = (/ cos(this%angles(1))*cos(this%angles(2)), -sin(this%angles(1)), -cos(this%angles(1))*sin(this%angles(2)) /)
	rot(2,:) = (/ sin(this%angles(1))*cos(this%angles(2)),  cos(this%angles(1)), -sin(this%angles(1))*sin(this%angles(2)) /)
	rot(3,:) = (/ sin(this%angles(2)), 0.0_p_double, cos(this%angles(2)) /)
	sqrt2 = sqrt( 2.0_p_double )

	! plane wave coordinates
	dk_long = (  this%lon_int_range(2) -  this%lon_int_range(1) ) / (this%num_waves(1) - 1)
	dk_per1 = ( this%per1_int_range(2) - this%per1_int_range(1) ) / (this%num_waves(2) - 1)
	dk_per2 = ( this%per2_int_range(2) - this%per2_int_range(1) ) / (this%num_waves(3) - 1)
	
	! calculate maximum amplitude for tolerance
	! assuming the largest contribution is the principal wave ...
	max_amp = abs( this%a0 * this%k_function( this%omega0, 0.0_p_double, 0.0_p_double ) * dk_long * dk_per1 * dk_per2 )

	! iterate for every plane wave
	do ik1 = 0, this%num_waves(1) - 1 
		! calculate direction of each wave - longitudinal 
		k_long =  this%lon_int_range(1) + real(ik1) * dk_long
	do ik2 = 0, this%num_waves(2) - 1
		! calculate direction of each wave - perpendicular 1 
		k_per1 = this%per1_int_range(1) + real(ik2) * dk_per1
	do ik3 = 0, this%num_waves(3) - 1
		! calculate direction of each wave - perpendicular 2 
		k_per2 = this%per2_int_range(1) + real(ik3) * dk_per2
	! --------------------------------------------------
			
		! calculate direction of wave in rotated space
		k1 = k_long * rot(1,1) + k_per1 * rot(1,2) + k_per2 * rot(1,3)
		k2 = k_long * rot(2,1) + k_per1 * rot(2,2) + k_per2 * rot(2,3)
		k3 = k_long * rot(3,1) + k_per1 * rot(3,2) + k_per2 * rot(3,3)

		! calculate omega corresponding to this wave
		omega = sqrt( k_long**2 + k_per1**2 + k_per2**2 )
		!omega = 2/dt * asin( dt * sqrt( sin(k1*ldx(1)/2)**2/ldx(1)**2 + sin(k2*ldx(2)/2)**2/ldx(2)**2 + sin(k3*ldx(3)/2)**2/ldx(3)**2 ) )

		! calculate intensity of each wave
		f_k = this%a0 * this%k_function( k_long, k_per1, k_per2 ) * dk_long * dk_per1 * dk_per2
		
		! tolerance for speep up
		if ( abs(f_k) < this%tolerance * max_amp ) then
			! testing, debug and control
			if ( control .and. mpi_node() == 0 ) then
				write(0,*) "<control> ", "cycled (tolerance) >", k_long, k_per1, k_per2
			endif
			cycle 
		endif

		! get perpendicular unitary directions
		x_per1 = (/ -k_per1/sqrt(k_long**2+k_per1**2), k_long/sqrt(k_long**2+k_per1**2), 0.0_p_double /)
		x_per2 = (/ -k_long*k_per2/sqrt(k_long**2+k_per1**2)/omega, -k_per1*k_per2/sqrt(k_long**2+k_per1**2)/omega, sqrt(k_long**2+k_per1**2)/omega /)

		! calculate direction of perpendicular unitary directions in rotated space
		xper1 = (/ x_per1(1)*rot(1,1)+x_per1(2)*rot(1,2)+x_per1(3)*rot(1,3), x_per1(1)*rot(2,1)+x_per1(2)*rot(2,2)+x_per1(3)*rot(2,3), x_per1(1)*rot(3,1)+x_per1(2)*rot(3,2)+x_per1(3)*rot(3,3) /)
		xper2 = (/ x_per2(1)*rot(1,1)+x_per2(2)*rot(1,2)+x_per2(3)*rot(1,3), x_per2(1)*rot(2,1)+x_per2(2)*rot(2,2)+x_per2(3)*rot(2,3), x_per2(1)*rot(3,1)+x_per2(2)*rot(3,2)+x_per2(3)*rot(3,3) /)

		! get jones vector for polarization
		call get_jones_vector( this, k_long, k_per1, k_per2, jones_vec1, jones_vec2 )

		! calculate component for each field
		e1_fac = ( ( xper1(1) + im_unit*xper2(1) ) * jones_vec1 + ( xper1(1) - im_unit*xper2(1) ) * jones_vec2 ) / sqrt2
		e2_fac = ( ( xper1(2) + im_unit*xper2(2) ) * jones_vec1 + ( xper1(2) - im_unit*xper2(2) ) * jones_vec2 ) / sqrt2
		e3_fac = ( ( xper1(3) + im_unit*xper2(3) ) * jones_vec1 + ( xper1(3) - im_unit*xper2(3) ) * jones_vec2 ) / sqrt2
		b1_fac = ( ( xper2(1) - im_unit*xper1(1) ) * jones_vec1 + ( xper2(1) + im_unit*xper1(1) ) * jones_vec2 ) / sqrt2
		b2_fac = ( ( xper2(2) - im_unit*xper1(2) ) * jones_vec1 + ( xper2(2) + im_unit*xper1(2) ) * jones_vec2 ) / sqrt2
		b3_fac = ( ( xper2(3) - im_unit*xper1(3) ) * jones_vec1 + ( xper2(3) + im_unit*xper1(3) ) * jones_vec2 ) / sqrt2

		! start integration scheme
		f_k = f_k * omega * exp( - im_unit * omega * (t + this%lon_position) + im_unit * this%phase0 )

		! boost pulse
		if (this%if_boost) then
			! save betas
			beta1 = sqrt( 1 - 1/this%gamma**2 )
			beta2 = 0.0_p_double
			beta3 = 0.0_p_double
			beta  = sqrt( beta1**2 + beta2**2 + beta3**2 )
			! boost k-four-vector
			w_b = omega; k1_b = k1; k2_b = k2; k3_b = k3
			omega = this%gamma * (w_b - beta1*k1_b - beta2*k2_b - beta2*k3_b)
			k1    = k1 + (this%gamma-1) * (beta1*k1_b + beta2*k2_b + beta3*k3_b) * beta1 / beta**2 - this%gamma * w_b * beta1
			k2    = k2 + (this%gamma-1) * (beta1*k1_b + beta2*k2_b + beta3*k3_b) * beta2 / beta**2 - this%gamma * w_b * beta2
			k3    = k3 + (this%gamma-1) * (beta1*k1_b + beta2*k2_b + beta3*k3_b) * beta3 / beta**2 - this%gamma * w_b * beta3
			! boost field factors
			e1_b = e1_fac; e2_b = e2_fac; e3_b = e3_fac; b1_b = b1_fac; b2_b = b2_fac; b3_b = b3_fac;
			e1_fac = this%gamma * ( e1_b + beta2*b3_b - beta3*b2_b ) - (this%gamma-1) * (e1_b*beta1 + e2_b*beta2 + e3_b*beta3) * beta1 / beta**2
			e2_fac = this%gamma * ( e2_b - beta1*b3_b + beta3*b1_b ) - (this%gamma-1) * (e1_b*beta1 + e2_b*beta2 + e3_b*beta3) * beta2 / beta**2
			e3_fac = this%gamma * ( e3_b + beta1*b2_b - beta2*b1_b ) - (this%gamma-1) * (e1_b*beta1 + e2_b*beta2 + e3_b*beta3) * beta3 / beta**2
			b1_fac = this%gamma * ( b1_b - beta2*e3_b + beta3*e2_b ) - (this%gamma-1) * (b1_b*beta1 + b2_b*beta2 + b3_b*beta3) * beta1 / beta**2
			b2_fac = this%gamma * ( b2_b + beta1*e3_b - beta3*e1_b ) - (this%gamma-1) * (b1_b*beta1 + b2_b*beta2 + b3_b*beta3) * beta2 / beta**2
			b3_fac = this%gamma * ( b3_b - beta1*e2_b + beta2*e1_b ) - (this%gamma-1) * (b1_b*beta1 + b2_b*beta2 + b3_b*beta3) * beta3 / beta**2
		endif

		! iterate for every space point
		do i1 = b%lbound(2), b%ubound(2)
			! calculate grid positions and exponentials for x1
			x1 = x1min + real(i1-1, p_double) * ldx(1)
			x1_2 = x1 + dx1_2
			exp_x1   = exp( im_unit * k1*x1   )
			exp_x1_2 = exp( im_unit * k1*x1_2 )
		do i2 = b%lbound(3), b%ubound(3)
			! calculate grid positions and exponentials for x1
			x2 = x2min + real(i2-1, p_double) * ldx(2)
			x2_2 = x2 + dx2_2
			exp_x2   = exp( im_unit * k2*x2   )
			exp_x2_2 = exp( im_unit * k2*x2_2 )
		do i3 = b%lbound(4), b%ubound(4)
			! calculate grid positions and exponentials for x1
			x3 = x3min + real(i3-1, p_double) * ldx(3)
			x3_2 = x3 + dx3_2
			exp_x3   = exp( im_unit * k3*x3   )
			exp_x3_2 = exp( im_unit * k3*x3_2 )
		! --------------------------------------------------
	
				! e1 at (x1_2, x2  , x3  )
				e%f3(1,i1,i2,i3) = e%f3(1,i1,i2,i3) + real( f_k * exp_x1_2 * exp_x2   * exp_x3   * e1_fac, p_double)
				! e2 at (x1  , x2_2, x3  )
				e%f3(2,i1,i2,i3) = e%f3(2,i1,i2,i3) + real( f_k * exp_x1   * exp_x2_2 * exp_x3   * e2_fac, p_double)
				! e3 at (x1  , x2  , x3_2)
				e%f3(3,i1,i2,i3) = e%f3(3,i1,i2,i3) + real( f_k * exp_x1   * exp_x2   * exp_x3_2 * e3_fac, p_double)
				! b1 at (x1  , x2_2, x3_2)
				b%f3(1,i1,i2,i3) = b%f3(1,i1,i2,i3) + real( f_k * exp_x1   * exp_x2_2 * exp_x3_2 * b1_fac, p_double)
				! b2 at (x1_2, x2  , x3_2)
				b%f3(2,i1,i2,i3) = b%f3(2,i1,i2,i3) + real( f_k * exp_x1_2 * exp_x2   * exp_x3_2 * b2_fac, p_double)
				! b3 at (x1_2, x2_2, x3  )
				b%f3(3,i1,i2,i3) = b%f3(3,i1,i2,i3) + real( f_k * exp_x1_2 * exp_x2_2 * exp_x3   * b3_fac, p_double)

		enddo
		enddo
		enddo

	enddo
	enddo
		! testing, debug and control
		!if ( control .and. mpi_node() == 0 ) then
		!	write(0,*) "<control> ", "injected long ", ik1, "of", this%num_waves(1)
		!endif
	enddo

end subroutine launch_zpulse_planewaves_3d_normal
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Launch pulse routine
subroutine launch_zpulse_planewaves( this, emf, g_space, nx_p_min, g_nx, t, dt, no_co )

	class( t_zpulse_planewaves ), intent(inout)  :: this
	class( t_emf ) ,  intent(inout)              :: emf

	type( t_space ), intent(in)                  :: g_space

	integer, intent(in), dimension(:)            :: nx_p_min, g_nx
	real(p_double), intent(in)                   :: t
	real(p_double), intent(in)                   :: dt
	class( t_node_conf ), intent(in)             :: no_co

	! local variables
	type( t_vdf )                                :: e_pulse, b_pulse
	type( t_space )                              :: lg_space
	integer, dimension(p_x_dim)                  :: lnx_p_min, lg_nx

	! testing, debug and control
	real(p_double) :: start_time, end_time

	! launch zpulse
	if ( this%if_launch .and. ( t >= this%launch_time ) ) then

		if ( control ) then
			call cpu_time(start_time)
		endif

		! initialize local variables
		lg_space = g_space
		lnx_p_min = nx_p_min(1:p_x_dim)
		lg_nx = g_nx(1:p_x_dim)

		! create electric and magnetic field vdfs
		call e_pulse % new( emf%e, zero = .true. )
		call b_pulse % new( emf%b, zero = .true. )

		! choose dimension and if paraxial
		if( p_x_dim == 2 .and. .not. this%paraxial ) then
			call launch_zpulse_planewaves_2d_normal( this, e_pulse, b_pulse, t, dt, nx_p_min, g_space )
		else if( p_x_dim == 2 .and. this%paraxial ) then
			!call launch_zpulse_planewaves_2d_paraxial( this, e_pulse, b_pulse, t, dt, nx_p_min, g_space )
			if ( mpi_node() == 0 ) then
				print *, ''
				print *, '   2D paraxial not yet implemented'
				print *, '   aborting...'
			endif
			stop
		else if( p_x_dim == 3 .and. .not. this%paraxial ) then
			call launch_zpulse_planewaves_3d_normal( this, e_pulse, b_pulse, t, dt, nx_p_min, g_space )
		else if( p_x_dim == 3 .and. this%paraxial ) then
			!call launch_zpulse_planewaves_3d_paraxial( this, e_pulse, b_pulse, t, dt, nx_p_min, g_space )
			if ( mpi_node() == 0 ) then
				print *, ''
				print *, '   3D paraxial not yet implemented'
				print *, '   aborting...'
			endif
			stop
		endif		

		! add the created pulses to the simulation
		call add( emf%e, e_pulse )
		call add( emf%b, b_pulse )

		! clear memory
		call e_pulse % cleanup()
		call b_pulse % cleanup()

		! zpulse has been launched
		this%if_launch = .false.

		! testing, debug and control
		if ( control ) then
			call cpu_time(end_time)
			write(0,*) "<control> ", "Elapsed time for injection in node", mpi_node(), "was", end_time-start_time
		endif

	endif

	! launch zpulse 
	! TODO: inner variable to control this part
	if ( this%wall_injection .and. t > this%launch_time ) then
		call launch_em_planewaves_wall_2d( this, emf%e, emf%b, g_space, nx_p_min, t, dt, no_co )
	endif

end subroutine launch_zpulse_planewaves
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Wave vector weight function
function k_function_zpulse_planewaves( this, k_lon, k_per1, k_per2 ) result( f_k )

	implicit none
	class( t_zpulse_planewaves ), intent(in)   :: this
	real(p_double), intent(in)                 :: k_lon, k_per1, k_per2
	complex(p_double)                          :: f_k

	! auxiliar variables
	real(p_double) 							   :: f_k_lon_arg, vg

	!f_k_lon_arg = sign( sqrt( abs( k_lon**2 + (k_per1**2+k_per2**2) * this%abouraddy_factor ) ), k_lon**2 + (k_per1**2+k_per2**2) * this%abouraddy_factor )
	if ( this%abouraddy_angle == pi / 2 ) then
		f_k_lon_arg = k_lon
	else 
		vg = tan( this%abouraddy_angle )
		f_k_lon_arg = ( -k_lon * vg + sqrt( k_lon**2 + (k_per1**2+k_per2**2) ) ) / ( 1 - vg )
	endif
	f_k = this%k_lon_function( f_k_lon_arg ) * this%k_per_function( k_per1, k_per2 ) 

end function k_function_zpulse_planewaves
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Transverse wave vector weight function
function k_per_function_zpulse_planewaves( this, k_per1, k_per2 ) result( f_k_per )

	implicit none
	class( t_zpulse_planewaves ), intent(in)   :: this
	real(p_double), intent(in)                 :: k_per1, k_per2
	complex(p_double)                          :: f_k_per

	! auxiliar variables
	real(p_double)                             :: sig, kr, ktheta, w0
	integer 								   :: p, l

	select case ( p_x_dim )
		case(2)
			select case ( this%per_type )
			case (k_per_plane)
				f_k_per = 1.0_p_double
			case (k_per_gaussian)

				sig = sqrt( 2.0_p_double ) / this%per_w0(1)
				f_k_per = 1.0 / sqrt( 2.0*pi ) / sig * exp( -k_per1**2/2/sig**2 )

			end select

		case(3)
			select case ( this%per_type )
			case (k_per_plane)
				f_k_per = 1.0_p_double
			case (k_per_gaussian)

				sig = sqrt( 2.0_p_double ) / this%per_w0(1)
				f_k_per = 1.0 / ( 2.0*pi ) / sig**2 * exp( -( k_per1**2 + k_per2**2 )/2/sig**2 )
			
			case (k_per_laguerre_gaussian)

				kr = sqrt( k_per1**2 + k_per2**2 )
				ktheta = atan2( k_per2, k_per1 )
				p = this%per_tem_mode(1)
				l = this%per_tem_mode(2)
				w0 = this%per_w0(1)
				
				f_k_per = w0**2 * (-im_unit)**l * (-1)**(p-l)
				f_k_per = f_k_per * (kr*w0/sqrt(2.0_p_double))**abs(l)
				f_k_per = f_k_per * laguerre( p, abs(l), kr**2*w0**2/2 )
				f_k_per = f_k_per * exp( -kr**2*w0**2/4 + im_unit*l*ktheta )

				f_k_per = f_k_per / this%normalization

			end select

	end select

end function k_per_function_zpulse_planewaves
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Longitudinal wave vector weight function
function k_lon_function_zpulse_planewaves( this, k ) result( f_k_lon )

	implicit none
	class( t_zpulse_planewaves ), intent(in)   :: this
	real(p_double), intent(in)                 :: k
	real(p_double) 							   :: k_lon
	complex(p_double)                          :: f_k_lon

	! auxiliar variables
	real(p_double)                             :: sig, real_part, imag_part

	k_lon = k

	select case ( this%lon_type )
	case (k_lon_plane)
		f_k_lon = 1.0_p_double
	case (k_lon_gaussian)

		sig = 2.0 * sqrt( 2.0 * log(2.0) ) / this%lon_duration
		f_k_lon = 1.0 / sqrt( 2.0*pi ) / sig * exp( -(k_lon-this%omega0)**2/2/sig**2 )

	case (k_lon_sin2)

		if ( k_lon == this%omega0 ) k_lon = k_lon + 0.001

		real_part = (pi**1.5*((-(this%lon_fall**2*(k_lon-this%omega0)**2) + pi**2)*sin(((k_lon-this%omega0)*(this%lon_fall + this%lon_flat - this%lon_rise))/2.) + (-(this%lon_fall**2*(k_lon-this%omega0)**2) + 2*pi**2 - (k_lon-this%omega0)**2*this%lon_rise**2)*sin(((k_lon-this%omega0)*(this%lon_fall + this%lon_flat + this%lon_rise))/2.) + (-pi**2 + (k_lon-this%omega0)**2*this%lon_rise**2)*sin(this%lon_fall*(k_lon-this%omega0) - ((k_lon-this%omega0)*(this%lon_fall + this%lon_flat + this%lon_rise))/2.)))/(2.*sqrt(2.0_p_double)*(k_lon-this%omega0)*(this%lon_fall*(k_lon-this%omega0) - pi)*(this%lon_fall*(k_lon-this%omega0) + pi)*(-pi + (k_lon-this%omega0)*this%lon_rise)*(pi + (k_lon-this%omega0)*this%lon_rise))
		imag_part = (pi**1.5*((this%lon_fall**2*(k_lon-this%omega0)**2 - pi**2)*cos(((k_lon-this%omega0)*(this%lon_fall + this%lon_flat - this%lon_rise))/2.) + (k_lon-this%omega0)**2*(this%lon_fall**2 - this%lon_rise**2)*cos(((k_lon-this%omega0)*(this%lon_fall + this%lon_flat + this%lon_rise))/2.) + (pi**2 - (k_lon-this%omega0)**2*this%lon_rise**2)*cos(this%lon_fall*(k_lon-this%omega0) - ((k_lon-this%omega0)*(this%lon_fall + this%lon_flat + this%lon_rise))/2.)))/(2.*sqrt(2.0_p_double)*(k_lon-this%omega0)*(this%lon_fall*(k_lon-this%omega0) - pi)*(this%lon_fall*(k_lon-this%omega0) + pi)*(-pi + (k_lon-this%omega0)*this%lon_rise)*(pi + (k_lon-this%omega0)*this%lon_rise))
		f_k_lon = (real_part + im_unit * imag_part) / sqrt(2.0_p_double * pi)

	end select

end function k_lon_function_zpulse_planewaves
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Get Jones Vector for plane wave
subroutine get_jones_vector( this, k_lon, k_per1, k_per2, jones_vec1, jones_vec2 )
	
	class( t_zpulse_planewaves ), intent(in)	:: this
	real(p_double), intent(in)					:: k_lon, k_per1, k_per2
	complex(p_double), intent(out)				:: jones_vec1, jones_vec2
  
	! NOT YET WORKING FOR POINCARE BEAMS
	jones_vec1 = cos( pi/4.0 * (1.0 - this%pol_type) )
	jones_vec2 = sin( pi/4.0 * (1.0 - this%pol_type) ) * exp( -im_unit*2.0*this%pol )
  
end subroutine get_jones_vector
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Longitudinal wave vector weight function
subroutine launch_em_planewaves_wall_2d( this, e, b, g_space, nx_p_min, t, dt, no_co )
	
	implicit none
	
	integer, parameter :: rank = 2
	
	class( t_zpulse_planewaves ), intent(inout) :: this
	type( t_vdf ), intent(inout)		  :: e           ! magnetic field
	type( t_vdf ), intent(inout)		  :: b           ! magnetic field
	type( t_space ), intent(in)			  :: g_space     ! global space information
	integer, intent(in), dimension(:)	  :: nx_p_min
	real( p_double ), intent(in)		  :: t           ! simulation time
	real( p_double ), intent(in)		  :: dt          ! time step
	class( t_node_conf ), intent(in)	  :: no_co       ! node configuration

	real(p_double), dimension(2,rank)	 :: g_x_range
	real(p_double), dimension(rank)      :: ldx

	integer 							 :: ip, inj
	real(p_double) 						 :: k_angle

	real(p_double)                       :: x1min, x1, x1_2, dx1_2
	real(p_double)                       :: x2min, x2, x2_2, dx2_2
	real(p_double) 						 :: max_amp

	! other variables
	real(p_double)                       :: cos_pol, sin_pol
	real(p_double), dimension(rank)      :: x_per, xper
	real(p_double), dimension(rank,rank) :: rot
	real(p_double)                       :: k_lon, k_per, dk_lon, dk_per, k1, k2, omega, aux_boost
	complex(p_double) 					 :: e1_fac, e2_fac, e3_fac, b1_fac, b2_fac, b3_fac
	complex(p_double)                    :: f_k
	real(p_double) 						 :: amp
	integer                              :: ik1, ik2, i1, i2, i1min, i1max, i2min, i2max
	complex(p_double)  					 :: exp_x1, exp_x1_2, exp_x2, exp_x2_2
	real(p_double) 						 :: beta, beta1, beta2, beta3
	complex(p_double) 					 :: w_b, k1_b, k2_b, e1_b, e2_b, e3_b, b1_b, b2_b, b3_b

	! executable statements
	call get_x_bnd( g_space, g_x_range )
	ldx(1) = b%dx_(1)
	ldx(2) = b%dx_(2)

	! spatial coordinates
	! XMIN ALREADY INCLUDES FOCUS POSITION
	x1min = g_x_range(p_lower,1) + real(nx_p_min(1)-1, p_double)*ldx(1) - this%focus(1)
	x2min = g_x_range(p_lower,2) + real(nx_p_min(2)-1, p_double)*ldx(2) - this%focus(2)
	dx1_2 = ldx(1)/2.0_p_double
	dx2_2 = ldx(2)/2.0_p_double

	! trignometry
	cos_pol = cos( this%pol )
	sin_pol = sin( this%pol )
	rot(1,:) = (/ cos( this%angles(1) ), -sin( this%angles(1) ) /)
	rot(2,:) = (/ sin( this%angles(1) ),  cos( this%angles(1) ) /)

	! plane wave coordinates
	if ( this%num_waves(1) > 1 ) then
		dk_lon = ( this%lon_int_range(2)  - this%lon_int_range(1)  ) / (this%num_waves(1) - 1)
	else; dk_lon = 1.0_p_double; endif
	if ( this%num_waves(2) > 1 ) then
		dk_per = ( this%per1_int_range(2) - this%per1_int_range(1) ) / (this%num_waves(2) - 1)
	else; dk_per = 1.0_p_double; endif
	! calculate maximum amplitude for tolerance
	! assuming the largest contribution is the principal wave ...
	max_amp = abs( this%a0 * this%k_function( this%omega0, 0.0_p_double, 0.0_p_double ) * dk_lon * dk_per )

	! iterate for every plane wave
	do ik1 = 0, this%num_waves(1) - 1 
		! calculate direction of each wave - longitudinal
		k_lon = this%lon_int_range(1)  + real(ik1) * dk_lon
	do ik2 = 0, this%num_waves(2) - 1
		! calculate direction of each wave - perpendicular
		k_per = this%per1_int_range(1) + real(ik2) * dk_per
	! --------------------------------------------------

		! calculate direction of wave in rotated space
		k1    = k_lon * rot(1,1) + k_per * rot(1,2)
		k2    = k_lon * rot(2,1) + k_per * rot(2,2)

		! calculate omega corresponding to this wave
		omega = sqrt( k_lon**2 + k_per**2 )
		!omega = 2/dt * asin( dt * sqrt( sin(k1*ldx(1)/2)**2/ldx(1)**2 + sin(k2*ldx(2)/2)**2/ldx(2)**2 ) )

		! calculate intensity of each wave
		f_k = this%a0 * this%k_function( k_lon, k_per, 0.0_p_double ) * dk_lon * dk_per
			
		! tolerance for speep up
		if ( abs(f_k) < this%tolerance * max_amp ) cycle

		! get perpendicular unitary direction
		x_per = (/ -k_per/sqrt(k_lon**2+k_per**2), k_lon/sqrt(k_lon**2+k_per**2) /)

		! calculate direction of perpendicular unitary direction in rotated space
		xper = (/ x_per(1)*rot(1,1)+x_per(2)*rot(1,2), x_per(1)*rot(2,1)+x_per(2)*rot(2,2) /)

		! calculate component for each field
		e1_fac = xper(1) * cos_pol - this%pol_type*im_unit * xper(1) * sin_pol
		e2_fac = xper(2) * cos_pol - this%pol_type*im_unit * xper(2) * sin_pol
		e3_fac = this%pol_type*im_unit * cos_pol + sin_pol
		b1_fac = - this%pol_type*im_unit * xper(1) * cos_pol - xper(1) * sin_pol
		b2_fac = - this%pol_type*im_unit * xper(2) * cos_pol - xper(2) * sin_pol
		b3_fac = cos_pol - this%pol_type*im_unit * sin_pol

		! boost pulse
		if (this%if_boost) then
			! save betas
			beta1 = sqrt( 1 - 1/this%gamma**2 )
			beta2 = 0.0_p_double
			beta3 = 0.0_p_double
			beta  = sqrt( beta1**2 + beta2**2 + beta3**2 )
			! boost k-four-vector
			w_b = omega; k1_b = k1; k2_b = k2
			omega = this%gamma * (w_b - beta1*k1_b - beta2*k2_b)
			k1    = k1 + (this%gamma-1) * (beta1*k1_b + beta2*k2_b) * beta1 / beta**2 - this%gamma * w_b * beta1
			k2    = k2 + (this%gamma-1) * (beta1*k1_b + beta2*k2_b) * beta2 / beta**2 - this%gamma * w_b * beta2
			! boost field factors
			e1_b = e1_fac; e2_b = e2_fac; e3_b = e3_fac; b1_b = b1_fac; b2_b = b2_fac; b3_b = b3_fac;
			e1_fac = this%gamma * ( e1_b + beta2*b3_b - beta3*b2_b ) - (this%gamma-1) * (e1_b*beta1 + e2_b*beta2 + e3_b*beta3) * beta1 / beta**2
			e2_fac = this%gamma * ( e2_b - beta1*b3_b + beta3*b1_b ) - (this%gamma-1) * (e1_b*beta1 + e2_b*beta2 + e3_b*beta3) * beta2 / beta**2
			e3_fac = this%gamma * ( e3_b + beta1*b2_b - beta2*b1_b ) - (this%gamma-1) * (e1_b*beta1 + e2_b*beta2 + e3_b*beta3) * beta3 / beta**2
			b1_fac = this%gamma * ( b1_b - beta2*e3_b + beta3*e2_b ) - (this%gamma-1) * (b1_b*beta1 + b2_b*beta2 + b3_b*beta3) * beta1 / beta**2
			b2_fac = this%gamma * ( b2_b + beta1*e3_b - beta3*e1_b ) - (this%gamma-1) * (b1_b*beta1 + b2_b*beta2 + b3_b*beta3) * beta2 / beta**2
			b3_fac = this%gamma * ( b3_b - beta1*e2_b + beta2*e1_b ) - (this%gamma-1) * (b1_b*beta1 + b2_b*beta2 + b3_b*beta3) * beta3 / beta**2
		endif

		! start integration scheme
		f_k = f_k * omega * exp( - im_unit * omega * (t + this%lon_position) + im_unit * this%phase0 )

		! angle of k
        k_angle = atan2( k2, k1 )

		! bottom wall
        if ( my_ngp( no_co, 2) == 1 ) then  ! calculate grid positions and exponentials for x2
            inj = 1
            x2 = x2min + real(inj-1, p_double) * ldx(2)
            x2_2 = x2 + dx2_2
			exp_x2   = exp( im_unit * k2*x2   )
			exp_x2_2 = exp( im_unit * k2*x2_2 )
            do ip = 0, b%nx_(1)+2   ! calculate grid positions and exponentials for x1
                x1 = x1min + real(ip-1, p_double) * ldx(1)
                x1_2 = x1 + dx1_2
				exp_x1   = exp( im_unit * k1*x1   )
				exp_x1_2 = exp( im_unit * k1*x1_2 )
            ! --------------------------------------------------
                e%f2(1, ip, inj) = e%f2(1, ip, inj) + dt/ldx(2) * max(0.0, sin(k_angle)) * real( f_k * exp_x1_2 * exp_x2   * e1_fac, p_double)
                e%f2(2, ip, inj) = e%f2(2, ip, inj) + dt/ldx(2) * max(0.0, sin(k_angle)) * real( f_k * exp_x1   * exp_x2_2 * e2_fac, p_double)
                e%f2(3, ip, inj) = e%f2(3, ip, inj) + dt/ldx(2) * max(0.0, sin(k_angle)) * real( f_k * exp_x1   * exp_x2   * e3_fac, p_double)
                b%f2(1, ip, inj) = b%f2(1, ip, inj) + dt/ldx(2) * max(0.0, sin(k_angle)) * real( f_k * exp_x1   * exp_x2_2 * b1_fac, p_double)
                b%f2(2, ip, inj) = b%f2(2, ip, inj) + dt/ldx(2) * max(0.0, sin(k_angle)) * real( f_k * exp_x1_2 * exp_x2   * b2_fac, p_double)
                b%f2(3, ip, inj) = b%f2(3, ip, inj) + dt/ldx(2) * max(0.0, sin(k_angle)) * real( f_k * exp_x1_2 * exp_x2_2 * b3_fac, p_double)
            enddo
        endif

		! front wall
       	if ( .false. .and. my_ngp( no_co, 1) == nx( no_co, 1 ) ) then ! calculate grid positions and exponentials for x1
			inj = b%nx_(1)
            x1 = x1min + real(inj-1, p_double) * ldx(1)
			x1_2 = x1 + dx1_2
			exp_x1   = exp( im_unit * k1*x1   )
			exp_x1_2 = exp( im_unit * k1*x1_2 )
			do ip = 0, b%nx_(2)+2 ! calculate grid positions and exponentials for x2
				x2 = x2min + real(ip-1, p_double) * ldx(2)
				x2_2 = x2 + dx2_2
				exp_x2   = exp( im_unit * k2*x2   )
				exp_x2_2 = exp( im_unit * k2*x2_2 )
			! --------------------------------------------------
				e%f2(1, inj, ip) = e%f2(1, inj, ip) + (1-cos(k_angle)) * real( f_k * exp_x1_2 * exp_x2   * e1_fac, p_double)
				e%f2(2, inj, ip) = e%f2(2, inj, ip) + (1-cos(k_angle)) * real( f_k * exp_x1   * exp_x2_2 * e2_fac, p_double)
				e%f2(3, inj, ip) = e%f2(3, inj, ip) + (1-cos(k_angle)) * real( f_k * exp_x1   * exp_x2   * e3_fac, p_double)
				b%f2(1, inj, ip) = b%f2(1, inj, ip) + (1-cos(k_angle)) * real( f_k * exp_x1   * exp_x2_2 * b1_fac, p_double)
				b%f2(2, inj, ip) = b%f2(2, inj, ip) + (1-cos(k_angle)) * real( f_k * exp_x1_2 * exp_x2   * b2_fac, p_double)
				b%f2(3, inj, ip) = b%f2(3, inj, ip) + (1-cos(k_angle)) * real( f_k * exp_x1_2 * exp_x2_2 * b3_fac, p_double)
 
			enddo
       	endif


	enddo
	enddo

end subroutine launch_em_planewaves_wall_2d
!-----------------------------------------------------------------------------------------



end module m_zpulse_planewaves







