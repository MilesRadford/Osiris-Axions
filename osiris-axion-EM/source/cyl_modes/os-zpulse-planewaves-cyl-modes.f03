#include "os-config.h"
#include "os-preprocess.fpp"

module m_zpulse_planewaves_cyl_modes

#include "memory/memory.h"

use m_parameters
use m_math
use m_fparser,       only: p_max_expr_len, p_k_fparse, eval
use m_node_conf,     only: t_node_conf
use m_zpulse_planewaves, only : t_zpulse_planewaves
use m_input_file,    only: t_input_file, get_namelist, p_max_nml_section_name
use m_emf_define,    only: t_emf_bound, t_emf
use m_emf_cyl_modes, only: t_emf_cyl_modes
use m_space,         only: t_space, x_bnd, get_x_bnd
use m_vdf_define,    only: t_vdf
use m_grid_define,   only: t_grid
use m_grid_cyl_modes,only: t_grid_cyl_modes
use m_space,         only: xmin, xmax 
use m_vdf_math,      only: add

use m_zpulse_planewaves, only: k_per_gaussian, k_per_laguerre_gaussian
use m_zpulse_planewaves, only: im_unit

! testing, debug and control
logical, parameter :: control = .true.

!-----------------------------------------------------------------------------------------
type, extends( t_zpulse_planewaves ) :: t_zpulse_planewaves_cyl_modes

	! new variables go here
	integer   :: m_pulse

contains
	procedure :: read_input => read_input_zpulse_planewaves_cyl_modes
	procedure :: launch => launch_zpulse_planewaves_cyl_modes
	procedure :: check_dimensionality => check_dimensionality_zpulse_planewaves_cyl_modes
	!procedure :: correct_axis_fields => correct_axis_fields_zpulse_planewaves_cyl_modes

	procedure :: k_function => k_function_zpulse_planewaves_cyl_modes
	procedure :: k_per_function => k_per_function_zpulse_planewaves_cyl_modes

end type t_zpulse_planewaves_cyl_modes
!-----------------------------------------------------------------------------------------

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Reading input
subroutine read_input_zpulse_planewaves_cyl_modes( this, input_file, g_space, bnd_con, periodic, grid, sim_options )

	implicit none

	class( t_zpulse_planewaves_cyl_modes ), intent(inout) :: this
	class( t_input_file ), intent(inout) :: input_file
	type( t_space ), intent(in) :: g_space
	class (t_emf_bound), intent(in) :: bnd_con
	logical, dimension(:), intent(in) :: periodic
	class( t_grid ), intent(in) :: grid
	type( t_options ), intent(in) :: sim_options

	! Call superclass method
	call this % t_zpulse_planewaves % read_input( input_file, g_space, bnd_con, periodic, grid, sim_options )

	! Adjust integration range
	this%per1_int_range(1) = 0.0_p_double

	! No paraxial method available
	if ( this%paraxial ) then
		if ( mpi_node() == 0 ) then
			write(0,*)  ""
			write(0,*)  "   Error reading zpulse_planewaves parameters"
			write(0,*)  "   paraxial approximation not usable when running with cylindrical modes"
			write(0,*)  "   aborting..."
		endif
		stop
	endif

	! Only longitudinal propagation
	if ( ( this%angles(1) /= 0 .and. this%angles(1) /= pi ) .or. this%angles(2) /= 0 ) then
		if ( mpi_node() == 0 ) then
			write(0,*)  ""
			write(0,*)  "   Error reading zpulse_planewaves parameters"
			write(0,*)  "   only longitudinal propagation allowed when running with cylindrical modes"
			write(0,*)  "   aborting..."
		endif
		stop
	endif

	! Focus on the axis
	if ( this%focus(2) /= 0 .or. this%focus(3) /= 0 ) then
		if ( mpi_node() == 0 ) then
			write(0,*)  ""
			write(0,*)  "   Error reading zpulse_planewaves parameters"
			write(0,*)  "   pulse focus must be on axis when running with cylindrical modes"
			write(0,*)  "   aborting..."
		endif
		stop
	endif

	! Take care of normalization
	this%normalization = 1.0_p_double

	! Calculate maximum mode needed
	select case ( this%per_type )
		case (k_per_gaussian)
			this%m_pulse = 1
		case (k_per_laguerre_gaussian)
			this%m_pulse = abs( this%per_tem_mode(2) ) + 1
	end select

end subroutine read_input_zpulse_planewaves_cyl_modes
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Check if parameters are good to use with current dimensionality
subroutine check_dimensionality_zpulse_planewaves_cyl_modes( this ) 

	implicit none 

	class( t_zpulse_planewaves_cyl_modes ), intent(inout) :: this

	! nothing to do

end subroutine check_dimensionality_zpulse_planewaves_cyl_modes
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Launch 2D pulse routine non paraxial
subroutine launch_zpulse_planewaves_cyl_modes_normal( this, e_re, b_re, e_im, b_im, t, nx_p_min, g_space )

	implicit none

	integer, parameter                   :: rank = 2

	class( t_zpulse_planewaves_cyl_modes ), intent(inout) :: this
	type( t_vdf ), dimension(this%m_pulse), intent(inout) :: e_re, b_re, e_im, b_im
	real(p_double), intent(in)           :: t
	integer, dimension(:), intent(in)    :: nx_p_min
	type( t_space ), intent(in)          :: g_space

	! local variables
	real(p_double), dimension(2,rank)    :: g_x_range
	real(p_double), dimension(rank)      :: ldx

	real(p_double)                       :: zmin, z, z_2, dz_2
	real(p_double)                       :: rmin, r, r_2, dr_2

	! other variables
	integer                              :: ikz, ikr, ir, iz, n, in, m, s
	integer, dimension (:), allocatable	 :: n_array 
	real(p_double)                       :: kz, kr, dkz, dkr, omega, aux_boost, phi0
	complex(p_double)                    :: f_k, amp 
	real(p_double)						 :: a0_amp, a0_arg, psi_z, psi_z_2

	! executable statements
	call get_x_bnd( g_space, g_x_range )
	ldx(1) = b_re(1)%dx_(1)
	ldx(2) = b_re(1)%dx_(2)

	! spatial coordinates
	! XMIN ALREADY INCLUDES FOCUS POSITION
	zmin = g_x_range(p_lower,1) + real(nx_p_min(1)-1, p_double)*ldx(1) - this%focus(1)
	rmin = g_x_range(p_lower,2) + real(nx_p_min(2)-1, p_double)*ldx(2)
	dz_2 = ldx(1)/2.0_p_double
	dr_2 = ldx(2)/2.0_p_double

	! store polarization
	phi0 = this%pol

	! get minimum and maximum m
	select case ( this%per_type )
		case ( k_per_gaussian ) 
			allocate ( n_array(1) ) 
			n_array(1) = 0
		case ( k_per_laguerre_gaussian ) 
			allocate ( n_array(1) ) 
			n_array(1) = this%per_tem_mode(2)
	end select

	! testing, debug and control
	if ( control .and. mpi_node() == 0 ) then
		write(0,*) "<control> ", "lon_int_range", this%lon_int_range
		write(0,*) "<control> ", "per1_int_range", this%per1_int_range
	endif

	! iterate for every plane wave
	dkz = ( this%lon_int_range(2)  - this%lon_int_range(1)  ) / (this%num_waves(1) - 1)
	dkr = ( this%per1_int_range(2) - this%per1_int_range(1) ) / (this%num_waves(2) - 1)
	do in = 1, size(n_array)
		! get m mode with respective sign
		n = abs( n_array(in) )
		s = sign( 1, n_array(in) )
	do ikz = 0, this%num_waves(1) - 1 
		! calculate direction of each wave - longitudinal
		kz = this%lon_int_range(1)  + real(ikz) * dkz
	do ikr = 0, this%num_waves(2) - 1
		! calculate direction of each wave - perpendicular
		kr = this%per1_int_range(1) + real(ikr) * dkr
	! --------------------------------------------------

		! calculate omega
		omega = sqrt( kz**2 + kr**2 )

		! boost pulse
		if (this%if_boost) then
			aux_boost = omega
			omega = this%gamma * omega - sqrt(this%gamma**2-1) * kz
			kz    = this%gamma * kz    - sqrt(this%gamma**2-1) * aux_boost
		endif

		! calculate intensity of each wave
		f_k = this%a0 * this%k_function( kz, kr, real(n, p_double) ) * omega * kr * dkz * dkr * s**n
		a0_amp = abs( f_k )
		a0_arg = atan2( dimag( f_k ), real( f_k ) )
	
		! testing, debug and control
		if ( control .and. mpi_node() == 0 ) then
			write(0,*) ikz, ikr
		endif 

		! iterate for every space point
		do iz = b_re(1)%lbound(2), b_re(1)%ubound(2)
			! calculate grid positions and exponentials for z
			z = zmin + real(iz-1, p_double) * ldx(1)
			z_2 = z + dz_2
			psi_z   = s * ( kz * z   - omega * (t + this%lon_position) + a0_arg + this%phase0 )
			psi_z_2 = s * ( kz * z_2 - omega * (t + this%lon_position) + a0_arg + this%phase0 )
		do ir = b_re(1)%lbound(3), b_re(1)%ubound(3)
			! calculate grid positions and exponentials for r
			r = rmin + real(ir-1, p_double) * ldx(2)
			r_2 = r + dr_2
			if (r_2 == 0) r_2 = 1e-6 		! TODO: not the most correct thing to do 
		! --------------------------------------------------

			if ( n == 0 ) then
			! for n = 0 and m = 1
				m = 1
				! e1 at (z_2, r  )
				amp = (a0_amp*kr*bessel_jn(1,kr*r)*sin(psi_z_2))/(exp(im_unit*phi0)*omega)
				e_re(m)%f2(1,iz,ir) = e_re(m)%f2(1,iz,ir) +  real( amp )
				e_im(m)%f2(1,iz,ir) = e_im(m)%f2(1,iz,ir) - dimag( amp )

				! e2 at (z  , r_2)
				amp = (a0_amp*kz*bessel_jn(0,kr*r_2)*cos(psi_z))/(exp(im_unit*phi0)*\
				omega)
				e_re(m)%f2(2,iz,ir) = e_re(m)%f2(2,iz,ir) +  real( amp )
				e_im(m)%f2(2,iz,ir) = e_im(m)%f2(2,iz,ir) - dimag( amp )

				! e3 at (z  , r  )
				amp = (im_unit*a0_amp*kz*bessel_jn(0,kr*r)*cos(psi_z))/(exp(im_unit*phi0)*\
				omega)
				e_re(m)%f2(3,iz,ir) = e_re(m)%f2(3,iz,ir) +  real( amp )
				e_im(m)%f2(3,iz,ir) = e_im(m)%f2(3,iz,ir) - dimag( amp )


				! b1 at (z  , r_2)
				amp = (-im_unit*a0_amp*kz*kr*bessel_jn(1,kr*r_2)*sin(psi_z))/(exp(im_unit*phi0)*omega**2)
				b_re(m)%f2(1,iz,ir) = b_re(m)%f2(1,iz,ir) +  real( amp )
				b_im(m)%f2(1,iz,ir) = b_im(m)%f2(1,iz,ir) - dimag( amp )

				! b2 at (z_2, r  )
				amp = (-im_unit*a0_amp*(kz**2*bessel_jn(0,kr*r) + \
				kr/r*bessel_jn(1,kr*r))*cos(psi_z_2))/(exp(im_unit*phi0)*omega**2)
				b_re(m)%f2(2,iz,ir) = b_re(m)%f2(2,iz,ir) +  real( amp )
				b_im(m)%f2(2,iz,ir) = b_im(m)%f2(2,iz,ir) - dimag( amp )

				! b3 at (z_2, r_2)
				amp = (a0_amp*((2*kz**2 + kr**2)*bessel_jn(0,kr*r_2) - \
				kr**2*bessel_jn(2,kr*r_2))*cos(psi_z_2))/(2.*exp(im_unit*phi0)*\
				omega**2)
				b_re(m)%f2(3,iz,ir) = b_re(m)%f2(3,iz,ir) +  real( amp )
				b_im(m)%f2(3,iz,ir) = b_im(m)%f2(3,iz,ir) - dimag( amp )
			
			else
			! for n > 0 and m = n - 1
				m = n - 1
				if ( m == 0 ) m = this%m_pulse + 1
				! e1 at (z_2, r  )
				amp = (s*im_unit/2*a0_amp*exp(im_unit*(phi0 + psi_z_2))*kr*bessel_jn(-1 + \
				n,kr*r))/omega
				e_re(m)%f2(1,iz,ir) = e_re(m)%f2(1,iz,ir) +  real( amp )
				e_im(m)%f2(1,iz,ir) = e_im(m)%f2(1,iz,ir) - dimag( amp )

				! e2 at (z  , r_2)
				amp = (a0_amp*exp(im_unit*(phi0 + \
				psi_z))*kz*bessel_jn(n,kr*r_2))/(2.*omega)
				e_re(m)%f2(2,iz,ir) = e_re(m)%f2(2,iz,ir) +  real( amp )
				e_im(m)%f2(2,iz,ir) = e_im(m)%f2(2,iz,ir) - dimag( amp )

				! e3 at (z  , r  )
				amp = (-im_unit/2*a0_amp*exp(im_unit*(phi0 + \
				psi_z))*kz*bessel_jn(n,kr*r))/omega
				e_re(m)%f2(3,iz,ir) = e_re(m)%f2(3,iz,ir) +  real( amp )
				e_im(m)%f2(3,iz,ir) = e_im(m)%f2(3,iz,ir) - dimag( amp )


				! b1 at (z  , r_2)
				amp = -0.5*s*(a0_amp*exp(im_unit*(phi0 + psi_z))*kz*kr*bessel_jn(-1 + \
				n,kr*r_2))/omega**2
				b_re(m)%f2(1,iz,ir) = b_re(m)%f2(1,iz,ir) +  real( amp )
				b_im(m)%f2(1,iz,ir) = b_im(m)%f2(1,iz,ir) - dimag( amp )

				! b2 at (z_2, r  )
				amp = (im_unit/2*a0_amp*exp(im_unit*(phi0 + psi_z_2))*(kr/r*(-1 + \
				n)*bessel_jn(-1 + n,kr*r) + kz**2*bessel_jn(n,kr*r)))/(omega**2)
				b_re(m)%f2(2,iz,ir) = b_re(m)%f2(2,iz,ir) +  real( amp )
				b_im(m)%f2(2,iz,ir) = b_im(m)%f2(2,iz,ir) - dimag( amp )

				! b3 at (z_2, r_2)
				amp = -0.25*(a0_amp*exp(im_unit*(phi0 + psi_z_2))*(kr**2*bessel_jn(-2 + \
				n,kr*r_2) - (2*kz**2 + kr**2)*bessel_jn(n,kr*r_2)))/omega**2
				b_re(m)%f2(3,iz,ir) = b_re(m)%f2(3,iz,ir) +  real( amp )
				b_im(m)%f2(3,iz,ir) = b_im(m)%f2(3,iz,ir) - dimag( amp )

			! for n > 0 and m = n + 1
				m = n + 1
				! e1 at (z_2, r  )
				amp = (-0.5*s*im_unit*a0_amp*exp(im_unit*(-phi0 + psi_z_2))*kr*bessel_jn(1 + \
				n,kr*r))/omega
				e_re(m)%f2(1,iz,ir) = e_re(m)%f2(1,iz,ir) +  real( amp )
				e_im(m)%f2(1,iz,ir) = e_im(m)%f2(1,iz,ir) - dimag( amp )

				! e2 at (z  , r_2)
				amp = (a0_amp*exp(im_unit*(-phi0 + \
				psi_z))*kz*bessel_jn(n,kr*r_2))/(2.*omega)
				e_re(m)%f2(2,iz,ir) = e_re(m)%f2(2,iz,ir) +  real( amp )
				e_im(m)%f2(2,iz,ir) = e_im(m)%f2(2,iz,ir) - dimag( amp )

				! e3 at (z  , r  )
				amp = (0.5*im_unit*a0_amp*exp(im_unit*(-phi0 + \
				psi_z))*kz*bessel_jn(n,kr*r))/omega
				e_re(m)%f2(3,iz,ir) = e_re(m)%f2(3,iz,ir) +  real( amp )
				e_im(m)%f2(3,iz,ir) = e_im(m)%f2(3,iz,ir) - dimag( amp )


				! b1 at (z  , r_2)
				amp = -0.5*s*(a0_amp*exp(im_unit*(-phi0 + psi_z))*kz*kr*bessel_jn(1 + \
				n,kr*r_2))/omega**2
				b_re(m)%f2(1,iz,ir) = b_re(m)%f2(1,iz,ir) +  real( amp )
				b_im(m)%f2(1,iz,ir) = b_im(m)%f2(1,iz,ir) - dimag( amp )

				! b2 at (z_2, r  )
				amp = (-im_unit/2*a0_amp*exp(im_unit*(-phi0 + \
				psi_z_2))*(kz**2*bessel_jn(n,kr*r) + kr/r*(1 + n)*bessel_jn(1 + \
				n,kr*r)))/(omega**2)
				b_re(m)%f2(2,iz,ir) = b_re(m)%f2(2,iz,ir) +  real( amp )
				b_im(m)%f2(2,iz,ir) = b_im(m)%f2(2,iz,ir) - dimag( amp )

				! b3 at (z_2, r_2)
				amp = (a0_amp*exp(im_unit*(-phi0 + psi_z_2))*((2*kz**2 + \
				kr**2)*bessel_jn(n,kr*r_2) - kr**2*bessel_jn(2 + \
				n,kr*r_2)))/(4.*omega**2)
				b_re(m)%f2(3,iz,ir) = b_re(m)%f2(3,iz,ir) +  real( amp )
				b_im(m)%f2(3,iz,ir) = b_im(m)%f2(3,iz,ir) - dimag( amp )

			
			endif

		enddo
		enddo

	enddo
	enddo
	enddo

	deallocate (n_array)  

end subroutine launch_zpulse_planewaves_cyl_modes_normal
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Launch pulse routine
subroutine launch_zpulse_planewaves_cyl_modes( this, emf, g_space, nx_p_min, g_nx, t, dt, no_co )

	implicit none

	class( t_zpulse_planewaves_cyl_modes ), intent(inout) :: this
	class( t_emf ), intent(inout) :: emf
	type( t_space ), intent(in) :: g_space
	integer, intent(in), dimension(:) :: nx_p_min, g_nx
	real(p_double), intent(in) :: t
	real(p_double), intent(in) :: dt
	class( t_node_conf ), intent(in) :: no_co

	! field variables
	type( t_vdf ), dimension(this%m_pulse+1) :: e_pulse_re, b_pulse_re, e_pulse_im, b_pulse_im

	! extra variables
	integer :: m

	select type( emf )
	class is( t_emf_cyl_modes )

		! Check if simulation has enough m modes for pulse
		if( ubound( emf%b_cyl_m%pf_re, 1 ) < this%m_pulse ) then 
			if ( mpi_node() == 0 ) then
				write(0,*)  ""
				write(0,*)  "   ERROR: num cyl for simulation modes must be >= pulse modes"
				write(0,*)  "   n_modes: ", ubound( emf%b_cyl_m%pf_re, 1 )
				write(0,*)  "   m_pulse: ", this%m_pulse 
			endif
			stop
		endif

		! Launch
		if ( this%if_launch .and. ( t >= this%launch_time ) ) then

			! create new vdfs to hold the pulse field
			do m = 1, this%m_pulse+1
				call e_pulse_re(m) % new( emf%e, zero = .true. )
				call b_pulse_re(m) % new( emf%b, zero = .true. )
				call e_pulse_im(m) % new( emf%e, zero = .true. )
				call b_pulse_im(m) % new( emf%b, zero = .true. )
			enddo

			call launch_zpulse_planewaves_cyl_modes_normal( this, e_pulse_re, b_pulse_re, e_pulse_im, b_pulse_im, t, nx_p_min, g_space )

			! take care of circular polarisation
			if (this%pol_type /= 0) then

				this%pol    = this%pol    + real( pi_2, p_double )
				this%phase0 = this%phase0 + real( sign( pi_2, real( this%pol_type, p_double )) , p_double )
				
				call launch_zpulse_planewaves_cyl_modes_normal( this, e_pulse_re, b_pulse_re, e_pulse_im, b_pulse_im, t, nx_p_min, g_space )

				this%pol    = this%pol    - real( pi_2, p_double )
				this%phase0 = this%phase0 - real( sign( pi_2, real( this%pol_type, p_double )) , p_double )

			endif
			
			! add laser pulse to fields
			do m = 1, this%m_pulse
				call add( emf%e_cyl_m%pf_re(m), e_pulse_re(m) )
				call add( emf%b_cyl_m%pf_re(m), b_pulse_re(m) )
				call add( emf%e_cyl_m%pf_im(m), e_pulse_im(m) )
				call add( emf%b_cyl_m%pf_im(m), b_pulse_im(m) )
			enddo
			! add component 0
			call add( emf%e_cyl_m%pf_re(0), e_pulse_re(this%m_pulse+1) )
			call add( emf%b_cyl_m%pf_re(0), b_pulse_re(this%m_pulse+1) )

			! free b_pulse and e_pulse vdfs
			do m = 1, this%m_pulse+1
				call e_pulse_re(m) % cleanup()
				call b_pulse_re(m) % cleanup()
				call e_pulse_im(m) % cleanup()
				call b_pulse_im(m) % cleanup()
			enddo

			! pulse has been launched, turn off if_launch
			this%if_launch = .false.

		endif 

	end select

end subroutine launch_zpulse_planewaves_cyl_modes
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Wave vector weight function
function k_function_zpulse_planewaves_cyl_modes( this, k_lon, k_per1, k_per2 ) result( f_k )

	implicit none
	class( t_zpulse_planewaves_cyl_modes ), intent(in)   :: this
	real(p_double), intent(in)                 :: k_lon, k_per1, k_per2
	real(p_double)			                   :: kz, kr, m
	complex(p_double)                          :: f_k

	! auxiliar variables
	real(p_double) 							   :: f_k_lon_arg, vg

	! parse data to new variables
	kz = k_lon 
	kr = k_per1
	m = k_per2

	if ( this%abouraddy_angle == pi / 2 ) then
		f_k_lon_arg = kz
	else 
		vg = tan( this%abouraddy_angle )
		f_k_lon_arg = ( -kz * vg + sqrt( kz**2 + kr**2 ) ) / ( 1 - vg )
	endif
	f_k = this%k_lon_function( f_k_lon_arg ) * this%k_per_function( kr, m ) 

end function k_function_zpulse_planewaves_cyl_modes
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
! Transverse wave vector weight function
function k_per_function_zpulse_planewaves_cyl_modes( this, k_per1, k_per2 ) result( f_k_per )

	implicit none
	class( t_zpulse_planewaves_cyl_modes ), intent(in)   :: this
	real(p_double), intent(in)                 :: k_per1, k_per2
	real(p_double)			                   :: kr
	integer									   :: m
	complex(p_double)                          :: f_k_per

	! auxiliar variables
	real(p_double)                             :: w0
	integer 								   :: p, l

	! parse data to new variables
	kr = k_per1
	m = int(k_per2)

	select case ( this%per_type )
		case (k_per_gaussian)

			if ( m == 0 ) then
				f_k_per = this%per_w0(1)**2 / 2.0 * exp( - kr**2*this%per_w0(1)**2 / 4 )
			else 
				f_k_per = 0.0_p_double
			endif

		case (k_per_laguerre_gaussian)

			p = this%per_tem_mode(1)
			l = this%per_tem_mode(2)
			w0 = this%per_w0(1)

			if ( m == abs(l) ) then
				!f_k_per = w0**2 * (-1)**p / 2
				!f_k_per = w0**2 * (-im_unit)**l * (-1)**(p)

				f_k_per = w0**2 * (-1)**(p)
				f_k_per = f_k_per * exp( -kr**2*w0**2/4 )
				f_k_per = f_k_per * (kr*w0/sqrt(2.0_p_double))**abs(l)
				f_k_per = f_k_per * laguerre( p, abs(l), kr**2*w0**2/2 )

				f_k_per = f_k_per / this%normalization

				if ( l < 0 ) f_k_per = f_k_per - 2 * im_unit * dimag(f_k_per)

			else 
				f_k_per = 0.0_p_double
			endif

	end select

end function k_per_function_zpulse_planewaves_cyl_modes
!-----------------------------------------------------------------------------------------


end module m_zpulse_planewaves_cyl_modes