!   Basic Parts:
!     Class overides to call cuda things/disable selectable via Osiris things using configuration state.
!
!     Events that notify Cuda:   FORTRAN-C intermitant calls
!     Events C To Fortran
!
!     Workshorse Cuda stuff called from FORTRAN
!     Serialazation of basic Osiris FORTRAN structures into C (setup cuda machine state)
!     Integration of Diagnosics
!
!     How to get stuuf like Hybrid solver to work.
!
module m_simulation_cuda_define

#include "os-config.h"
#include "os-preprocess.fpp"

use m_simulation
use m_particles_define

implicit none
private

type, extends( t_emf ) :: t_emf_cuda
  contains
    procedure :: update_boundary    => update_boundary_emf_cuda

end type t_emf_cuda

type, extends( t_particles ) :: t_particles_cuda
  contains
    procedure :: update_boundary    => update_boundary_particles_cuda

end type t_particles_cuda

type, extends( t_current ) :: t_current_cuda
  contains
    procedure :: update_boundary    => update_boundary_current_cuda
end type t_current_cuda

type, extends( t_simulation ) :: t_simulation_cuda
    ! No member-specific data
  contains
    procedure :: allocate_objs      => allocate_objs_cuda
    procedure :: read_input         => read_input_cuda
    procedure :: list_algorithm     => list_algorithm_cuda
    procedure :: test_input         => test_input_cuda
    procedure :: init               => init_cuda
end type t_simulation_cuda

! a pointer to the cuda simualtion state to be used within this module
type(t_simulation_cuda),pointer ::  cuda_sim_state

! actual code is in os-simulation-cuda.f03 to prevent circular dependnecies.
interface
  subroutine init_cuda( sim )
    import t_simulation_cuda
    class(t_simulation_cuda), intent(inout) ::  sim
  end subroutine
end interface

public :: t_simulation_cuda, t_emf_cuda, t_particles_cuda, t_current_cuda, cuda_sim_state

contains

!*******************************************************************
!**************** For overriding Simualtion methods ****************
!-----------------------------------------------------------------------------------------
! save a copy of the incoming simulation object for use within this module.
!-----------------------------------------------------------------------------------------
subroutine set_local_cuda_state(sim)
  implicit none
  class( t_simulation_cuda ), intent(inout),target :: sim

  ! this could be done in the one of the overidden 'simulation' methods, BUT Fortran
  !   being Fortran requires the 'target' attribute on anything a pointer points to...
  !   and this attribute throws off the polymorphic overriding of the methods. Thus is small
  !   stub function is here to add the 'target' attribute to the sim object and make a pointer to it
  cuda_sim_state=>sim
end subroutine set_local_cuda_state

!-----------------------------------------------------------------------------------------
! Initialize simulation objects
!-----------------------------------------------------------------------------------------
subroutine allocate_objs_cuda( sim )

  !use m_emf_pgc
  !use m_particles_pgc
  use iso_c_binding
  implicit none

  class( t_simulation_cuda ), intent(inout) :: sim

  allocate( t_emf_cuda :: sim%emf )
  allocate( t_particles_cuda :: sim%part )
  allocate( t_current_cuda :: sim%jay )

  ! allocate any remaining objects in superclass
  call sim % t_simulation % allocate_objs()


  !cuda_sim_state = sim
  call set_local_cuda_state(sim)
  !SCR_ROOT('Allocating Cuda emf and particle objects.')

end subroutine allocate_objs_cuda


!-----------------------------------------------------------------------------------------
! Read simulation parameters from input file
!-----------------------------------------------------------------------------------------
subroutine read_input_cuda( sim, input_file )

  use m_particles

  implicit none

  class ( t_simulation_cuda ), intent(inout) :: sim
  class( t_input_file ), intent(inout) :: input_file

  !SCR_ROOT('Reading parallel node configuration... ')
  !call read_nml(  sim%no_co, input_file, p_x_dim, sim%options%algorithm  )

  ! call super class
  call sim % t_simulation % read_input(input_file)
  !SCR_ROOT('CUDA: read_input_sim')

end subroutine read_input_cuda

subroutine cool_callback1(n,size)
  integer :: n,size
  print *, "EVENT!!! "
  print *, "EVENT!!! ", n,size

end subroutine cool_callback1

subroutine init_spec_buffer_event(n,size)
  integer :: n,size
  print *, "init_spec_buffer_event EVENT!!! "
  print *, "init_spec_buffer_event EVENT!!! ", n,size

end subroutine init_spec_buffer_event

subroutine inject_area_spec_event(n,size)
  integer :: n,size
  print *, "inject_area_spec_event EVENT!!! "
  print *, "inject_area_spec_event EVENT!!! ", n,size

end subroutine inject_area_spec_event



!-----------------------------------------------------------------------------------------
! Get the minimum number of guard cells in all directions / boundaries
!-----------------------------------------------------------------------------------------
subroutine list_algorithm_cuda( sim )

  use m_particles

  implicit none

  class( t_simulation_cuda ), intent(in) :: sim

  ! call super class
  call sim % t_simulation % list_algorithm( )
  !SCR_ROOT('CUDA: list_algorithm_cuda')

end subroutine list_algorithm_cuda


!-----------------------------------------------------------------------------------------
! Check the input information for errors
!-----------------------------------------------------------------------------------------
subroutine test_input_cuda( sim )

  implicit none

  class( t_simulation_cuda ), intent(inout) :: sim

  ! call super class
  call sim % t_simulation % test_input( )
  !SCR_ROOT('CUDA: test_input_cuda')

end subroutine test_input_cuda






!************************************************************
!**************** For overriding EMF methods ****************
!-----------------------------------------------------------------------------------------
!       update boundaries of the electro-magnetic field
!-----------------------------------------------------------------------------------------
subroutine update_boundary_emf_cuda( this, g_space, no_co )

  use m_emf_define
  use m_space
  use m_node_conf
  use m_emf_bound

  implicit none

  class( t_emf_cuda ),intent( inout )  ::  this
  type( t_space ),     intent(in) :: g_space
  class( t_node_conf ), intent(in) :: no_co


  ! Check that we should use this or the CUDA

  ! update boundaries on emf fields
  !SCR_ROOT("CUDA: START EMF BOUNDY")
  !call update_boundary( this%bnd_con, this%b, this%e,  g_space, no_co )
  call this % t_emf % update_boundary( g_space, no_co )
  !SCR_ROOT("CUDA: END   EMF BOUNDY")

  ! no update is required for external fields since these are constant or calculated at every time
  ! step

end subroutine update_boundary_emf_cuda






!*****************************************************************
!**************** For overriding Patricle methods ****************
!-----------------------------------------------------------------------------------------
! updates the particle data at the boundaries
!-----------------------------------------------------------------------------------------
subroutine update_boundary_particles_cuda( this, jay, g_space, no_co, dt )

  use m_system
  use m_particles_define
  use m_species
  use m_neutral
  use m_vdf_define
  use m_space
  use m_node_conf
  use m_logprof

  implicit none

  class( t_particles_cuda ), intent(inout) :: this
  type( t_vdf ),       intent(inout) :: jay

  type( t_space ),     intent(in) :: g_space
  class( t_node_conf ), intent(in) :: no_co
  real(p_double),     intent(in) :: dt

! local variables

  integer :: species_id, neutral_id

! executable statements
  !SCR_ROOT('START: CUDA Update PArticles Boundry')
  call this % t_particles % update_boundary( jay, g_space, no_co, dt )
  !SCR_ROOT('END:   CUDA Update PArticles Boundry')

end subroutine update_boundary_particles_cuda
!-----------------------------------------------------------------------------------------







!****************************************************************
!**************** For overriding Current methods ****************
!-----------------------------------------------------------------------------------------
! Update boundary of electrical current
!-----------------------------------------------------------------------------------------
subroutine update_boundary_current_cuda( this, no_co )

  use m_vdf_comm
  use m_node_conf

  implicit none

  class( t_current_cuda ), intent(inout) :: this
  class( t_node_conf ),        intent(in) :: no_co

  integer :: i

  !SCR_ROOT('START: CUDA Update Current Boundry')
  call this % t_current % update_boundary( no_co )
  !SCR_ROOT('END:   CUDA Update Current Boundry')

end subroutine update_boundary_current_cuda
!---------------------------------------------------------------------------------------------------


end module m_simulation_cuda_define



!*****************************************************************
!**************** FORTRAN callbacks called from C ****************
!-----------------------------------------------------------------------------------------
! TODO: make these use ISO C Bindings
!-----------------------------------------------------------------------------------------

! Cuda/C code pusher needs to resize (i.e. increase) the particle buffer size because too manyy
!   particles flowed from other nodes via MPI. (The CUDA code currently handles its own MPI..)
!   one possiibility to better integrate so the CUDA stuff uses the basic OSIRIS MPI code.
subroutine c_asks__grow_species_buffer( species_index1, num_par_new1 )
  use m_species_utilities
  use m_particles_define
  use m_particles_define
  use m_simulation_cuda_define, only: cuda_sim_state

  implicit none
  integer, intent(in):: species_index1, num_par_new1

  call grow_buffer( cuda_sim_state%part%species(species_index1), num_par_new1 )

end subroutine c_asks__grow_species_buffer


! Cuda/C code calls this when it needs to inject new particles into the simulation grid
!   This is odne, for example, when using moving window.
subroutine inject_particles_into_area ( species_index, ig_xbnd_inj, g_nx )

  use m_parameters
  use m_species_utilities
  use m_particles_define
  use m_particles_define
  use m_simulation_cuda_define, only: cuda_sim_state

  implicit none

  integer, intent(in):: species_index
  integer, dimension(3, p_x_dim),   intent(in):: ig_xbnd_inj
  integer, dimension(p_x_dim),   intent(in):: g_nx

  !call inject_area( this, ig_xbnd_inj, get_g_nx(grid) )

  !call inject_area( cuda_sim_state%part%species(species_index), ig_xbnd_inj,g_nx )

end subroutine inject_particles_into_area
