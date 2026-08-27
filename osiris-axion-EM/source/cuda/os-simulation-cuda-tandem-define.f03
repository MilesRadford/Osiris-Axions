!#define DEBUG_FILE 1

module m_simulation_cuda_tandem_define

#include "os-config.h"
#include "os-preprocess.fpp"

use m_system
use m_parameters

use m_node_conf
use m_grid_define
use m_grid
use m_grid_parallel

use m_time_step
use m_restart
use m_system

use m_space
use m_time
use m_emf_define
use m_particles_define

use m_zpulse
use m_current_define
use m_antenna_array

use m_input_file
use m_simulation
use m_random
use m_cosimulation
use m_handoff

implicit none

private

! This should and needs to go away.. 
!   each thing listed here should be tied back to specific simulation instances.
type t_global_state
  class(t_random), pointer, public :: rng => null()
end type t_global_state

type, extends( t_simulation ) :: t_simulation_cuda_tandem
  
  class( t_simulation ),pointer :: sim1 => null()
  class( t_simulation ),pointer :: sim2 => null()

  type ( t_global_state ) :: sim1_global_state
  type ( t_global_state ) :: sim2_global_state

contains

   procedure :: iter               => iter_sim_tandem

   procedure :: allocate_objs      => allocate_objs_sim_tandem
   procedure :: init               => init_sim_tandem
   procedure :: cleanup            => cleanup_sim_tandem
   procedure :: read_input         => read_input_sim_tandem

   procedure :: write_checkpoint   => write_checkpoint_sim_tandem
   procedure :: list_algorithm     => list_algorithm_sim_tandem
   procedure :: test_input         => test_input_sim_tandem

   procedure :: set_int_load       => set_int_load_sim_tandem

   procedure :: set_vars_to_set_1
   procedure :: set_vars_to_set_2

end type t_simulation_cuda_tandem


class( t_simulation_cuda_tandem ), pointer :: the_sim

character(len=50) :: context_return_value1
character(len=50) :: context_return_value2
character(len=50) :: context_return_temp


! exported symbols
public :: t_simulation_cuda_tandem, the_sim

contains

!-----------------------------------------------------------------------------------------
! Read simulation parameters from input file
!-----------------------------------------------------------------------------------------
subroutine read_input_sim_tandem( sim, input_file )

  use m_particles

  implicit none

  class ( t_simulation_cuda_tandem ), intent(inout) :: sim
  type( t_input_file ), intent(inout) :: input_file

  type( t_input_file ) :: input_file2
  integer :: ierr

  ! propagate the 'options' to the sub-objects
  sim%sim1%options = sim%options
  sim%sim2%options = sim%options

  ! 
  !   must be read in order.. so create a 2nd 
  ! The incoming 'input file' we shall say goes with SIM1. So create a duplicate to satify the initization of SIM2
  call open_input( input_file2, trim( sim%sim2%options%input_file ), mpi_world_comm() )
  ! One hitch: the incoming 'input file' had already had the 'simualtion' section read. As of now, these things
  !   namelist sections must all be read and in order, so pretend to read that section in for our new 'input file'
  call get_namelist( input_file2, "nl_simulation", ierr )

  SCR_ROOT("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
  SCR_ROOT("%%%%%%%%%%%%%% READING INPUT SIM 1 %%%%%%%%%%%%%%")
  call sim%sim1%read_input( input_file  )
  SCR_ROOT("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
  SCR_ROOT("%%%%%%%%%%%%%% READING INPUT SIM 2 %%%%%%%%%%%%%%")
  call sim%sim2%read_input( input_file2  )
  SCR_ROOT("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
  call sys_chdir('sim1', ierr)

end subroutine read_input_sim_tandem


!-----------------------------------------------------------------------------------------
! Initialize simulation objects
!-----------------------------------------------------------------------------------------
subroutine allocate_objs_sim_tandem( sim )

  implicit none

  interface
    function global_simulation_factory( algorithm  )
      import t_simulation
      implicit none

      integer, intent(in) :: algorithm

      class( t_simulation ), pointer :: global_simulation_factory
    end function
  end interface 

  class( t_simulation_cuda_tandem ), intent(inout) :: sim

  class( t_simulation ), pointer :: ptr

  ! hack: 'global_simulation_factory' does not even make use of it's input_file
  !       paramter, so just pass in a junk one.
  ! TODO: take input_file out of the args for 'global_simulation_factory'
  sim%options%algorithm = p_standard
  sim%sim1  => global_simulation_factory( p_standard )

  ! propagate the 'options' to the sub-objects
  !sim%sim1%options = sim%options
  !sim%sim2%options = sim%options

  sim%options%algorithm = p_sim_cuda
  sim%sim2 => global_simulation_factory( p_sim_cuda )

  ! propagate the 'options' to the sub-objects
  !sim%sim1%options = sim%options
  !sim%sim2%options = sim%options

  call sim%sim1%allocate_objs( )
  call sim%sim2%allocate_objs( )

  call sim%set_vars_to_set_1()

end subroutine allocate_objs_sim_tandem

!-----------------------------------------------------------------------------------------
! Initialize simulation objects
!-----------------------------------------------------------------------------------------

subroutine save_sim_local( sim )
  implicit none
  class ( t_simulation_cuda_tandem ), target, intent(inout) ::  sim
  the_sim => sim
end subroutine


subroutine init_sim_tandem( sim )

  use m_grid_parallel
  use m_random
  use m_logprof
  use m_particles
  use m_emf
  use m_diagnostic_utilities
  use m_node_conf 
  use m_species_define

  implicit none

  ! dummy variables

  class ( t_simulation_cuda_tandem ), intent(inout) ::  sim
  integer :: i
  class( t_species ), pointer                 :: species

  SCR_ROOT("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
  SCR_ROOT("%%%%%%%%%%%%%% INIT SIM 1 %%%%%%%%%%%%%%")
  call sim%sim1%init( )
  sim%sim1_global_state%rng => rng

  SCR_ROOT("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
  SCR_ROOT("%%%%%%%%%%%%%% INIT SIM 2 %%%%%%%%%%%%%%")
  call sim%sim2%init( )
  sim%sim2_global_state%rng => rng
  SCR_ROOT("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

  !sim%tstep = sim%sim1%tstep
  !sim%time = sim%sim1%time
  call sim%set_vars_to_set_1()

  !print *, root(sim%sim1%no_co)
  !print *, root(sim%sim2%no_co)
  !print *, root(sim%no_co)
  call save_sim_local( sim )
  
  ! ADAM COROUTINE 
  call setup_cosimluation_magic( sim )

  do i=1, sim%part%num_species
      species => sim%sim2%part%get_spec_by_idx( i )
      species%n_sort = 0
  enddo



end subroutine init_sim_tandem

!-----------------------------------------------------------------------------------------
! Cleanup simulation objects
!-----------------------------------------------------------------------------------------
subroutine cleanup_sim_tandem( sim )

  use m_vdf_comm
  use m_wall_comm
  use m_particles
  use m_diagnostic_utilities

  implicit none

  class( t_simulation_cuda_tandem ), intent(inout) :: sim

  call sim%sim1%cleanup( )
  call sim%sim2%cleanup( )

end subroutine cleanup_sim_tandem

!-----------------------------------------------------------------------------------------
! Simulation iteration
!-----------------------------------------------------------------------------------------
subroutine iter_sim_tandem_old( sim )
  use m_emf
  use m_emf_marder
  use m_particles
  use m_logprof
  use m_particles_charge
  use m_space
  use m_vdf_define

  implicit none

  class ( t_simulation_cuda_tandem ), intent(inout) :: sim

  ! coming into 'iter_sim_tandem', we have set the active variables to SIM1 (in the init method)
  call sim1_iter() !call sim%sim1%iter( )
  call sim%set_vars_to_set_2()
  
  call sim2_iter()  !call sim%sim2%iter( )

  ! set back to SIM1 for the next iteration of the loop.
  call sim%set_vars_to_set_1()


end subroutine



subroutine do_some_stuff() 

  use m_time_step
  use m_species_define
  
  implicit none
  integer :: i, str_idx, int_tag_value, j, k, d
  integer :: intval1, intval2, n_step
  real(p_k_fld) :: temp, tol, tt, tt1, tt2
  class(t_species), pointer :: species

  n_step =  n(the_sim%tstep)
  !if( n(the_sim%tstep) < 3290) then
  !  return
  !endif

  tol = 1e-5
  if(  context_return_value1 /= context_return_value2 ) then
    return
  endif
  
  str_idx = len(context_return_value1); int_tag_value = 0;
  do while ( int_tag_value > -1 )
    int_tag_value = get_trailing_int_in_str( context_return_value1,str_idx)  
  end do
  
  context_return_temp = context_return_value1(1:str_idx)
 
#if 0
  select case( trim( context_return_temp ) )


    !case( "PUSH_START" )
    !  !intval1 = get_trailing_int_in_str( context_return_value1 )
    !  if(intval1 > 0) then
    !    print *, "     PUSH_START on species ", intval1
    !    ! here, if needed, we copy all the particales from the CPU HOST to the GPU
    !    call COPY_PARTICLES_TO_GPU( intval1 );
    !  endif

    case( "AFTER_EMF_SOLVER" )
      print *, "      DOING FIELD COMPARE!!!!"
      if( .not. the_sim%sim1%emf%e%equals( the_sim%sim2%emf%e,tol ) ) then
        
        do i=1,size(the_sim%sim1%emf%e % buffer)
          temp = ABS( the_sim%sim1%emf%e%buffer(i) - the_sim%sim2%emf%e%buffer(i) )
          if(temp > tol) then
            print *, "             ", temp 
          end if
        enddo
        stop
      end if

    case("before_update_boundary_particles" )
      print *, "    Doing the CPU sort only for testing!"
      ! jay is already copied over..
      !species => the_sim%sim2%part%species
      !print *, "         copying particles from to gpu to cpu"
      !call SYNCH_PARTICLE_BUFFERS_CUDA( 1 )

    case("PUSH_START" )
      tt = 0.0
      do i=1,size(the_sim%sim1%emf%e % buffer)
        temp = ABS( the_sim%sim1%emf%e%buffer(i) - the_sim%sim2%emf%e%buffer(i) )
        if(temp > tt) then
          tt = temp
        end if
      enddo
       print *, "    Comparing Fields. Max error:", tt

    ! Test to see if Jay (the main point of all this) is in good shape!
    case( "AFTER_advance_deposit" )
      print *, "    Valiodating final Current!"
      !call VALIDATE_DEPOSIT(  the_sim%sim2%part%jay_scratch%f2 )
      tt = 0.0
      tt1 = 0.0
      tt2 = 0.0

      !if( .not. the_sim%sim1%jay%pf(1)%equals( the_sim%sim2%jay%pf(1),1e-2 ) ) then
        do j=lbound( the_sim%sim1%jay%pf(1)%f2, 3), ubound( the_sim%sim1%jay%pf(1)%f2, 3)
          do i=lbound( the_sim%sim1%jay%pf(1)%f2, 2), ubound( the_sim%sim1%jay%pf(1)%f2, 2)
            do d=1,3
              temp = ABS( the_sim%sim1%jay%pf(1)%f2(d,i,j) - the_sim%sim2%jay%pf(1)%f2(d,i,j ) )
              if(temp > tt) then
                tt = temp
              endif
              if(  the_sim%sim1%jay%pf(1)%f2(d,i,j) > tt1 ) then
                 tt1 = the_sim%sim1%jay%pf(1)%f2(d,i,j)
              endif
              if(  the_sim%sim2%jay%pf(1)%f2(d,i,j) > tt2 ) then
                 tt2 = the_sim%sim2%jay%pf(1)%f2(d,i,j)
              endif

              !if(temp > tol) then
              !  print *, i, j, d, "     JAY Over Tolerance!: ", temp
              !endif
            enddo
          enddo
        enddo
        !print *, tt
      !endif
      print *, "    Valiodating final Current. Max error:", tt
      print *, "    ", tt1, tt2

      ! if( .not. the_sim%sim1%jay%pf(1)%equals( the_sim%sim2%jay%pf(1),tol ) ) then
      !   do i=1,size(the_sim%sim1%jay%pf(1) % buffer)
      !     temp = ABS( the_sim%sim1%jay%pf(1)%buffer(i) - the_sim%sim2%jay%pf(1)%buffer(i) )
      !     if(temp > tol) then
            
      !     end if
      !   enddo
      !   stop
      ! endif
    case default
      !print *,"    Event ",  trim( context_return_temp ), ' has no action.'
  end select
#endif

end subroutine


subroutine iter_sim_tandem( sim )

  use m_emf
  use m_emf_marder
  use m_particles
  use m_logprof
  use m_particles_charge
  use m_space
  use m_vdf_define

  implicit none

  class ( t_simulation_cuda_tandem ), intent(inout) :: sim

  ! Quick debugging test:
  !   Change the particles data for sim1 and see if we get different diagnostic results 
  !     between sim1 and sim2
  ! if( n(sim%tstep) .eq. 0) then
  !   SCR_ROOT("%%%%%%%%%%%%%% FIRST TIME STEP %%%%%%%%%%%%%%%")
  !   sim%sim1%part%species%p(1,1) = 4.0
  ! endif
  
  call start_sim_itteration()

  context_return_value1 = ''
  context_return_value2 = ''

  do while( context_return_value1 .ne. '__END'  .AND. context_return_value2 .ne. '__END')

    ! coming into 'iter_sim_tandem', we have set the active variables to SIM1 (in the init method)
    !print *, "********************> Switch to SIM 1 **************"
    call execute_sim1_iter()
    context_return_value1 = context_return_value
    !print *, "********************> SIM 1 paused as ", context_return_value
    !print *, "     (SIM 1 gave context:", trim(context_return_value1), ")"

    call sim%set_vars_to_set_2()

    !print *, "********************> Switch to SIM 2 **************"
    call execute_sim2_iter()
    context_return_value2 = context_return_value
    !print *, "********************> SIM 2 paused as ", context_return_value
    !print *, "     (SIM 2 gave context:", trim(context_return_value2), ")"

    ! set back to SIM1 for the next iteration of the loop.
    call sim%set_vars_to_set_1()

    call do_some_stuff( )
  end do

end subroutine iter_sim_tandem

!-----------------------------------------------------------------------------------------
! Get the minimum number of guard cells in all directions / boundaries
!-----------------------------------------------------------------------------------------
subroutine list_algorithm_sim_tandem( sim )

  use m_particles
  use m_emf
  use m_random

  implicit none

  class( t_simulation_cuda_tandem ), intent(in) :: sim

  SCR_ROOT("%%%%%%%%%%%%%% FINAL SETUP OF SIM 1 %%%%%%%%%%%%%%")
  call sim%sim1%list_algorithm( )
  SCR_ROOT("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
  SCR_ROOT("%%%%%%%%%%%%%% FINAL SETUP OF SIM 2 %%%%%%%%%%%%%%")
  call sim%sim2%list_algorithm( )
  SCR_ROOT("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

end subroutine list_algorithm_sim_tandem


!---------------------------------------------------------------------------------------
! Write checkpoint data
!---------------------------------------------------------------------------------------
subroutine cheat_set_var_to_1( sim )
  implicit none
  class ( t_simulation_cuda_tandem ), intent(inout) ::  sim
  call sim%set_vars_to_set_1()
end subroutine

subroutine write_checkpoint_sim_tandem( sim, restart_handle )

  use m_random
  use m_particles

  implicit none

  class ( t_simulation_cuda_tandem ), intent(in) ::  sim
  type ( t_restart_handle ), intent(inout)    ::  restart_handle

  SCR_ROOT("%%%%%%%%%%%%%% RESTARTS NOT IMPLEMENTED/TESTED YET %%%%%%%%%%%%%%")

  SCR_ROOT("---> Writing restart data for SIM 1")
  call sim%sim1%write_checkpoint( restart_handle )
  call the_sim%set_vars_to_set_2()
  SCR_ROOT("---> Writing restart data for SIM 2")
  call sim%sim2%write_checkpoint( restart_handle )
  call the_sim%set_vars_to_set_1()
  SCR_ROOT("---> Finished writing Restart data.")

end subroutine write_checkpoint_sim_tandem

!-----------------------------------------------------------------------------------------
! Check the input information for errors
!-----------------------------------------------------------------------------------------
subroutine test_input_sim_tandem( sim )

  use m_emf
  use m_grid

  implicit none

  class( t_simulation_cuda_tandem ), intent(inout) :: sim

  call sim%sim1%test_input( )
  call sim%sim2%test_input( )

end subroutine test_input_sim_tandem

!-----------------------------------------------------------------------------------------
! Sets the integral of load on all directions for node partition / dynamic load balance
! - The result is saved on a grid object, which may not be the main grid object
!-----------------------------------------------------------------------------------------
subroutine set_int_load_sim_tandem( sim, grid, n )

  use m_grid_parallel
  use m_particles
  use m_logprof

  implicit none

  class( t_simulation_cuda_tandem ), intent(inout) :: sim
  class( t_grid ), intent(inout) :: grid
  integer, intent(in) :: n

  call sim%sim1%set_int_load( grid, n )
  call sim%sim2%set_int_load( grid, n )

end subroutine set_int_load_sim_tandem
!-----------------------------------------------------------------------------------------

subroutine set_vars_to_set_1( sim )
  implicit none
  class( t_simulation_cuda_tandem ), intent(inout) :: sim
  integer :: ierr

  sim%antenna_array = sim%sim1%antenna_array
  sim%restart = sim%sim1%restart
  sim%tstep = sim%sim1%tstep
  sim%time = sim%sim1%time
  !sim%sim1%tstep = sim%tstep
  !sim%sim1%time = sim%time

  sim%grid => sim%sim1%grid
  sim%no_co => sim%sim1%no_co
  sim%g_space => sim%sim1%g_space
  sim%jay => sim%sim1%jay
  sim%emf => sim%sim1%emf
  sim%part => sim%sim1%part

  rng => sim%sim1_global_state%rng

  call sys_chdir('../sim1', ierr)
end subroutine set_vars_to_set_1

subroutine set_vars_to_set_2( sim )
  implicit none
  class( t_simulation_cuda_tandem ), intent(inout) :: sim
  integer :: ierr
  sim%antenna_array = sim%sim2%antenna_array
  sim%restart = sim%sim2%restart
  sim%tstep = sim%sim2%tstep
  sim%time = sim%sim2%time
  !sim%sim2%tstep = sim%tstep
  !sim%sim2%time = sim%time

  sim%grid => sim%sim2%grid
  sim%no_co => sim%sim2%no_co
  sim%g_space => sim%sim2%g_space
  sim%jay => sim%sim2%jay
  sim%emf => sim%sim2%emf
  sim%part => sim%sim2%part

  rng => sim%sim2_global_state%rng

  call sys_chdir('../sim2', ierr)

end subroutine set_vars_to_set_2

end module m_simulation_cuda_tandem_define


!
! These functions are the routines called by the backing C++ Coroutine code
!   upon switching contexts.. The C++ has no idea what these routines do.. they are blind calls.
!
! TODO: these functions are here only becasue I am too lazy right now to get the proper function pointer
!  for use by C++ directly.
subroutine sim1_iter()
  use m_simulation_cuda_tandem_define
  implicit none
  call the_sim%sim1%iter( )
end subroutine

subroutine sim2_iter()
  use m_simulation_cuda_tandem_define
  implicit none
  call the_sim%sim2%iter( )
end subroutine
