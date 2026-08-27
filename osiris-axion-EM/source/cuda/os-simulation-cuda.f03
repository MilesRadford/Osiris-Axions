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

subroutine grow_buffer_spec_event(n,size)
  !use iso_c_binding

  implicit none

  integer :: n,size

  !C_FUNLOC
  ! call c_func_with_cb_to_fortran( 5, 7, c_funloc ( callback_in_fortran ) )
  !
  ! Callback for 'C' to use in order to communicate
  !function callback_in_fortran (arg) bind(c)
  !integer(c_int), intent(in), value :: arg
  !integer(c_int) :: callback_in_fortran
  !
  !callback_in_fortran = arg + arg
  !end function callback_in_fortran

  !void c_func_with_cb_to_fortran (int i, int*j, void* func_ptr) {
  !    int tempi, tempj;
  !    int (*func)(int);
  !
  !    func = ( int(*)(int) ) func_ptr;
  !    tempi = func( i );
  !    tempj = func( *j );
  !
  !  printf("C called Fortran and it spit back i=%d j=%d \n", i , *j);
  !  printf("C called Fortran and it spit back i=%d j=%d \n", tempi , tempj);
  !}

  print *, "!!!!!!!!!Cgrow_buffer_spec_event EVENT!!! "
  print *, "grow_buffer_spec_event EVENT!!! ", n,size

end subroutine grow_buffer_spec_event

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
! Initialize simulation objects
!-----------------------------------------------------------------------------------------
subroutine init_cuda( sim )

  use m_grid_parallel
  use m_random
  use m_logprof
  use m_particles
  use m_simulation_cuda_define
  use m_species_define
  use m_particles_define

  use m_walk_data_structures
  
  implicit none

  interface
    subroutine grow_buffer_spec_event(n,size)
      integer :: n,size
    end subroutine
  end interface

  interface
    subroutine init_spec_buffer_event(n,size)
      integer :: n,size
    end subroutine init_spec_buffer_event
  end interface

  interface
    subroutine inject_area_spec_event(n,size)
      integer :: n,size
    end subroutine inject_area_spec_event
  end interface

  ! dummy variables

  class ( t_simulation_cuda ), intent(inout) ::  sim
  integer :: k

  ! setup the function pointers to 'listen' to respond to various particle events
  do k=1, sim%part%num_species
    sim%part%species(k)%events%grow_buffer_spec_event => grow_buffer_spec_event
    sim%part%species(k)%events%init_spec_buffer_event => init_spec_buffer_event
    sim%part%species(k)%events%inject_area_spec_event => inject_area_spec_event
  end do

  ! call super class
  call sim % t_simulation % init()

  ! call 'setup_cuda' which resides in os-cuda-serialize-state.f90
  !   It's first mian function is to serialize the state and create C/C++ structure
  !   that mirrors Osiris's Fortran description objects.
  call setup_cuda( sim )  

end subroutine init_cuda

