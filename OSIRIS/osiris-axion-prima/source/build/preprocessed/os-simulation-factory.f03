# 0 "os-simulation-factory.f03"
# 0 "<built-in>"
# 0 "<command-line>"
# 1 "os-simulation-factory.f03"
!
! Creating a new simluation type involves:
!
! 1. a 'use' statement to include the new code
! 2. a new entry in the 'case' statement within the 'create_core_simulation' function.
!
module m_simulation_factory

  use m_system
  use m_simulation






    use m_simulation_qed
# 45 "os-simulation-factory.f03"
    use m_simulation_neutral_spin



    use m_simulation_xxfel
# 61 "os-simulation-factory.f03"
    use m_simulation_gr






  implicit none

  private

  interface create
    module procedure create_simulation
  end interface

  public :: create, get_algorithm_from_string

  contains

!---------------------------------------------------------------------------------------------------
! Externally called routine that ultimatly passes back a t_simulation pointer.
!---------------------------------------------------------------------------------------------------
  function create_simulation( opts )
!---------------------------------------------------------------------------------------------------
    use m_input_file
    use m_parameters
    implicit none

    type( t_options ), intent(in) :: opts
    class( t_simulation ), pointer :: create_simulation

    class( t_simulation ), pointer :: sim

    sim => create_core_simulation( opts )

    create_simulation => sim
  end function
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Utility function (called from os-main.f03) that maps simulation type strings
! to the internally used integer parameters
!---------------------------------------------------------------------------------------------------
function get_algorithm_from_string(simulation_label)
  use m_parameters
  implicit none
  character(len=20) :: simulation_label
  integer :: get_algorithm_from_string

  select case( trim( simulation_label ) )
    case( "standard" )
      get_algorithm_from_string = p_standard
    case( "pgc" )
      get_algorithm_from_string = p_sim_pgc
    case( "qed" )
      get_algorithm_from_string = p_sim_qed
    case( "shear" )
      get_algorithm_from_string = p_sim_shear
    case( "quasi-3D" )
      get_algorithm_from_string = p_sim_cyl_modes
    case( "radiative", "radiat" )
      get_algorithm_from_string = p_sim_rad
    case( "tiles" )
      get_algorithm_from_string = p_sim_tiles
    case( "split-damp", "overdense" )
      get_algorithm_from_string = p_sim_overdense
    case( "split-damp-q3d", "overdense-cyl" )
      get_algorithm_from_string = p_sim_overdense_cyl
    case( "neutral-spin" )
      get_algorithm_from_string = p_sim_neut_spin
    case( "exact-ext-emf", "xxfel" )
      get_algorithm_from_string = p_sim_xxfel
    case( "qed-q3d", "qed-cyl" )
      get_algorithm_from_string = p_sim_qed_cyl
    case( "radiative-q3d", "radiat-cyl" )
      get_algorithm_from_string = p_sim_rad_cyl
    case( "gr" )
      get_algorithm_from_string = p_sim_gr
    case( "cuda" )
      get_algorithm_from_string = p_sim_cuda
    case default
      call error_and_die( "Error reading global simulation parameters"//NEW_LINE('A')//&
                         &"algorithm: \'"//trim( simulation_label )//"\' is invalid or&
                         & the code is not compiled to support it. Aborting.")
  end select
end function

!---------------------------------------------------------------------------------------------------
! Core that creates the basic simulation types.
! This is the code that used to be in 'os-main%initialize'
!---------------------------------------------------------------------------------------------------
function create_core_simulation( opts )
!---------------------------------------------------------------------------------------------------
  use m_input_file
  use m_parameters
  implicit none

  type( t_options ), intent(in) :: opts
  class( t_simulation ), pointer :: create_core_simulation

  class( t_simulation ), pointer :: sim

    ! Allocate simulation object based on algorithm
  select case ( opts%algorithm )
    case( p_standard )
      ! Standard EM-PIC simulation
      allocate( t_simulation :: sim )
    case( p_sim_pgc )
      ! Ponderomotive guiding center simulation



        call error_and_die('(*error*) You requested a simulation using PGC but the code is not compiled with PGC support.')

    case( p_sim_qed )
      ! QED simulation

        allocate( t_simulation_qed :: sim )



    case( p_sim_shear )
      ! shear framework



        call error_and_die('(*error*) You requested a simulation using shear but the code is not compiled with shear support.')

    case( p_sim_cyl_modes )
      ! quasi-3D (cylindrical modes) simulation



        call error_and_die('(*error*) You requested a simulation using quasi-3D but the code is not compiled with quasi-3D support.')

    case( p_sim_rad )



        call error_and_die('(*error*) You requested a simulation using the radiative algorithm but the code is not&
          & compiled with radiative support.')

    case( p_sim_tiles )



        call error_and_die('(*error*) You requested a simulation using tiles but the code is not compiled with&
    & tiles support.')

    case( p_sim_overdense )



        call error_and_die('(*error*) You requested a simulation using the overdense algorithm but the code is&
    & not compiled with overdense support.')

    case( p_sim_overdense_cyl )



        call error_and_die('(*error*) You requested a simulation using the overdense-q3d algorithm but the code is&
    & not compiled with overdense-q3d support.')

    case( p_sim_neut_spin )

        allocate( t_simulation_neutral_spin:: sim )




    case( p_sim_xxfel )

        allocate( t_simulation_xxfel :: sim )




    case( p_sim_qed_cyl )



        call error_and_die('(*error*) You requested a simulation using the qed-q3d algorithm but the code is&
    & not compiled with qed-q3d support (must also have quasi-3D and QED support).')

    case( p_sim_rad_cyl )



        call error_and_die('(*error*) You requested a simulation using the radiative-q3d algorithm but the code is not&
          & compiled with radiative-q3d support.')

    case( p_sim_gr )

        allocate( t_simulation_gr :: sim )




    case( p_sim_cuda )



        call error_and_die('(*error*) You requested a simulation using the cuda &
                            & algorithm but the code is not compiled with cuda support')

    case default
      call error_and_die( '(*error*) Invalid simulation mode, only standard, pgc, qed, shear, quasi-3D,&
                          & radiative, tiles, overdense, overdense-q3d, neutral-spin, xxfel, qed-q3d, gr,&
                          & radiative-q3d and cuda are currently implemented.')
  end select
  create_core_simulation => sim
end function
!---------------------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------------------
! Utility function to print out fatal error conditions.
!---------------------------------------------------------------------------------------------------
subroutine error_and_die( message )
!---------------------------------------------------------------------------------------------------
    use m_parameters
    implicit none
    character(len=*) :: message

    if ( mpi_node() == 0 ) then
      write(0,*) message
      call abort_program( p_err_invalid )
    endif

  end subroutine
!---------------------------------------------------------------------------------------------------
end module


!---------------------------------------------------------------------------------------------------
! Version of the simulation factory that rests outside a module
! to avoid circular depencies in cases where a simulation wants to,
! inside it's own code, use the factory to create other simulation objects
! ( an example of this is the 'tandem' simulation class )
!---------------------------------------------------------------------------------------------------
! function global_simulation_factory( opts )

! use m_system, only: t_options
! use m_simulation, only: t_simulation
! use m_simulation_factory

! implicit none

! !integer, intent(in) :: algorithm
! type( t_options ), intent(inout) :: opts

! class( t_simulation ), pointer :: global_simulation_factory


! global_simulation_factory => create( opts )

! end function
