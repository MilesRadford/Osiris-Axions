#include "os-preprocess.fpp"
#include "os-config.h"

module m_axions_define

#include "memory/memory.h"

!-----------------------------------------------------------------------------------------

use m_system
use m_parameters
use m_species_define, only: t_species, t_part_idx, t_spec_msg
use m_species_define_qed, only: t_species_qed
use m_psource_std
use m_species_udist
use m_mergedel, only: t_merge_info

use m_space
use m_grid_define
use m_node_conf
use m_restart
use m_time_step

use m_emf_define
use m_vdf_define
use m_vdf_comm, only : t_vdf_msg
use m_input_file

!-----------------------------------------------------------------------------------------
! need to be defined

use m_a_brems
use m_a_compt
use m_a_prima

!----------------------------------------------------------------------------------------

implicit none

private

! string to id restart data
character(len=*), parameter :: p_axions_rst_id = "axions rst data - 0x0000"

!-----------------------------------------------------------------------------------------
! t_axions
!  - axion particle species
!-----------------------------------------------------------------------------------------

type, extends( t_species_qed ) :: t_axions

! replace pp with decay to photons
! procedures

end type t_axions

!-----------------------------------------------------------------------------------------
! Read axion species information from input file
!-----------------------------------------------------------------------------------------

subroutine read_input_axions( this, input_file, def_name, periodic, if_move, grid, &
  dt, read_prof, sim_options )

    use m_system
    use m_parameters

    use m_qed
    use m_input_file
    use m_species_define
    use m_grid_define

    use m_profile, only : new_source

    implicit none

    class( t_axions ),  intent(inout) :: this
    class( t_input_file ), intent(inout) :: input_file

    character(len = *),    intent(in) :: def_name
    logical, dimension(:), intent(in) :: periodic, if_move
    class( t_grid ),               intent(in) :: grid

    real(p_double), intent(in) :: dt
    logical, intent(in) :: read_prof

    type(t_options), intent(in) :: sim_options
    
end subroutine read_input_axions

!-----------------------------------------------------------------------------------------
! Printout the algorithm used by the pusher
!-----------------------------------------------------------------------------------------

subroutine list_algorithm_axions( this )

  implicit none
  class( t_axions ), intent(in) :: this

  SCR_ROOT("Species : ", trim(this%name))
  print *, '   - axion QED pusher'
  

end subroutine list_algorithm_axions

!-----------------------------------------------------------------------------------------
! Initialize object data
!-----------------------------------------------------------------------------------------

subroutine init_axions( this, sp_id, interpolation, grid_center, grid, g_space, emf, jay, &
  no_co, send_vdf, recv_vdf, bnd_cross, node_cross, send_spec, &
  recv_spec, ndump_fac, restart, restart_handle, sim_options, tstep, tmax )

    use m_grid_define, only : t_grid
    use m_space, only : t_space
    use m_emf_define, only : t_emf
    use m_current_define, only : t_current
    use m_node_conf

    use m_restart

    use m_species_define
    use m_system, only : t_options
    use m_vdf_comm, only : t_vdf_msg

    implicit none

    ! dummy variables
    class( t_axions ), intent(inout) :: this

    integer,     intent(in) :: sp_id
    integer, intent(in) :: interpolation
    logical, intent(in) :: grid_center
    class( t_grid ), intent(in)     :: grid
    type( t_space ),     intent(in) :: g_space
    class( t_emf ), intent(inout) :: emf
    class( t_current ), intent(inout) :: jay
    class( t_node_conf ), intent(in) :: no_co
    type(t_vdf_msg), dimension(2), intent(inout) :: send_vdf, recv_vdf
    type( t_part_idx ), dimension(2), intent(inout) :: bnd_cross
    type( t_part_idx ), intent(inout) :: node_cross
    type( t_spec_msg ), dimension(2), intent(inout) :: send_spec, recv_spec
    integer, intent(in) :: ndump_fac
    logical, intent(in) :: restart
    type( t_restart_handle ), intent(in) :: restart_handle
    type(t_options), intent(in) :: sim_options
    type( t_time_step ), intent(in) :: tstep
    real(p_double), intent(in) :: tmax 

    ! local variables
    integer, dimension(2, p_x_dim) :: gc_num
    real( p_double ), dimension(p_x_dim) :: dx

    ! call superclass initialization
    call this % t_species % init( sp_id, interpolation, grid_center, grid, g_space, emf, jay, &
        no_co, send_vdf, recv_vdf, bnd_cross, node_cross, send_spec, &
        recv_spec, ndump_fac, restart, restart_handle, sim_options, tstep, tmax )
        
! set up emission events-------------



!------------------------------------   

    ! setup local data
    if ( restart ) then
        call this % read_checkpoint( restart_handle )
    endif
        
end subroutine init_axions

!-----------------------------------------------------------------------------------------
! Return axion energy on local node, calculating it if need be
!-----------------------------------------------------------------------------------------

! set up once physics included

!-----------------------------------------------------------------------------------------
! Return pair energy on local node  (replace pairprod with decay to photons?)
!-----------------------------------------------------------------------------------------

! most likely photon pair energies









































