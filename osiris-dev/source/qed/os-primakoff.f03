#include "os-preprocess.fpp"
#include "os-config.h"

module m_primakoff

#include "memory/memory.h"
  !use m_random_class

  use m_parameters
  use m_system
  use m_node_conf
  use m_random
  use m_vdf_define
  use m_vdf
  use m_species_define
  use m_species_define_qed 	! Not sure this is needed

implicit none

private

! -----
! CLASS DEFINITIONS
! -----

type :: t_primakoff_info
  logical :: if_primakoff
  integer :: Z_ion		! Not sure this is needed
  type( t_vdf ) :: rho_ion	! Not sure this is needed 


  ! multiplicative factor applied on the probability
  real( p_double ) :: proba_mult

  ! Type of cross section
  integer :: cs_model
  !logical :: energy_damp 	! I don't think we need this for primakoff. We might need it for axion bremsstrahlung (See ReadME)
  logical :: decay_photons
  real( p_double ) :: p_emit_cutoff	! I don't think we need this for primakoff. We might need it for axion bremsstrahlung (See ReadME)


  ! pointer to the ion species for primakoff
  class(t_species), pointer :: ion

  integer, dimension(:), allocatable :: i_ion


  contains


  procedure :: read_input => read_input_primakoff

end type

public :: t_primakoff_info

contains


!-----------------------------------------------------------------------------------------
  ! Reads info for primakoff
  !-----------------------------------------------------------------------------------------
  subroutine read_input_primakoff( this, input_file, num_qed)
  use m_input_file, only : t_input_file, get_namelist

  implicit none

  class( t_primakoff_info ), intent(inout) :: this
  type( t_input_file ), intent(inout) :: input_file
  integer, intent(in) :: num_qed

  logical :: if_primakoff, decay_photons
  integer :: i, ierr, Z_ion !, nb_couples !I don't think we need this anywhere. 
! I guess we don't neither need "i" defined in the line just above this. 
  integer :: cs_model
  real(p_double) :: proba_mult
  integer, dimension(:), allocatable :: i_ion


  namelist /nl_qed_primakoff/ if_primakoff, Z_ion, proba_mult, cs_model, decay_photons, i_ion

  if_primakoff = .false.
  Z_ion = -1
  proba_mult = 1.0
  cs_model = -1
  decay_photons = .true.
  allocate(i_ion(num_qed))
  i_ion = 0


  ! Get namelist text from input file
  call get_namelist( input_file, "nl_qed_primakoff", ierr )
  
  if ( ierr == 0 ) then
      read (input_file%nml_text, nml = nl_qed_primakoff, iostat = ierr)
      if (ierr /= 0) then
          print *, ""
          print *, "   Error reading qed_primakoff parameters "
          print *, "   aborting..."
          stop
      endif
  else
      SCR_ROOT("   - no particle primakoff specified")
  endif


  this%if_primakoff = if_primakoff
  this%Z_ion = Z_ion
  this%proba_mult = proba_mult
  this%cs_model = cs_model 	!   this%cs_model = int(cs_model)
  this%decay_photons = decay_photons
  this%i_ion = i_ion


  end subroutine read_input_primakoff

end module m_primakoff
