# 1 "qed-cyl/os-spec-define-qed-cyl.f03"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "qed-cyl/os-spec-define-qed-cyl.f03"
!-----------------------------------------------------------------------------------------

! QED species definition module

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
# 6 "qed-cyl/os-spec-define-qed-cyl.f03" 2
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
# 7 "qed-cyl/os-spec-define-qed-cyl.f03" 2

module m_species_define_qedcyl

# 1 "./memory/memory.h" 1
! Include file for the memory module
! Files must include this file using #include "memory/memory.h" rather than just using the module
! the #include must be placed where the module use statement would usually be i.e.
!
! module module2
!
! use module1
! #include "memory/memory.h"
!
# 19 "./memory/memory.h"
use memory
# 11 "qed-cyl/os-spec-define-qed-cyl.f03" 2

use m_system
use m_parameters
use m_species_define
use m_particles_cyl_modes, only: t_particles_cyl_modes
use m_species_cyl_modes_define, only: t_species_cyl_modes
use m_input_file, only: t_input_file, get_namelist
use m_vdf_report, only: new
use m_species_collisions, only: read_nml
use m_grid_define, only: t_grid
use m_space, only: t_space
use m_emf_define, only: t_emf
use m_node_conf, only: t_node_conf
use m_current_define, only: t_current
use m_bnd, only: t_bnd
use m_species, only: move_window
use m_vdf_comm, only : t_vdf_msg
use m_vdf_define
use m_time_step
use m_restart
use m_emf_cyl_modes, only : t_emf_cyl_modes
use m_bremsstrahlung_qedcyl, only : t_bremsstrahlung_info_qedcyl
use m_cyl_modes, only: t_cyl_modes, setup

implicit none

private

!-----------------------------------------------------------------------------------------
! t_species_qed
! - species with qed pusher and merging
!-----------------------------------------------------------------------------------------

type, extends( t_species_cyl_modes ) :: t_species_qedcyl

  real(p_k_part) :: n0 !reference density [cm^-3]
  real(p_k_part) :: omega_p0 !reference electron plasma frequency
  real(p_double) :: qed_g_cutoff, p_emit_cutoff ! cutoff parameters for qed emissions

  integer :: ndump_fac_rad

  ! emission flag
  logical :: if_damp_classical
  logical :: if_damp_qed

  ! radiated energy
  real(p_double), dimension(5) :: energy_rad = 0.0_p_double

  ! pointer to qed group photons
  class( t_species_cyl_modes ), pointer :: photons ! photons, always associated

  ! emission events diagnostic
  integer :: ndump_fac_n_emit
  type( t_vdf ) :: n_emit

  ! chi parameter diagnostic
  integer :: ndump_fac_chi
  type( t_vdf ) :: weight, weight_chi

  ! radiation spectrum diagnostic ndump
  integer :: ndump_fac_radspect, ndump_fac_chi_emit

  ! radiation and chi spectrum parameters
  real(p_double) :: radspect_emin, chi_emit_min
  real(p_double) :: radspect_emax, chi_emit_max
  integer :: radspect_bins, chi_emit_nbins

  ! radiation spectrum arrays
  real(p_double), dimension(:), pointer :: radspect => null()
  real(p_double), dimension(:), pointer :: radspect_oor => null()
  real(p_double), dimension(:), pointer :: chispect_photons => null()
  real(p_double), dimension(:), pointer :: chispect_leptons => null()

  ! spherical radiation spectrum diagnostic ndump
  integer :: ndump_fac_radanglespect

  ! spherical radiation spectrum parameters
  real(p_double) :: radanglespect_thresh
  integer :: radanglespect_bins

  ! spherical radiation spectrum arrays
  real(p_double), dimension(:,:,:), pointer :: radanglespect => null()

  ! spherical radiation spectrum diagnostic ndump
  integer :: ndump_fac_radsphere

  ! spherical radiation spectrum parameters
  real(p_double) :: radsphere_emin, radsphere_emax
  integer :: radsphere_nbins
  logical :: if_add_classical_radsphere

  ! spherical radiation spectrum arrays
  real(p_double), dimension(:,:,:), pointer :: radsphere => null()

  ! angular radiation detector ndump
  integer :: ndump_fac_raddetector

  ! radiation detector parameters
  real(p_double) :: raddetector_emin, raddetector_emax, raddetector_phimin, raddetector_phimax, &
                    raddetector_thetamin, raddetector_thetamax
  integer :: raddetector_nbins, raddetector_dir

  ! radiation detector array
  real(p_double), dimension(:,:), pointer :: raddetector => null()

  ! Bremsstrahlung photon creation
  type( t_bremsstrahlung_info_qedcyl ) :: bremsstrahlung_info

contains

  procedure :: allocate_objs => allocate_objs_spec_qedcyl
  procedure :: read_input => read_input_spec_qedcyl
  procedure :: set_dudt => set_dudt_spec_qedcyl
  procedure :: init => init_spec_qedcyl
  procedure :: list_algorithm => list_algorithm_spec_qedcyl
  procedure :: report => report_spec_qedcyl
  procedure :: report_energy_rad => report_energy_rad_spec_qedcyl
  procedure :: cleanup => cleanup_spec_qedcyl

end type t_species_qedcyl

interface
  subroutine report_spec_qedcyl( this, emf, g_space, grid, no_co, tstep, t, send_msg, recv_msg )

    import t_species_qedcyl, t_emf, t_space, t_grid, t_node_conf, t_time_step, p_double, t_vdf_msg

    class( t_species_qedcyl ), intent(inout) :: this

    class( t_emf ), intent(inout) :: emf
    type( t_space ), intent(in) :: g_space
    class( t_grid ), intent(in) :: grid
    class( t_node_conf ), intent(in) :: no_co
    type( t_time_step ), intent(in) :: tstep
    real(p_double), intent(in) :: t
    type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  end subroutine
end interface

interface
  subroutine set_dudt_spec_qedcyl( this )

    import t_species_qedcyl

    class( t_species_qedcyl ), intent(inout) :: this

  end subroutine
end interface

interface
  subroutine report_energy_rad_spec_qedcyl( this, no_co, tstep, t )

    import t_species_qedcyl, t_node_conf, t_time_step, p_double

    class( t_species_qedcyl ), intent(in) :: this

    class( t_node_conf ), intent(in) :: no_co
    type( t_time_step ), intent(in) :: tstep
    real(p_double), intent(in) :: t

  end subroutine
end interface

public t_species_qedcyl

! type diad QED
type, extends( t_diag_species ) :: t_diag_species_qedcyl

contains

  !procedure :: avail_report_quants => avail_report_quants_species_qed
  !procedure :: init_report_quants => init_report_quants_species_qed
  procedure :: read_input => read_input_diag_spec_qedcyl

end type t_diag_species_qedcyl

contains

!-----------------------------------------------------------------------------------------
! Allocates algorithm specific objects
!-----------------------------------------------------------------------------------------
subroutine allocate_objs_spec_qedcyl( this )

  implicit none

  class( t_species_qedcyl ), intent(inout) :: this

  ! local variables

  ! allocate algorithm specific t_diag_species object
  if ( .not. associated( this%diag ) ) then
      allocate( t_diag_species_qedcyl :: this%diag )
  endif

  ! setup radiated energy spectrum diagnostic
  if (this%ndump_fac_radspect > 0) then
      ! allocate arrays
      call alloc(this%radspect, (/ this%radspect_bins /),"qed-cyl/os-spec-define-qed-cyl.f03",208)
      call alloc(this%radspect_oor, (/ 2 /),"qed-cyl/os-spec-define-qed-cyl.f03",209)

      ! set all values to zero
      this%radspect = 0.0_p_double
      this%radspect_oor = 0
  endif

  ! setup chi spectrum diagnostic
  if (this%ndump_fac_chi_emit > 0) then
      ! allocate arrays
      call alloc(this%chispect_leptons, (/ this%chi_emit_nbins /),"qed-cyl/os-spec-define-qed-cyl.f03",219)
      call alloc(this%chispect_photons, (/ this%chi_emit_nbins /),"qed-cyl/os-spec-define-qed-cyl.f03",220)

      ! set all values to zero
      this%chispect_leptons = 0.0_p_double
      this%chispect_photons = 0.0_p_double

  endif

  ! setup a spectrum differential in energy and angle
  if (this%ndump_fac_radanglespect > 0) then
      ! allocate arrays
      call alloc(this%radanglespect, (/ this%radanglespect_bins, this%radanglespect_bins, 6 /),"qed-cyl/os-spec-define-qed-cyl.f03",231)

      ! set all values to zero
      this%radspect = 0.0_p_double
  endif

  ! setup spherical radiated energy spectrum diagnostic
  if (this%ndump_fac_radsphere > 0) then
      ! allocate arrays
      call alloc(this%radsphere, (/ this%radsphere_nbins, this%radsphere_nbins, 6 /),"qed-cyl/os-spec-define-qed-cyl.f03",240)

      ! set all values to zero
      this%radsphere = 0.0_p_double
  endif

    ! setup the detector diagnostic
  if (this%ndump_fac_raddetector > 0) then
      ! allocate arrays
      call alloc(this%raddetector, (/ this%raddetector_nbins, this%raddetector_nbins /),"qed-cyl/os-spec-define-qed-cyl.f03",249)

      ! set all values to zero
      this%raddetector = 0.0_p_double
  endif

  ! no need to allocate superclass objects, this is done on t_species % read_input()

end subroutine allocate_objs_spec_qedcyl
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Read qed species information from input file
!-----------------------------------------------------------------------------------------
subroutine read_input_spec_qedcyl( this, input_file, def_name, periodic, if_move, grid, &
  dt, read_prof, sim_options )

    class( t_species_qedcyl ), intent(inout) :: this
    class( t_input_file ), intent(inout) :: input_file

    character(len = *), intent(in) :: def_name
    logical, dimension(:), intent(in) :: periodic, if_move
    class( t_grid ), intent(in) :: grid

    real(p_double), intent(in) :: dt
    logical, intent(in) :: read_prof

    type(t_options), intent(in) :: sim_options

    ! call superclass read_input
    call this % t_species_cyl_modes % read_input( input_file, def_name, periodic, if_move, grid, &
                            dt, read_prof, sim_options )

    ! allocate algorithm specific objects
    call this % allocate_objs()

    ! set the omega_p0 and n0 variables
    this % omega_p0 = real( sim_options % omega_p0, p_k_part )
    this % n0 = real( sim_options % n0, p_k_part )

end subroutine read_input_spec_qedcyl
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Initialize QED species specific diagnostics
!-----------------------------------------------------------------------------------------
subroutine read_input_diag_spec_qedcyl( this, input_file, gamma )

  use m_species_diagnostics

  class ( t_diag_species_qedcyl ), intent(inout) :: this
  class( t_input_file ), intent(inout) :: input_file
  real(p_double), intent(in) :: gamma

  call this % t_diag_species % read_input ( input_file, gamma )

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Printout the algorithm used by the pusher
!-----------------------------------------------------------------------------------------
subroutine list_algorithm_spec_qedcyl( this )

  implicit none
  class( t_species_qedcyl ), intent(in) :: this

  if (mpi_node()==0) print *,"Species : ", trim(this%name)

  ! Push type
  select case ( this%push_type )
  case (p_std)
      print *, '   - QED standard (Boris) pusher'
  case(p_exact)
      print *, '   - QED exact pusher'
  case(p_exact_rr)
      print *, '   - QED exact pusher with radiation reaction'
  end select

  if ( this%free_stream ) then
      print *, '- Free streaming particles (no dudt)'
      write(0,*) '- (*warning*) ', trim(this%name), ' are free streaming!'
  endif

end subroutine list_algorithm_spec_qedcyl
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Cleans up algorithm specific objects
!-----------------------------------------------------------------------------------------
subroutine cleanup_spec_qedcyl( this )

  implicit none

  class(t_species_qedcyl) , intent(inout) :: this

  ! cleanup local data structures
  ! call this % emit % cleanup()

  ! cleanup superclass
  call this % t_species_cyl_modes % cleanup()

  ! cleanup radiated energy spectrum diagnostic arrays
  if (this%ndump_fac_radspect > 0) then
    ! allocate arrays
    call freemem(this%radspect,"qed-cyl/os-spec-define-qed-cyl.f03",354)
    call freemem(this%radspect_oor,"qed-cyl/os-spec-define-qed-cyl.f03",355)
  endif

  ! cleanup chi spectrum diagnostic arrays
  if (this%ndump_fac_chi_emit > 0) then
    ! allocate arrays
    call freemem(this%chispect_leptons,"qed-cyl/os-spec-define-qed-cyl.f03",361)
    call freemem(this%chispect_photons,"qed-cyl/os-spec-define-qed-cyl.f03",362)
  endif

  ! cleanup spherical radiated energy spectrum diagnostic arrays
  if (this%ndump_fac_radanglespect > 0) then
    ! allocate arrays
    call freemem(this%radanglespect,"qed-cyl/os-spec-define-qed-cyl.f03",368)
  endif

  ! cleanup spherical radiated energy spectrum diagnostic arrays
  if (this%ndump_fac_radsphere > 0) then
    ! deallocate arrays
    call freemem(this%radsphere,"qed-cyl/os-spec-define-qed-cyl.f03",374)
  endif

  ! cleanup detector diagnostic array
  if (this%ndump_fac_raddetector > 0) then
    ! deallocate arrays
    call freemem(this%raddetector,"qed-cyl/os-spec-define-qed-cyl.f03",380)
  endif

  ! cleanup emission events diagnostic vdf
  if (this%ndump_fac_n_emit > 0) then
    call this % n_emit % cleanup()
  endif

  ! cleanup chi parameter diagnostic vdf
  if (this%ndump_fac_chi > 0) then
    call this % weight_chi % cleanup()
    call this % weight % cleanup()
  endif

end subroutine cleanup_spec_qedcyl
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Initializes qed-cyl species
!-----------------------------------------------------------------------------------------
subroutine init_spec_qedcyl( this, sp_id, interpolation, grid_center, grid, g_space, emf, jay, &
  no_co, send_vdf, recv_vdf, bnd_cross, node_cross, send_spec, &
  recv_spec, ndump_fac, restart, restart_handle, sim_options, tstep, tmax )
!-----------------------------------------------------------------------------------------
! sets up this data structure from the given information
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_species_qedcyl ), intent(inout) :: this

  integer, intent(in) :: sp_id
  integer, intent(in) :: interpolation
  logical, intent(in) :: grid_center
  class( t_grid ), intent(in) :: grid
  type( t_space ), intent(in) :: g_space
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
  integer :: n_modes

  n_modes = this % n_cyl_modes

  ! call superclass init
  call this % t_species_cyl_modes % init( sp_id, interpolation, grid_center, grid, g_space, emf, jay, &
      no_co, send_vdf, recv_vdf, bnd_cross, node_cross, send_spec, &
      recv_spec, ndump_fac, restart, restart_handle, sim_options, tstep, tmax )

  ! overwrite pusher option
  call this % set_dudt()

  ! reset radiated energy diagnostic
  this % energy_rad = 0.0_p_double

  ! setup emission events diagnostic
  if (this%ndump_fac_n_emit > 0) then
      gc_num = 0
      dx = 0.
      call this % n_emit % new( p_x_dim, 1, this%my_nx_p(3,1:p_x_dim), gc_num, dx, .true. )
  endif

  ! setup chi parameter diagnostic
  if (this%ndump_fac_chi > 0) then
      gc_num = 0
      dx = 0.
      call this % weight % new( p_x_dim, 1, this%my_nx_p(3,1:p_x_dim), gc_num, dx, .true. )
      call this % weight_chi % new( p_x_dim, 1, this%my_nx_p(3,1:p_x_dim), gc_num, dx, .true. )
  endif

  ! setup vdf with ion density for bremsstrahlung
  if (this%bremsstrahlung_info%if_bremsstrahlung) then
      ! The number of guard cells required depends on the interpolation level
      gc_num(p_lower,:) = this%interpolation
      gc_num(p_upper,:) = this%interpolation + 1
      dx = 0.
      call this % bremsstrahlung_info % rho_ion % new( p_x_dim, 1, this%my_nx_p(3,1:p_x_dim), gc_num, dx, .true. )
      call setup( this % bremsstrahlung_info % rho_ion_cyl_m, n_modes, this % bremsstrahlung_info % rho_ion )
  endif

end subroutine init_spec_qedcyl
!-----------------------------------------------------------------------------------------

end module m_species_define_qedcyl
