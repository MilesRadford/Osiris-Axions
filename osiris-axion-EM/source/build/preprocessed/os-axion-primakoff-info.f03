# 1 "qed/os-axion-primakoff-info.f03"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "qed/os-axion-primakoff-info.f03"
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
# 2 "qed/os-axion-primakoff-info.f03" 2
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
# 3 "qed/os-axion-primakoff-info.f03" 2

module m_axion_primakoff_info

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
# 7 "qed/os-axion-primakoff-info.f03" 2

  use m_parameters
  use m_vdf_define, only: t_vdf
  use m_input_file, only: t_input_file, get_namelist
  use m_species_define, only: t_species
  use m_system, only: mpi_node

  implicit none
  private

  type, public :: t_axion_primakoff_info

    logical :: if_primakoff = .false.
    logical :: if_inverse = .false.

    ! coupling g_{aγ} [GeV^-1] (user supplies; operator handles units consistently with rest of module)
    real( p_double ) :: g_agamma = 0.0_p_double

    ! plasma controls
    real( p_double ) :: Te_keV = 0.0_p_double
    logical :: include_ee = .true.
    logical :: use_screening = .true.

    logical :: assume_neutral = .true.
    real( p_double ) :: ne_factor = 1.0_p_double

    ! variance reduction
    real( p_double ) :: w_axion = 1.0_p_double
    real( p_double ) :: lambda_cap = 20.0_p_double

    ! optional fixed screening scale
    real( p_double ) :: kappa_keV = -1.0_p_double

    ! ion charge (σ ∝ Z^2)
    integer :: Z_ion = 1

    ! optional debug prints
    logical :: debug = .false.

! optional photon initialization override (testing / validation)
! photon_init: 'none' (default) or 'planck' to resample photon momenta once at start
character(len=16) :: photon_init = 'none'
real( p_double ) :: photon_Te_keV = -1.0_p_double
logical :: photon_quiet_pairs= .true.
integer :: photon_nmax = 200
logical :: photon_init_done = .false.


    ! ion coupling (per QED group)
    integer, allocatable :: i_ion(:)

    ! runtime wiring
    class( t_species ), pointer :: ion => null()

    ! ion density VDF (allocated/deposited elsewhere; mirrored BH pattern)
    type( t_vdf ) :: rho_ion

  contains
    procedure, public :: set_defaults => primakoff_info_set_defaults
    procedure, public :: sanity_check => primakoff_info_sanity_check
    procedure, public :: read_input => primakoff_info_read_input
  end type t_axion_primakoff_info

contains

  subroutine primakoff_info_set_defaults(this)
    class( t_axion_primakoff_info ), intent(inout) :: this
    this%if_primakoff = .false.
    this%if_inverse = .false.
    this%g_agamma = 0.0_p_double
    this%Te_keV = 0.0_p_double
    this%include_ee = .true.
    this%use_screening = .true.
    this%assume_neutral = .true.
    this%ne_factor = 1.0_p_double
    this%w_axion = 0.0_p_double
    this%lambda_cap = 0.0_p_double
    this%kappa_keV = -1.0_p_double
    this%Z_ion = 1
    this%debug = .false.
    this%photon_init = 'none'
    this%photon_Te_keV = -1.0_p_double
    this%photon_quiet_pairs= .true.
    this%photon_nmax = 200
    this%photon_init_done = .false.
    if (allocated(this%i_ion)) this%i_ion = 1
    nullify(this%ion)
  end subroutine primakoff_info_set_defaults

  subroutine primakoff_info_sanity_check(this)
    class( t_axion_primakoff_info ), intent(inout) :: this
    if (this%w_axion <= 0.0_p_double) this%w_axion = 1.0_p_double
    if (this%lambda_cap <= 0.0_p_double) this%lambda_cap = 20.0_p_double
    if (this%ne_factor <= 0.0_p_double) this%ne_factor = 1.0_p_double
    if (this%Z_ion < 1) this%Z_ion = 1
    if (this%photon_nmax < 1) this%photon_nmax = 1
  end subroutine primakoff_info_sanity_check

  subroutine primakoff_info_read_input(this, input_file, num_qed)
    class( t_axion_primakoff_info ), intent(inout) :: this
    type( t_input_file ), intent(inout) :: input_file
    integer, intent(in) :: num_qed

    integer :: ierr

    logical :: if_primakoff, if_inverse, include_ee, use_screening, assume_neutral, debug
    logical :: photon_quiet_pairs
    real( p_double ) :: g_agamma, Te_keV, ne_factor, w_axion, lambda_cap, kappa_keV
    real( p_double ) :: photon_Te_keV
    integer :: Z_ion
    integer :: photon_nmax
    character(len=16) :: photon_init
    integer, allocatable :: i_ion(:)

    ! IMPORTANT:
    ! - follow the working axion-brems convention: namelist group name starts with "nl_"
    namelist /nl_qed_axion_primakoff/ &
      if_primakoff, if_inverse, g_agamma, Te_keV, include_ee, use_screening, &
      assume_neutral, ne_factor, w_axion, lambda_cap, kappa_keV, Z_ion, debug, i_ion, &
      photon_init, photon_Te_keV, photon_quiet_pairs, photon_nmax

    call this%set_defaults()

    if (num_qed < 1) then
      call this%sanity_check()
      return
    end if

    allocate(i_ion(num_qed))
    i_ion = 1

    ! seed locals from defaults
    if_primakoff = this%if_primakoff
    if_inverse = this%if_inverse
    g_agamma = this%g_agamma
    Te_keV = this%Te_keV
    include_ee = this%include_ee
    use_screening = this%use_screening
    assume_neutral = this%assume_neutral
    ne_factor = this%ne_factor
    w_axion = this%w_axion
    lambda_cap = this%lambda_cap
    kappa_keV = this%kappa_keV
    Z_ion = this%Z_ion
    debug = this%debug


    photon_init = this%photon_init
    photon_Te_keV = this%photon_Te_keV
    photon_quiet_pairs = this%photon_quiet_pairs
    photon_nmax = this%photon_nmax
    call get_namelist( input_file, "nl_qed_axion_primakoff", ierr )

    if ( ierr == 0 ) then
      read (input_file%nml_text, nml = nl_qed_axion_primakoff, iostat = ierr)
      if ( ierr /= 0 ) then
        print *, ""
        print *, "   Error reading qed_axion_primakoff parameters "
        print *, "   aborting..."
        stop
      endif
    else
      if (mpi_node()==0) print *,"   - no particle axion primakoff specified"
    endif

    ! commit
    this%if_primakoff = if_primakoff
    this%if_inverse = if_inverse
    this%g_agamma = g_agamma
    this%Te_keV = Te_keV
    this%include_ee = include_ee
    this%use_screening = use_screening
    this%assume_neutral = assume_neutral
    this%ne_factor = ne_factor
    this%w_axion = w_axion
    this%lambda_cap = lambda_cap
    this%kappa_keV = kappa_keV
    this%Z_ion = Z_ion
    this%debug = debug


    this%photon_init = photon_init
    this%photon_Te_keV = photon_Te_keV
    this%photon_quiet_pairs = photon_quiet_pairs
    this%photon_nmax = photon_nmax
    this%photon_init_done = .false.
    if (.not. allocated(this%i_ion)) then
      allocate(this%i_ion(num_qed))
    else if (size(this%i_ion) /= num_qed) then
      deallocate(this%i_ion)
      allocate(this%i_ion(num_qed))
    end if
    this%i_ion = i_ion

    call this%sanity_check()

    if ( mpi_node() == 0 ) then
      write(*,*) "qed_axion_primakoff read: if_prim=", this%if_primakoff, " if_inv=", this%if_inverse, &
                 " g=", this%g_agamma, " w=", this%w_axion, " photon_init=", trim(this%photon_init)
    end if

    deallocate(i_ion)
  end subroutine primakoff_info_read_input

end module m_axion_primakoff_info
