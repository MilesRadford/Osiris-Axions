# 1 "cyl_modes/os-boosted-diag-cyl-modes.f03"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "cyl_modes/os-boosted-diag-cyl-modes.f03"
! this module defines data structures that wrap the t_boosted_diag class
! so that it can be used for dumping boosted quasi-3d data.

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
# 5 "cyl_modes/os-boosted-diag-cyl-modes.f03" 2
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
# 6 "cyl_modes/os-boosted-diag-cyl-modes.f03" 2


module m_boosted_diag_cyl_modes

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
# 11 "cyl_modes/os-boosted-diag-cyl-modes.f03" 2

use m_system
use m_parameters

use m_boosted_grid_mgr
use m_boosted_diag
use m_cyl_modes


implicit none


!-------------------------------------------------------------------------------
type, extends( t_boosted_diag ) :: t_boosted_diag_cyl_modes_helper

    integer :: mode_number
    logical :: mode_is_real

contains

    procedure :: get_dir_fname_label => get_dir_fname_label_boosted_diag_cyl_modes_helper

end type t_boosted_diag_cyl_modes_helper
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
type, extends( t_boosted_diag ) :: t_boosted_diag_cyl_modes

    integer :: n_cyl_modes
    class( t_boosted_diag_cyl_modes_helper ), dimension(:), pointer :: pf_re => null()
    class( t_boosted_diag_cyl_modes_helper ), dimension(:), pointer :: pf_im => null()

    contains

    procedure :: restart_write => restart_write_boosted_diag_cyl_modes

    procedure :: allocate_objs => allocate_objs_boosted_diag_cyl_modes
    procedure :: init => init_boosted_diag_cyl_modes
    procedure :: init_helper => init_helper_boosted_diag_cyl_modes
    procedure :: cleanup => cleanup_boosted_diag_cyl_modes

    procedure :: create_files => create_files_boosted_diag_cyl_modes
    procedure :: report => report_boosted_diag_cyl_modes

    ! these procedures aren't overwritten because they have different args,
    ! specifically a t_cyl_modes instead of t_vdf as input
    procedure :: extract_cyl_modes => extract_boosted_diag_cyl_modes

end type t_boosted_diag_cyl_modes
!-------------------------------------------------------------------------------


contains


!-------------------------------------------------------------------------------
subroutine get_dir_fname_label_boosted_diag_cyl_modes_helper( this, component, t2_idx, dir, fname, label )

    implicit none

    class(t_boosted_diag_cyl_modes_helper), intent(in) :: this
    integer, intent(in) :: t2_idx
    integer, intent(in) :: component
    character(len=1024 ), intent(out) :: dir
    character(len=128 ), intent(out) :: fname
    character(len=128 ), intent(out) :: label

    character(len=128) :: name
    character(len=4) :: cyl_m_type_dir, cyl_m_type_file

    name = this%component_names( component )

    if( this%mode_is_real ) then
        cyl_m_type_dir = 'RE'
        cyl_m_type_file = 're'
    else
        cyl_m_type_dir = 'IM'
        cyl_m_type_file = 'im'
    endif

    write( dir, '(A, A, A,I0,A,A,A,A,A)' ) trim(this%basepath), p_dir_sep, 'MODE-', this%mode_number, &
            '-', trim(cyl_m_type_dir), p_dir_sep, trim(name), p_dir_sep

    write( fname, '(A,A,I0,A,A,A,A)' ) trim(name), '-', this%mode_number, '-', trim(cyl_m_type_file), '-', idx_string( t2_idx, 6 )

    label = trim(name)

end subroutine get_dir_fname_label_boosted_diag_cyl_modes_helper
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------------------------
subroutine restart_write_boosted_diag_cyl_modes( this, restart_handle )

    use m_restart

    implicit none

    class( t_boosted_diag_cyl_modes ), intent(in) :: this
    type( t_restart_handle ), intent(inout) :: restart_handle

    integer :: i

    if (mpi_node()==0) print *,'calling restart_write_boosted_diag'

    call this % pf_re(0) % restart_write( restart_handle )

    do i = 1, this % n_cyl_modes
        call this % pf_re(i) % restart_write( restart_handle )
        call this % pf_im(i) % restart_write( restart_handle )
    enddo

end subroutine restart_write_boosted_diag_cyl_modes
!-------------------------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
subroutine allocate_objs_boosted_diag_cyl_modes( this )

    implicit none

    class( t_boosted_diag_cyl_modes ), intent(inout) :: this

    call this % t_boosted_diag % allocate_objs()

    allocate( this%pf_re( 0 : this%n_cyl_modes ) )
    if (this%n_cyl_modes>0) allocate( this%pf_im( 1 : this%n_cyl_modes ) )

    ! the superclass allocate_objs is called in init for each of these boosted diags.

end subroutine allocate_objs_boosted_diag_cyl_modes
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
! assumes alloc_objs has already been called
!-------------------------------------------------------------------------------
subroutine init_boosted_diag_cyl_modes( this, grid, g_space, dt1, tmax1, &
                                        ndump_fac_global, basepath, ref_vdf, &
                                        gc_num, restart, restart_handle, init_grid_mgr_ )

    implicit none

    class( t_boosted_diag_cyl_modes ), intent(inout) :: this
    class( t_grid ), intent(in) :: grid
    type( t_space ), intent(in) :: g_space
    real(p_double), intent(in) :: dt1, tmax1
    integer, intent(in) :: ndump_fac_global
    class(t_vdf), intent(in) :: ref_vdf
    character(len=*), intent(in) :: basepath
    integer, dimension (2, p_x_dim), intent(in) :: gc_num
    logical, intent(in) :: restart
    type( t_restart_handle ), intent(in) :: restart_handle
    logical, intent(in), optional :: init_grid_mgr_

    integer :: i, ndump_total

    ndump_total = this%ndump_fac * ndump_fac_global

    call this % allocate_objs()

    ! init our own grid manager, but don't allocate our own storage buffs
    ! this one grid manager will be shared by all the member boosted_diags
    call this % init_grid_mgr( grid, g_space, dt1, tmax1 )
    call this % compute_data_buf_size( ndump_total )

    call this % init_components_to_buffer()

    call this % init_helper( this%pf_re(0), 0, .true. )

    call this % pf_re(0) % init( grid, g_space, dt1, tmax1, &
                 ndump_fac_global, &
                 basepath, ref_vdf, gc_num, restart, restart_handle, .false. )

    do i = 1, this % n_cyl_modes

        call this % init_helper( this%pf_re(i), i, .true. )

        call this % pf_re(i) % init( grid, g_space, dt1, tmax1, &
                 ndump_fac_global, &
                 basepath, ref_vdf, gc_num, restart, restart_handle, .false. )

        call this % init_helper( this%pf_im(i), i, .false. )

        call this % pf_im(i) % init( grid, g_space, dt1, tmax1, &
                 ndump_fac_global, &
                 basepath, ref_vdf, gc_num, restart, restart_handle, .false. )

    enddo

end subroutine init_boosted_diag_cyl_modes
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
subroutine cleanup_boosted_diag_cyl_modes( this )

    implicit none

    class( t_boosted_diag_cyl_modes ), intent(inout) :: this

    integer :: i

    ! cleanup our grid manager first
    ! after this it won't be associated for any of the pf_re or pf_im
    if( associated(this%grid_mgr) ) then
        call this%grid_mgr%cleanup()
        deallocate( this%grid_mgr )
    endif

    ! No storage buffers to cleanup here

    if (associated(this%pf_re)) then
        this % pf_re(0) % grid_mgr => null()
        call this % pf_re(0) % cleanup()

        do i = 1, this % n_cyl_modes
            this % pf_re(i) % grid_mgr => null()
            this % pf_im(i) % grid_mgr => null()
            call this % pf_re(i) % cleanup()
            call this % pf_im(i) % cleanup()
        enddo

        deallocate( this%pf_re )
        if (this%n_cyl_modes>0) deallocate( this%pf_im )
    endif

end subroutine cleanup_boosted_diag_cyl_modes
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
subroutine init_helper_boosted_diag_cyl_modes( this, boosted_diag, mode_number, mode_is_real )

    implicit none

    class( t_boosted_diag_cyl_modes ), intent(inout) :: this
    class( t_boosted_diag_cyl_modes_helper ), intent(inout) :: boosted_diag
    integer, intent(in) :: mode_number
    logical, intent(in) :: mode_is_real

    ! these are different for the different angular diags
    boosted_diag % mode_number = mode_number
    boosted_diag % mode_is_real = mode_is_real

    ! these are the same
    boosted_diag % diag_type = this % diag_type
    boosted_diag % num_diag_components = this % num_diag_components
    boosted_diag % diag_components = this % diag_components
    boosted_diag % num_diag_components_to_buffer = this % num_diag_components_to_buffer
    boosted_diag % diag_components_to_buffer = this % diag_components_to_buffer
    boosted_diag % component_names = this % component_names
    boosted_diag % ndump_fac = this % ndump_fac
    boosted_diag % data_buf_size = this % data_buf_size
    boosted_diag % grid_mgr => this % grid_mgr

end subroutine init_helper_boosted_diag_cyl_modes
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
subroutine create_files_boosted_diag_cyl_modes( this, grid )

    implicit none

    class( t_boosted_diag_cyl_modes ), intent(in) :: this
    class( t_grid ), intent(in) :: grid

    integer :: i

    call this % pf_re(0) % create_files( grid )

    do i = 1, this % n_cyl_modes
        call this % pf_re(i) % create_files( grid )
        call this % pf_im(i) % create_files( grid )
    enddo

end subroutine create_files_boosted_diag_cyl_modes
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------

subroutine report_boosted_diag_cyl_modes( this, tstep, grid )

    implicit none

    class( t_boosted_diag_cyl_modes ), intent(inout) :: this
    type( t_time_step ), intent(in) :: tstep
    class( t_grid ), intent(in) :: grid

    integer :: i

    if( n( tstep ) .eq. 0 ) then

        if (mpi_node()==0) print *,'Creating boosted files...'
        call this % create_files( grid )
        if (mpi_node()==0) print *,'Done.'

    else

        call this % pf_re(0) % report( tstep, grid )

        do i = 1, this % n_cyl_modes
            call this % pf_re(i) % report( tstep, grid )
            call this % pf_im(i) % report( tstep, grid )
        enddo

    endif

end subroutine report_boosted_diag_cyl_modes
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
subroutine extract_boosted_diag_cyl_modes( this, tstep, t1, grid, g_space, cyl_modes1, cyl_modes2 )

    implicit none

    class( t_boosted_diag_cyl_modes ), intent(inout) :: this
    type(t_time_step), intent(in) :: tstep
    real(p_double), intent(in) :: t1
    class( t_grid ), intent(in) :: grid
    type( t_space ), intent(in) :: g_space
    class(t_cyl_modes), intent(in) :: cyl_modes1
    class(t_cyl_modes), intent(in), optional :: cyl_modes2

    integer :: i

    if( present( cyl_modes2 ) ) then

        call this % pf_re(0) % extract( tstep, t1, grid, g_space, cyl_modes1%pf_re(0), cyl_modes2%pf_re(0) )

        do i = 1, this % n_cyl_modes
            call this % pf_re(i) % extract( tstep, t1, grid, g_space, cyl_modes1%pf_re(i), cyl_modes2%pf_re(i) )
            call this % pf_im(i) % extract( tstep, t1, grid, g_space, cyl_modes1%pf_re(i), cyl_modes2%pf_re(i) )
        enddo

    else

        call this % pf_re(0) % extract( tstep, t1, grid, g_space, cyl_modes1%pf_re(0) )

        do i = 1, this % n_cyl_modes
            call this % pf_re(i) % extract( tstep, t1, grid, g_space, cyl_modes1%pf_re(i) )
            call this % pf_im(i) % extract( tstep, t1, grid, g_space, cyl_modes1%pf_re(i) )
        enddo

    endif

end subroutine extract_boosted_diag_cyl_modes
!-------------------------------------------------------------------------------



end module m_boosted_diag_cyl_modes
