# 0 "spec/diagnostics/os-spec-rawdiag.f03"
# 0 "<built-in>"
# 0 "<command-line>"
# 1 "spec/diagnostics/os-spec-rawdiag.f03"
# 1 "./os-config.h" 1
! Configuration file for osiris

! ----------------------------------------------------------------------------------------
! Algorithm options
! ----------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------------
! System options
! ----------------------------------------------------------------------------------------

! MPI supports MPI_IN_PLACE in global operations
!#define 1

! Use OpenMP
!#define 1

! Use SIMD optimized code
!#define SIMD
!#define SIMD_SSE
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
# 2 "spec/diagnostics/os-spec-rawdiag.f03" 2
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
# 3 "spec/diagnostics/os-spec-rawdiag.f03" 2

module m_species_rawdiag

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
# 7 "spec/diagnostics/os-spec-rawdiag.f03" 2

use m_species_define
use m_fparser

use stringutil
use m_diagnostic_utilities
use m_node_conf

use m_space
use m_grid_define

use m_parameters

implicit none

private

interface write_raw
  module procedure write_raw_
end interface

public :: write_raw, write_raw_select_data, write_raw_create_file
public :: create_datasets, close_datasets

contains

!-------------------------------------------------------------------------------
! Get index of selected particles for output
!-------------------------------------------------------------------------------
subroutine write_raw_select_data( spec, t, idx, num_par )
!-------------------------------------------------------------------------------

  use m_random

  implicit none

  class( t_species ), intent(inout) :: spec
  real(p_double), intent(in) :: t
  integer, dimension(:), intent(out) :: idx
  integer, intent(out) :: num_par

  ! local variables
  real(p_k_part) :: gamma_limit_sqr, gamma_par_sqr, rnd
  real(p_k_part), dimension(4) :: pos ! was 3
  real(p_k_fparse) :: raw_eval
  real(p_k_fparse), dimension(p_p_dim + spec%get_n_x_dims() + 2) :: raw_var
  logical :: has_raw_math_expr
  integer :: l, lbuf, n_x_dim

  n_x_dim = spec%get_n_x_dims()

  ! initialize gamma limit
  gamma_limit_sqr = spec%diag%raw_gamma_limit**2

  ! loop over all particles

  l = 1
  lbuf = 0

  ! Copy all the particles to the temp buffer

  ! time for raw_math evaluation
  raw_var(n_x_dim + p_p_dim + 2) = t

  has_raw_math_expr = (spec%diag%raw_math_expr /= '')

  do
     if (l > spec%num_par) exit
     gamma_par_sqr = 1.0_p_k_part + spec%p(1,l)**2 + spec%p(2,l)**2 + spec%p(3,l)**2

     if ( gamma_par_sqr >= gamma_limit_sqr ) then

        raw_eval = 1
        if (has_raw_math_expr) then

           ! fill evaluation variables
           call spec % get_position( l, pos ) ! x1-x3
           raw_var(1:n_x_dim) = real( pos(1:n_x_dim), p_k_fparse )
           raw_var(n_x_dim+1) = spec%p(1,l) ! p1
           raw_var(n_x_dim+2) = spec%p(2,l) ! p2
           raw_var(n_x_dim+3) = spec%p(3,l) ! p3
           raw_var(n_x_dim+p_p_dim + 1) = sqrt(gamma_par_sqr) ! g

           ! evaluate
           raw_eval = eval( spec%diag%raw_func, raw_var )

        endif

        if (raw_eval > 0) then

              call rng % harvest_real2( rnd, spec%ix(:,l) )
              if ( rnd <= spec%diag%raw_fraction ) then

                 lbuf = lbuf + 1
                 idx(lbuf) = l

              endif

         endif
     endif
     l = l+1
  enddo

  num_par = lbuf

end subroutine write_raw_select_data
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
subroutine write_raw_create_file( spec, total_par, n, t, n_aux, diagFile, comm )
!-------------------------------------------------------------------------------

  implicit none

  class( t_species ), intent(inout) :: spec
  integer, intent(in) :: n
  integer, intent(in) :: total_par
  real(p_double), intent(in) :: t
  integer, intent(in) :: n_aux
  integer, intent(in) :: comm

  class( t_diag_file ), allocatable, intent(inout) :: diagFile

  ! local variables
  integer :: i, j, n_x_dim

  n_x_dim = spec%get_n_x_dims()
  ! create file

  ! Prepare path and file name
  call create_diag_file( diagFile )
  diagFile % ftype = p_diag_particles

  diagFile%filepath = trim(path_mass)//'RAW'//p_dir_sep//replace_blanks(trim(spec%name))//p_dir_sep
  diagFile%filename = get_filename( n_aux, '' // 'RAW'  // '-' // replace_blanks(trim(spec%name)) )

  diagFile%name = spec%name

  diagFile%iter%n = n
  diagFile%iter%t = t
  diagFile%iter%time_units = '1 / \omega_p'

  diagFile%particles%name = spec%name

  ! Species object doesn't have a label (yet) so just set it to the name
  diagFile%particles%label = spec%name

  diagFile%particles%np = total_par
  diagFile%particles%has_tags = spec%add_tag

  j = 0

  ! positions
  do i = 1, n_x_dim
    j = j+1
    diagFile%particles%quants(j) = 'x'//(char(iachar('0')+i))
    diagFile%particles%qlabels(j) = 'x_'//(char(iachar('0')+i))
    diagFile%particles%qunits(j) = 'c/\omega_p'
    diagFile%particles%offset_t(j) = 0.0_p_double
  enddo

  ! momenta
  do i = 1, p_p_dim
    j = j+1
    diagFile%particles%quants(j) = 'p'//(char(iachar('0')+i))
    diagFile%particles%qlabels(j) = 'p_'//(char(iachar('0')+i))
    diagFile%particles%qunits(j) = 'm_e c'
    diagFile%particles%offset_t(j) = -0.5_p_double
  enddo

  ! charge
  j = j+1
  diagFile%particles%quants(j) = 'q'
  diagFile%particles%qlabels(j) = 'q'
  diagFile%particles%qunits(j) = 'e'
  diagFile%particles%offset_t(j) = 0.0_p_double

  ! energy
  j = j+1
  diagFile%particles%quants(j) = 'ene'
  diagFile%particles%qlabels(j) = 'Ene'
  diagFile%particles%qunits(j) = 'm_e c^2'
  diagFile%particles%offset_t(j) = -0.5_p_double
# 205 "spec/diagnostics/os-spec-rawdiag.f03"
  if ( diagFile%particles%has_tags ) then
    j = j+1
    diagFile%particles%quants(j) = 'tag'
    diagFile%particles%qlabels(j) = 'Tag'
    diagFile%particles%qunits(j) = ''
    diagFile%particles%offset_t(j) = 0.0_p_double
  endif

  diagFile%particles%nquants = j

  ! create the file
  call diagFile % open( p_diag_create, comm )

end subroutine write_raw_create_file
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! Create datasets in the file
!-------------------------------------------------------------------------------
subroutine create_datasets( diagFile, total_par, datasets )
!-------------------------------------------------------------------------------

  implicit none

  class( t_diag_file ), intent(inout) :: diagFile
  integer, intent(in) :: total_par

  type( t_diag_dataset ), dimension(:), intent(inout) :: datasets

  integer :: i, j, data_type

  select case (p_diag_prec)
  case ( p_single )
    data_type = p_diag_float32
  case ( p_double )
    data_type = p_diag_float64
  end select

  if ( diagFile%particles%has_tags ) then
    j = diagFile % particles % nquants - 1
  else
    j = diagFile % particles % nquants
  endif

  ! particle quantities
  do i = 1, j
    call diagFile % start_cdset( diagFile % particles % quants(i), 1, [total_par], &
                            data_type, datasets(i) )
  enddo

  ! particle tags
  if ( diagFile%particles%has_tags ) then
    i = diagFile % particles % nquants
    call diagFile % start_cdset( 'tag', 2, [2,total_par], p_diag_int32, &
                            datasets(i) )
  endif

end subroutine create_datasets
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Close datasets in the file
!-------------------------------------------------------------------------------
subroutine close_datasets( diagFile, datasets )
!-------------------------------------------------------------------------------
  implicit none

  class( t_diag_file ), intent(inout) :: diagFile
  type( t_diag_dataset ), dimension(:), intent(inout) :: datasets

  integer :: i

  do i = 1, diagFile % particles % nquants
    call diagFile % end_cdset( datasets(i) )
  enddo

end subroutine close_datasets
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
! write raw particle data
!-------------------------------------------------------------------------------

subroutine write_raw_( spec, no_co, n, t, n_aux )

  use m_units
  !use mpi

  implicit none

  class( t_species ), intent(inout) :: spec
  class( t_node_conf ), intent(in) :: no_co
  integer, intent(in) :: n
  real(p_double), intent(in) :: t
  integer, intent(in) :: n_aux

  class( t_diag_file ), allocatable :: diagFile



  type( t_diag_dataset ), dimension( spec%get_n_x_dims() + p_p_dim + 3 ) :: datasets

  type( t_diag_chunk ) :: chunk, chunk_tag
  integer :: start, rank, num_par, total_par, color, comm, ierr
  integer :: ndims, i, j, k
  integer, dimension(:), pointer :: idx
  real( p_diag_prec ), dimension(:), pointer :: write_buffer
  integer, dimension( :, : ), pointer :: write_tag_buffer
  real( p_k_part ) :: u2

  ndims = spec%get_n_x_dims()

  ! get index of particles to output
  if ( spec%num_par > 0 ) then
    call alloc(idx, [spec%num_par],"spec/diagnostics/os-spec-rawdiag.f03",322)
    call write_raw_select_data( spec, t, idx, num_par )
  else
    num_par = 0
  endif


  if ( no_co % comm_size() > 1 ) then
    comm = MPI_COMM_NULL

    call MPI_ALLREDUCE( num_par, total_par, 1, MPI_INTEGER, MPI_SUM, no_co % comm, ierr )

    if ( total_par > 0 ) then
      if ( num_par > 0 ) then
        color = 1
      else
        color = MPI_UNDEFINED
      endif
      call MPI_COMM_SPLIT( no_co % comm, color, 0, comm, ierr )

      if ( num_par > 0 ) then
        call MPI_COMM_RANK( comm, rank, ierr )
        call MPI_EXSCAN( num_par, start, 1, MPI_INTEGER, MPI_SUM, comm, ierr )
        if ( rank == 0 ) start = 0
      endif
    endif
  else
    comm = MPI_COMM_SELF
    total_par = num_par
    start = 0
  endif

  ! Write data, if any
  if ( total_par > 0 ) then

    ! Only nodes with particles are involved in this section
    if ( num_par > 0 ) then

      ! Create the file
      call write_raw_create_file( spec, total_par, n, t, n_aux, diagFile, comm )

      ! Create the datasets
      call create_datasets( diagFile, total_par, datasets )

      call alloc(write_buffer, [ num_par ],"spec/diagnostics/os-spec-rawdiag.f03",366)
      chunk % start(1) = start
      chunk % count(1) = num_par
      chunk % stride(1) = 1
      chunk % data = c_loc( write_buffer )

      ! Positions
      k = 1
      do j = 1, ndims
        call spec % get_position( j, idx, num_par, write_buffer, spec%diag%raw_ref_box(j))
        call diagFile % write_par_cdset( datasets(k), chunk, start )
        k = k+1
      enddo

      ! Momenta
      do j = 1, p_p_dim
        do i = 1, num_par
          write_buffer(i) = real(spec%p(j,idx(i)), p_diag_prec)
        enddo
        call diagFile % write_par_cdset( datasets(k), chunk, start )
        k = k+1
      enddo

      ! Charge
      do i = 1, num_par
        write_buffer(i) = real(spec%q(idx(i)), p_diag_prec)
      enddo
      call diagFile % write_par_cdset( datasets(k), chunk, start )
      k = k+1

      ! Energy
      do i = 1, num_par
        u2 = spec%p( 1, idx(i) )**2 + spec%p( 2, idx(i) )**2 + spec%p( 3, idx(i) )**2
        write_buffer(i) = real(u2 / ( sqrt(1 + u2) + 1) , p_diag_prec)
      enddo
      call diagFile % write_par_cdset( datasets(k), chunk, start )
      k = k+1
# 417 "spec/diagnostics/os-spec-rawdiag.f03"
      call freemem(write_buffer,"spec/diagnostics/os-spec-rawdiag.f03",417)

      ! Tags
      if ( spec%add_tag ) then
        call alloc(write_tag_buffer, [2, num_par],"spec/diagnostics/os-spec-rawdiag.f03",421)
        chunk_tag % count(1) = 2
        chunk_tag % count(2) = num_par
        chunk_tag % start(1) = 0
        chunk_tag % start(2) = start
        chunk_tag % stride(1) = 1
        chunk_tag % stride(2) = 1
        chunk_tag % data = c_loc(write_tag_buffer)

        do i = 1, num_par
          write_tag_buffer(1,i) = spec%tag(1,idx(i))
          write_tag_buffer(2,i) = spec%tag(2,idx(i))
        enddo

        call diagFile % write_par_cdset( datasets(k), chunk_tag, 2*start )

        call freemem(write_tag_buffer,"spec/diagnostics/os-spec-rawdiag.f03",437)

      endif

      call close_datasets( diagFile, datasets )

      ! close file
      call diagFile % close( )
      deallocate( diagFile )
    endif
  else
    ! No particles, root node writes an empty file
    if ( no_co % local_rank() == 0 ) then
      call write_raw_create_file( spec, 0, n, t, n_aux, diagFile, MPI_COMM_NULL )
      call create_datasets( diagFile, num_par, datasets )
      call close_datasets( diagFile, datasets )
      call diagFile % close( )
      deallocate( diagFile )
    endif
  endif

  ! Free particle indexes
  if ( spec%num_par > 0 ) then
    call freemem(idx,"spec/diagnostics/os-spec-rawdiag.f03",460)
  endif

  if ( comm /= MPI_COMM_NULL .and. comm /= MPI_COMM_SELF ) then
    call MPI_COMM_FREE( comm, ierr )
  endif

end subroutine write_raw_


end module m_species_rawdiag
