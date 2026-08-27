# 1 "fft/os-fft.f03"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "fft/os-fft.f03"
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
# 2 "fft/os-fft.f03" 2
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
# 3 "fft/os-fft.f03" 2

! Only one-dimensional fft wrappers are currently implemented in this module
module m_fft

use, intrinsic :: iso_c_binding
use m_system

implicit none

! define data type of fft/ifft wrappers
# 30 "fft/os-fft.f03"
integer, parameter :: type_real = C_DOUBLE
integer, parameter :: type_complex = C_DOUBLE_COMPLEX
# 61 "fft/os-fft.f03"
interface fft_init
  module procedure fftw_init_1d
  module procedure fftw_init2_1d
end interface

interface fft
  module procedure fftw_r2c_1d2
  module procedure fftw_r2c_1d3
end interface

interface ifft
  module procedure fftw_idft_1d0
  module procedure fftw_idft_1d1
  module procedure fftw_idft_1d2
  module procedure fftw_c2r_1d2
  module procedure fftw_c2r_1d3
end interface

interface fft_cleanup
  module procedure fftw_cleanup_1d
end interface


type(C_PTR) :: plan_fft_x, plan_fft_y, plan_fft_z, &
                plan_ifft_x, plan_ifft_y, plan_ifft_z, &
                ptr_real_x, ptr_real_y, ptr_real_z, &
                ptr_complex_x, ptr_complex_y, ptr_complex_z

real(type_real), dimension(:), pointer :: fftw_real_x, fftw_real_y, fftw_real_z
complex(type_complex), dimension(:), pointer :: fftw_complex_x, fftw_complex_y, fftw_complex_z

private
public :: fft_init, fft, ifft, fft_cleanup, type_real, type_complex

save

contains
# 555 "fft/os-fft.f03"
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! NOTE: The following 6 routines are a stopgap fix to allow this code to compile
! on systems without FFTW. It will be removed once a permanent FFT
! library in introduced.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

subroutine err_out( name )
  implicit none
  character(*) :: name

  write(0,*) "An FFT routine was called but the Osiris was not compiled with FFT "
  write(0,*) "support via the FFTW library. "
  write(0,*) "Please see the FFTW_DIR parameter in your system's configuration file. "
  write(0,*) "The routine called was: ", trim(name)
  call abort_program( )

end subroutine
!-------------------------------------------------------------------------------
subroutine fftw_init_1d( n )
!-------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  implicit none

  ! dummy variables
  integer, intent(in) :: n

  call err_out("fftw_init_1d")
end subroutine
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine fftw_init2_1d( n, is_backward )
!-------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  implicit none

  ! dummy variables
  integer, intent(in) :: n
  logical, intent(in) :: is_backward

  call err_out("fftw_init2_1d")
end subroutine fftw_init2_1d
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine fftw_idft_1d0( f, n )
!-------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  implicit none

  real(p_double), intent(inout), dimension(:, :) :: f
  integer, intent(in) :: n

  call err_out("fftw_idft_1d0")
end subroutine fftw_idft_1d0
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine fftw_idft_1d1( f, n )
!-------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  implicit none

  real(p_double), intent(inout), dimension(:, :, :) :: f
  integer, dimension(2), intent(in) :: n

  call err_out("fftw_idft_1d1")
end subroutine fftw_idft_1d1
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine fftw_idft_1d2( f, n )
!-------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  implicit none

  real(p_double), intent(inout), dimension(:, :, :, :) :: f
  integer, dimension(2), intent(in) :: n

  call err_out("fftw_idft_1d2")
end subroutine fftw_idft_1d2
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine fftw_r2c_1d2( f, n1, n2 )
!-------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  implicit none

  ! dummy variables
  real(type_real), intent(inout), dimension(:,:,:), pointer :: f
  integer, intent(in) :: n1, n2

  call err_out("fftw_r2c_1d2")
end subroutine
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine fftw_c2r_1d2( f, n1, n2 )
!-------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  implicit none

  ! dummy variables
  real(type_real), intent(inout), dimension(:,:,:), pointer :: f
  integer, intent(in) :: n1, n2

  call err_out("fftw_c2r_1d2")
end subroutine
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine fftw_r2c_1d3( f, n1, n2, n3 )
!-------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  implicit none

  ! dummy variables
  real(type_real), intent(inout), dimension(:,:,:, :), pointer :: f
  integer, intent(in) :: n1, n2, n3

  call err_out("fftw_r2c_1d3")
end subroutine
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine fftw_c2r_1d3( f, n1, n2, n3 )
!-------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  implicit none

  ! dummy variables
  real(type_real), intent(inout), dimension(:,:,:, :), pointer :: f
  integer, intent(in) :: n1, n2, n3

  call err_out("fftw_c2r_1d3")
end subroutine
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine fftw_cleanup_1d( )
!-------------------------------------------------------------------------------
  use, intrinsic :: iso_c_binding
  implicit none

  call err_out("fftw_cleanup_1d")
end subroutine
!-------------------------------------------------------------------------------



end module m_fft
