# 1 "os-axion-fps-context.f03"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "os-axion-fps-context.f03"
module m_axion_fps_context

# 1 "os-config.h" 1
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
# 4 "os-axion-fps-context.f03" 2
# 1 "os-preprocess.fpp" 1
!
! File: os-preprocess.fpp
!
! A set of preprocessing macros for osiris, and the fpp/cpp preprocessor
!
!



! Macros for the IBM CPP/GNU preprocessor




! Assertion macro



! Debug functions
# 33 "os-preprocess.fpp"
! SCR hostname functions







! LOG functions




! ERROR functions






! WARNING functions




! functions for restart io
# 5 "os-axion-fps-context.f03" 2

  use m_vdf_define, only : t_vdf

  implicit none
  private

  logical, save :: fps_active = .false.
  logical, save :: fps_part_grid_center = .false.
  type(t_vdf), pointer, save :: fps_e1 => null()
  type(t_vdf), pointer, save :: fps_b1 => null()

  public :: set_axion_fps_context
  public :: clear_axion_fps_context
  public :: axion_fps_is_active
  public :: get_axion_fps_fields
  public :: axion_fps_part_grid_center

contains

subroutine set_axion_fps_context( e1, b1, part_grid_center )
  class(t_vdf), target, intent(in) :: e1, b1
  logical, intent(in), optional :: part_grid_center

  fps_e1 => e1
  fps_b1 => b1
  fps_active = .true.
  if ( present(part_grid_center) ) then
    fps_part_grid_center = part_grid_center
  else
    fps_part_grid_center = .false.
  endif
end subroutine set_axion_fps_context

subroutine clear_axion_fps_context()
  fps_active = .false.
  fps_part_grid_center = .false.
  nullify( fps_e1 )
  nullify( fps_b1 )
end subroutine clear_axion_fps_context

logical function axion_fps_is_active()
  axion_fps_is_active = fps_active .and. associated(fps_e1) .and. associated(fps_b1)
end function axion_fps_is_active

subroutine get_axion_fps_fields( e1, b1 )
  type(t_vdf), pointer, intent(out) :: e1, b1
  e1 => fps_e1
  b1 => fps_b1
end subroutine get_axion_fps_fields

logical function axion_fps_part_grid_center()
  axion_fps_part_grid_center = fps_part_grid_center
end function axion_fps_part_grid_center

end module m_axion_fps_context
