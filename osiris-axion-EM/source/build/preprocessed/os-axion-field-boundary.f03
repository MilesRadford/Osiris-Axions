# 1 "axion_field/os-axion-field-boundary.f03"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "axion_field/os-axion-field-boundary.f03"
module m_axion_field_boundary

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
# 4 "axion_field/os-axion-field-boundary.f03" 2
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
# 5 "axion_field/os-axion-field-boundary.f03" 2

  use m_system
  use m_parameters
  use m_vdf_define
  use m_vdf
  use m_axion_field_define, only : t_axion_field

  implicit none
  private

  public :: canonical_axion_bc
  public :: apply_axion_bc_scalar, apply_axion_bc_vector
  public :: periodic_fill_axion_field, periodic_fill_axion_vecfield

contains

! Normalize user boundary-condition strings to the internal axion-field boundary names.
character(len=32) function canonical_axion_bc( bc_in ) result( bc_out )
  character(len=*), intent(in) :: bc_in
  character(len=32) :: tmp
  integer :: i, ic

  tmp = '                                '
  tmp = adjustl( bc_in )
  do i = 1, len(tmp)
    ic = iachar( tmp(i:i) )
    if ( ic >= iachar('A') .and. ic <= iachar('Z') ) tmp(i:i) = achar( ic + 32 )
  enddo

  select case ( trim(tmp) )
  case ( '', 'periodic' )
    bc_out = 'periodic'
  case ( 'reflect', 'reflecting', 'mirror', 'neumann' )
    bc_out = 'reflect'
  case ( 'open', 'outflow', 'vpml', 'absorbing', 'absorb', 'sponge' )
    bc_out = 'open'
  case default
    write(err_buf__,*) 'Unsupported axion boundary condition = ', trim(bc_in);call err__("axion_field/os-axion-field-boundary.f03",42)
    print *, 'Valid axion boundary conditions are: periodic, reflect, open(outflow/vpml)'
    call abort_program( p_err_invalid )
  end select
end function canonical_axion_bc

! Apply the configured scalar-field boundary conditions on all active dimensions.
subroutine apply_axion_bc_scalar( this, fld )
  class(t_axion_field), intent(in) :: this
  type(t_vdf), intent(inout) :: fld

  call apply_axion_bc_scalar_dim1( fld, this%nx, this%gc_num, this%bc_type(p_lower,1), this%bc_type(p_upper,1) )
  if ( p_x_dim >= 2 ) call apply_axion_bc_scalar_dim2( fld, this%nx, this%gc_num, this%bc_type(p_lower,2), this%bc_type(p_upper,2) )
  if ( p_x_dim >= 3 ) call apply_axion_bc_scalar_dim3( fld, this%nx, this%gc_num, this%bc_type(p_lower,3), this%bc_type(p_upper,3) )
end subroutine apply_axion_bc_scalar

! Apply the configured vector-field boundary conditions on all active dimensions.
subroutine apply_axion_bc_vector( this, fld )
  class(t_axion_field), intent(in) :: this
  type(t_vdf), intent(inout) :: fld

  call apply_axion_bc_vector_dim1( fld, this%nx, this%gc_num, this%bc_type(p_lower,1), this%bc_type(p_upper,1) )
  if ( p_x_dim >= 2 ) call apply_axion_bc_vector_dim2( fld, this%nx, this%gc_num, this%bc_type(p_lower,2), this%bc_type(p_upper,2) )
  if ( p_x_dim >= 3 ) call apply_axion_bc_vector_dim3( fld, this%nx, this%gc_num, this%bc_type(p_lower,3), this%bc_type(p_upper,3) )
end subroutine apply_axion_bc_vector

! Fill scalar-field guard cells along x1 according to the selected boundary condition.
subroutine apply_axion_bc_scalar_dim1( fld, nx, gc_num, bc_lower, bc_upper )
  type(t_vdf), intent(inout) :: fld
  integer, dimension(p_max_dim), intent(in) :: nx
  integer, dimension(2,p_max_dim), intent(in) :: gc_num
  character(len=*), intent(in) :: bc_lower, bc_upper

  integer :: i1, i2, i3, gcl, gcu

  gcl = gc_num(p_lower,1)
  gcu = gc_num(p_upper,1)

  select case ( p_x_dim )
  case (1)
    do i1 = 1, gcl
      select case ( trim(bc_lower) )
      case ( 'periodic' )
        fld%f1(1,1-i1) = fld%f1(1,nx(1)-i1+1)
      case ( 'reflect' )
        fld%f1(1,1-i1) = fld%f1(1,i1)
      case ( 'open' )
        fld%f1(1,1-i1) = 2.0_p_k_fld * fld%f1(1,2-i1) - fld%f1(1,3-i1)
      end select
    enddo
    do i1 = 1, gcu
      select case ( trim(bc_upper) )
      case ( 'periodic' )
        fld%f1(1,nx(1)+i1) = fld%f1(1,i1)
      case ( 'reflect' )
        fld%f1(1,nx(1)+i1) = fld%f1(1,nx(1)-i1+1)
      case ( 'open' )
        fld%f1(1,nx(1)+i1) = 2.0_p_k_fld * fld%f1(1,nx(1)+i1-1) - fld%f1(1,nx(1)+i1-2)
      end select
    enddo
  case (2)
    do i2 = 1, nx(2)
      do i1 = 1, gcl
        select case ( trim(bc_lower) )
        case ( 'periodic' )
          fld%f2(1,1-i1,i2) = fld%f2(1,nx(1)-i1+1,i2)
        case ( 'reflect' )
          fld%f2(1,1-i1,i2) = fld%f2(1,i1,i2)
        case ( 'open' )
          fld%f2(1,1-i1,i2) = 2.0_p_k_fld * fld%f2(1,2-i1,i2) - fld%f2(1,3-i1,i2)
        end select
      enddo
      do i1 = 1, gcu
        select case ( trim(bc_upper) )
        case ( 'periodic' )
          fld%f2(1,nx(1)+i1,i2) = fld%f2(1,i1,i2)
        case ( 'reflect' )
          fld%f2(1,nx(1)+i1,i2) = fld%f2(1,nx(1)-i1+1,i2)
        case ( 'open' )
          fld%f2(1,nx(1)+i1,i2) = 2.0_p_k_fld * fld%f2(1,nx(1)+i1-1,i2) - fld%f2(1,nx(1)+i1-2,i2)
        end select
      enddo
    enddo
  case (3)
    do i3 = 1, nx(3)
      do i2 = 1, nx(2)
        do i1 = 1, gcl
          select case ( trim(bc_lower) )
          case ( 'periodic' )
            fld%f3(1,1-i1,i2,i3) = fld%f3(1,nx(1)-i1+1,i2,i3)
          case ( 'reflect' )
            fld%f3(1,1-i1,i2,i3) = fld%f3(1,i1,i2,i3)
          case ( 'open' )
            fld%f3(1,1-i1,i2,i3) = 2.0_p_k_fld * fld%f3(1,2-i1,i2,i3) - fld%f3(1,3-i1,i2,i3)
          end select
        enddo
        do i1 = 1, gcu
          select case ( trim(bc_upper) )
          case ( 'periodic' )
            fld%f3(1,nx(1)+i1,i2,i3) = fld%f3(1,i1,i2,i3)
          case ( 'reflect' )
            fld%f3(1,nx(1)+i1,i2,i3) = fld%f3(1,nx(1)-i1+1,i2,i3)
          case ( 'open' )
            fld%f3(1,nx(1)+i1,i2,i3) = 2.0_p_k_fld * fld%f3(1,nx(1)+i1-1,i2,i3) - fld%f3(1,nx(1)+i1-2,i2,i3)
          end select
        enddo
      enddo
    enddo
  end select
end subroutine apply_axion_bc_scalar_dim1

! Fill scalar-field guard cells along x2 according to the selected boundary condition.
subroutine apply_axion_bc_scalar_dim2( fld, nx, gc_num, bc_lower, bc_upper )
  type(t_vdf), intent(inout) :: fld
  integer, dimension(p_max_dim), intent(in) :: nx
  integer, dimension(2,p_max_dim), intent(in) :: gc_num
  character(len=*), intent(in) :: bc_lower, bc_upper

  integer :: i1, i2, i3, gcl, gcu
  integer :: l1, u1

  gcl = gc_num(p_lower,2)
  gcu = gc_num(p_upper,2)

  select case ( p_x_dim )
  case (2)
    l1 = lbound(fld%f2,2); u1 = ubound(fld%f2,2)
    do i1 = l1, u1
      do i2 = 1, gcl
        select case ( trim(bc_lower) )
        case ( 'periodic' )
          fld%f2(1,i1,1-i2) = fld%f2(1,i1,nx(2)-i2+1)
        case ( 'reflect' )
          fld%f2(1,i1,1-i2) = fld%f2(1,i1,i2)
        case ( 'open' )
          fld%f2(1,i1,1-i2) = 2.0_p_k_fld * fld%f2(1,i1,2-i2) - fld%f2(1,i1,3-i2)
        end select
      enddo
      do i2 = 1, gcu
        select case ( trim(bc_upper) )
        case ( 'periodic' )
          fld%f2(1,i1,nx(2)+i2) = fld%f2(1,i1,i2)
        case ( 'reflect' )
          fld%f2(1,i1,nx(2)+i2) = fld%f2(1,i1,nx(2)-i2+1)
        case ( 'open' )
          fld%f2(1,i1,nx(2)+i2) = 2.0_p_k_fld * fld%f2(1,i1,nx(2)+i2-1) - fld%f2(1,i1,nx(2)+i2-2)
        end select
      enddo
    enddo
  case (3)
    l1 = lbound(fld%f3,2); u1 = ubound(fld%f3,2)
    do i3 = 1, nx(3)
      do i1 = l1, u1
        do i2 = 1, gcl
          select case ( trim(bc_lower) )
          case ( 'periodic' )
            fld%f3(1,i1,1-i2,i3) = fld%f3(1,i1,nx(2)-i2+1,i3)
          case ( 'reflect' )
            fld%f3(1,i1,1-i2,i3) = fld%f3(1,i1,i2,i3)
          case ( 'open' )
            fld%f3(1,i1,1-i2,i3) = 2.0_p_k_fld * fld%f3(1,i1,2-i2,i3) - fld%f3(1,i1,3-i2,i3)
          end select
        enddo
        do i2 = 1, gcu
          select case ( trim(bc_upper) )
          case ( 'periodic' )
            fld%f3(1,i1,nx(2)+i2,i3) = fld%f3(1,i1,i2,i3)
          case ( 'reflect' )
            fld%f3(1,i1,nx(2)+i2,i3) = fld%f3(1,i1,nx(2)-i2+1,i3)
          case ( 'open' )
            fld%f3(1,i1,nx(2)+i2,i3) = 2.0_p_k_fld * fld%f3(1,i1,nx(2)+i2-1,i3) - fld%f3(1,i1,nx(2)+i2-2,i3)
          end select
        enddo
      enddo
    enddo
  end select
end subroutine apply_axion_bc_scalar_dim2

! Fill scalar-field guard cells along x3 according to the selected boundary condition.
subroutine apply_axion_bc_scalar_dim3( fld, nx, gc_num, bc_lower, bc_upper )
  type(t_vdf), intent(inout) :: fld
  integer, dimension(p_max_dim), intent(in) :: nx
  integer, dimension(2,p_max_dim), intent(in) :: gc_num
  character(len=*), intent(in) :: bc_lower, bc_upper

  integer :: i1, i2, i3, gcl, gcu
  integer :: l1, u1, l2, u2

  if ( p_x_dim /= 3 ) return

  gcl = gc_num(p_lower,3)
  gcu = gc_num(p_upper,3)
  l1 = lbound(fld%f3,2); u1 = ubound(fld%f3,2)
  l2 = lbound(fld%f3,3); u2 = ubound(fld%f3,3)

  do i2 = l2, u2
    do i1 = l1, u1
      do i3 = 1, gcl
        select case ( trim(bc_lower) )
        case ( 'periodic' )
          fld%f3(1,i1,i2,1-i3) = fld%f3(1,i1,i2,nx(3)-i3+1)
        case ( 'reflect' )
          fld%f3(1,i1,i2,1-i3) = fld%f3(1,i1,i2,i3)
        case ( 'open' )
          fld%f3(1,i1,i2,1-i3) = 2.0_p_k_fld * fld%f3(1,i1,i2,2-i3) - fld%f3(1,i1,i2,3-i3)
        end select
      enddo
      do i3 = 1, gcu
        select case ( trim(bc_upper) )
        case ( 'periodic' )
          fld%f3(1,i1,i2,nx(3)+i3) = fld%f3(1,i1,i2,i3)
        case ( 'reflect' )
          fld%f3(1,i1,i2,nx(3)+i3) = fld%f3(1,i1,i2,nx(3)-i3+1)
        case ( 'open' )
          fld%f3(1,i1,i2,nx(3)+i3) = 2.0_p_k_fld * fld%f3(1,i1,i2,nx(3)+i3-1) - fld%f3(1,i1,i2,nx(3)+i3-2)
        end select
      enddo
    enddo
  enddo
end subroutine apply_axion_bc_scalar_dim3

! Fill vector-field guard cells along x1 according to the selected boundary condition.
subroutine apply_axion_bc_vector_dim1( fld, nx, gc_num, bc_lower, bc_upper )
  type(t_vdf), intent(inout) :: fld
  integer, dimension(p_max_dim), intent(in) :: nx
  integer, dimension(2,p_max_dim), intent(in) :: gc_num
  character(len=*), intent(in) :: bc_lower, bc_upper

  integer :: i1, i2, i3, ic, gcl, gcu

  gcl = gc_num(p_lower,1)
  gcu = gc_num(p_upper,1)

  select case ( p_x_dim )
  case (1)
    do ic = 1, 3
      do i1 = 1, gcl
        select case ( trim(bc_lower) )
        case ( 'periodic' )
          fld%f1(ic,1-i1) = fld%f1(ic,nx(1)-i1+1)
        case ( 'reflect' )
          fld%f1(ic,1-i1) = fld%f1(ic,i1)
        case ( 'open' )
          fld%f1(ic,1-i1) = 2.0_p_k_fld * fld%f1(ic,2-i1) - fld%f1(ic,3-i1)
        end select
      enddo
      do i1 = 1, gcu
        select case ( trim(bc_upper) )
        case ( 'periodic' )
          fld%f1(ic,nx(1)+i1) = fld%f1(ic,i1)
        case ( 'reflect' )
          fld%f1(ic,nx(1)+i1) = fld%f1(ic,nx(1)-i1+1)
        case ( 'open' )
          fld%f1(ic,nx(1)+i1) = 2.0_p_k_fld * fld%f1(ic,nx(1)+i1-1) - fld%f1(ic,nx(1)+i1-2)
        end select
      enddo
    enddo
  case (2)
    do ic = 1, 3
      do i2 = 1, nx(2)
        do i1 = 1, gcl
          select case ( trim(bc_lower) )
          case ( 'periodic' )
            fld%f2(ic,1-i1,i2) = fld%f2(ic,nx(1)-i1+1,i2)
          case ( 'reflect' )
            fld%f2(ic,1-i1,i2) = fld%f2(ic,i1,i2)
          case ( 'open' )
            fld%f2(ic,1-i1,i2) = 2.0_p_k_fld * fld%f2(ic,2-i1,i2) - fld%f2(ic,3-i1,i2)
          end select
        enddo
        do i1 = 1, gcu
          select case ( trim(bc_upper) )
          case ( 'periodic' )
            fld%f2(ic,nx(1)+i1,i2) = fld%f2(ic,i1,i2)
          case ( 'reflect' )
            fld%f2(ic,nx(1)+i1,i2) = fld%f2(ic,nx(1)-i1+1,i2)
          case ( 'open' )
            fld%f2(ic,nx(1)+i1,i2) = 2.0_p_k_fld * fld%f2(ic,nx(1)+i1-1,i2) - fld%f2(ic,nx(1)+i1-2,i2)
          end select
        enddo
      enddo
    enddo
  case (3)
    do ic = 1, 3
      do i3 = 1, nx(3)
        do i2 = 1, nx(2)
          do i1 = 1, gcl
            select case ( trim(bc_lower) )
            case ( 'periodic' )
              fld%f3(ic,1-i1,i2,i3) = fld%f3(ic,nx(1)-i1+1,i2,i3)
            case ( 'reflect' )
              fld%f3(ic,1-i1,i2,i3) = fld%f3(ic,i1,i2,i3)
            case ( 'open' )
              fld%f3(ic,1-i1,i2,i3) = 2.0_p_k_fld * fld%f3(ic,2-i1,i2,i3) - fld%f3(ic,3-i1,i2,i3)
            end select
          enddo
          do i1 = 1, gcu
            select case ( trim(bc_upper) )
            case ( 'periodic' )
              fld%f3(ic,nx(1)+i1,i2,i3) = fld%f3(ic,i1,i2,i3)
            case ( 'reflect' )
              fld%f3(ic,nx(1)+i1,i2,i3) = fld%f3(ic,nx(1)-i1+1,i2,i3)
            case ( 'open' )
              fld%f3(ic,nx(1)+i1,i2,i3) = 2.0_p_k_fld * fld%f3(ic,nx(1)+i1-1,i2,i3) - fld%f3(ic,nx(1)+i1-2,i2,i3)
            end select
          enddo
        enddo
      enddo
    enddo
  end select
end subroutine apply_axion_bc_vector_dim1

! Fill vector-field guard cells along x2 according to the selected boundary condition.
subroutine apply_axion_bc_vector_dim2( fld, nx, gc_num, bc_lower, bc_upper )
  type(t_vdf), intent(inout) :: fld
  integer, dimension(p_max_dim), intent(in) :: nx
  integer, dimension(2,p_max_dim), intent(in) :: gc_num
  character(len=*), intent(in) :: bc_lower, bc_upper

  integer :: i1, i2, i3, ic, gcl, gcu
  integer :: l1, u1

  gcl = gc_num(p_lower,2)
  gcu = gc_num(p_upper,2)

  select case ( p_x_dim )
  case (2)
    l1 = lbound(fld%f2,2); u1 = ubound(fld%f2,2)
    do ic = 1, 3
      do i1 = l1, u1
        do i2 = 1, gcl
          select case ( trim(bc_lower) )
          case ( 'periodic' )
            fld%f2(ic,i1,1-i2) = fld%f2(ic,i1,nx(2)-i2+1)
          case ( 'reflect' )
            fld%f2(ic,i1,1-i2) = fld%f2(ic,i1,i2)
          case ( 'open' )
            fld%f2(ic,i1,1-i2) = 2.0_p_k_fld * fld%f2(ic,i1,2-i2) - fld%f2(ic,i1,3-i2)
          end select
        enddo
        do i2 = 1, gcu
          select case ( trim(bc_upper) )
          case ( 'periodic' )
            fld%f2(ic,i1,nx(2)+i2) = fld%f2(ic,i1,i2)
          case ( 'reflect' )
            fld%f2(ic,i1,nx(2)+i2) = fld%f2(ic,i1,nx(2)-i2+1)
          case ( 'open' )
            fld%f2(ic,i1,nx(2)+i2) = 2.0_p_k_fld * fld%f2(ic,i1,nx(2)+i2-1) - fld%f2(ic,i1,nx(2)+i2-2)
          end select
        enddo
      enddo
    enddo
  case (3)
    l1 = lbound(fld%f3,2); u1 = ubound(fld%f3,2)
    do ic = 1, 3
      do i3 = 1, nx(3)
        do i1 = l1, u1
          do i2 = 1, gcl
            select case ( trim(bc_lower) )
            case ( 'periodic' )
              fld%f3(ic,i1,1-i2,i3) = fld%f3(ic,i1,nx(2)-i2+1,i3)
            case ( 'reflect' )
              fld%f3(ic,i1,1-i2,i3) = fld%f3(ic,i1,i2,i3)
            case ( 'open' )
              fld%f3(ic,i1,1-i2,i3) = 2.0_p_k_fld * fld%f3(ic,i1,2-i2,i3) - fld%f3(ic,i1,3-i2,i3)
            end select
          enddo
          do i2 = 1, gcu
            select case ( trim(bc_upper) )
            case ( 'periodic' )
              fld%f3(ic,i1,nx(2)+i2,i3) = fld%f3(ic,i1,i2,i3)
            case ( 'reflect' )
              fld%f3(ic,i1,nx(2)+i2,i3) = fld%f3(ic,i1,nx(2)-i2+1,i3)
            case ( 'open' )
              fld%f3(ic,i1,nx(2)+i2,i3) = 2.0_p_k_fld * fld%f3(ic,i1,nx(2)+i2-1,i3) - fld%f3(ic,i1,nx(2)+i2-2,i3)
            end select
          enddo
        enddo
      enddo
    enddo
  end select
end subroutine apply_axion_bc_vector_dim2

! Fill vector-field guard cells along x3 according to the selected boundary condition.
subroutine apply_axion_bc_vector_dim3( fld, nx, gc_num, bc_lower, bc_upper )
  type(t_vdf), intent(inout) :: fld
  integer, dimension(p_max_dim), intent(in) :: nx
  integer, dimension(2,p_max_dim), intent(in) :: gc_num
  character(len=*), intent(in) :: bc_lower, bc_upper

  integer :: i1, i2, i3, ic, gcl, gcu
  integer :: l1, u1, l2, u2

  if ( p_x_dim /= 3 ) return

  gcl = gc_num(p_lower,3)
  gcu = gc_num(p_upper,3)
  l1 = lbound(fld%f3,2); u1 = ubound(fld%f3,2)
  l2 = lbound(fld%f3,3); u2 = ubound(fld%f3,3)

  do ic = 1, 3
    do i2 = l2, u2
      do i1 = l1, u1
        do i3 = 1, gcl
          select case ( trim(bc_lower) )
          case ( 'periodic' )
            fld%f3(ic,i1,i2,1-i3) = fld%f3(ic,i1,i2,nx(3)-i3+1)
          case ( 'reflect' )
            fld%f3(ic,i1,i2,1-i3) = fld%f3(ic,i1,i2,i3)
          case ( 'open' )
            fld%f3(ic,i1,i2,1-i3) = 2.0_p_k_fld * fld%f3(ic,i1,i2,2-i3) - fld%f3(ic,i1,i2,3-i3)
          end select
        enddo
        do i3 = 1, gcu
          select case ( trim(bc_upper) )
          case ( 'periodic' )
            fld%f3(ic,i1,i2,nx(3)+i3) = fld%f3(ic,i1,i2,i3)
          case ( 'reflect' )
            fld%f3(ic,i1,i2,nx(3)+i3) = fld%f3(ic,i1,i2,nx(3)-i3+1)
          case ( 'open' )
            fld%f3(ic,i1,i2,nx(3)+i3) = 2.0_p_k_fld * fld%f3(ic,i1,i2,nx(3)+i3-1) - fld%f3(ic,i1,i2,nx(3)+i3-2)
          end select
        enddo
      enddo
    enddo
  enddo
end subroutine apply_axion_bc_vector_dim3


! Fill scalar-field periodic guard cells after boundary communication.
subroutine periodic_fill_axion_field( fld, nx, gc_num )
  type(t_vdf), intent(inout) :: fld
  integer, dimension(p_max_dim), intent(in) :: nx
  integer, dimension(2,p_max_dim), intent(in) :: gc_num

  integer :: i1, i2, i3, gcl, gcu
  integer :: l1,u1,l2,u2,l3,u3

  select case ( p_x_dim )
  case (1)
    gcl = gc_num(p_lower,1); gcu = gc_num(p_upper,1)
    l1 = lbound(fld%f1,2); u1 = ubound(fld%f1,2)
    do i1 = 1, gcl
      fld%f1(1,1-i1) = fld%f1(1,nx(1)-i1+1)
    enddo
    do i1 = 1, gcu
      fld%f1(1,nx(1)+i1) = fld%f1(1,i1)
    enddo
  case (2)
    l1 = lbound(fld%f2,2); u1 = ubound(fld%f2,2)
    l2 = lbound(fld%f2,3); u2 = ubound(fld%f2,3)
    gcl = gc_num(p_lower,1); gcu = gc_num(p_upper,1)
    do i2 = 1, nx(2)
      do i1 = 1, gcl
        fld%f2(1,1-i1,i2) = fld%f2(1,nx(1)-i1+1,i2)
      enddo
      do i1 = 1, gcu
        fld%f2(1,nx(1)+i1,i2) = fld%f2(1,i1,i2)
      enddo
    enddo
    gcl = gc_num(p_lower,2); gcu = gc_num(p_upper,2)
    do i1 = l1, u1
      do i2 = 1, gcl
        fld%f2(1,i1,1-i2) = fld%f2(1,i1,nx(2)-i2+1)
      enddo
      do i2 = 1, gcu
        fld%f2(1,i1,nx(2)+i2) = fld%f2(1,i1,i2)
      enddo
    enddo
  case (3)
    l1 = lbound(fld%f3,2); u1 = ubound(fld%f3,2)
    l2 = lbound(fld%f3,3); u2 = ubound(fld%f3,3)
    l3 = lbound(fld%f3,4); u3 = ubound(fld%f3,4)
    gcl = gc_num(p_lower,1); gcu = gc_num(p_upper,1)
    do i3 = 1, nx(3)
      do i2 = 1, nx(2)
        do i1 = 1, gcl
          fld%f3(1,1-i1,i2,i3) = fld%f3(1,nx(1)-i1+1,i2,i3)
        enddo
        do i1 = 1, gcu
          fld%f3(1,nx(1)+i1,i2,i3) = fld%f3(1,i1,i2,i3)
        enddo
      enddo
    enddo
    gcl = gc_num(p_lower,2); gcu = gc_num(p_upper,2)
    do i3 = 1, nx(3)
      do i1 = l1, u1
        do i2 = 1, gcl
          fld%f3(1,i1,1-i2,i3) = fld%f3(1,i1,nx(2)-i2+1,i3)
        enddo
        do i2 = 1, gcu
          fld%f3(1,i1,nx(2)+i2,i3) = fld%f3(1,i1,i2,i3)
        enddo
      enddo
    enddo
    gcl = gc_num(p_lower,3); gcu = gc_num(p_upper,3)
    do i2 = l2, u2
      do i1 = l1, u1
        do i3 = 1, gcl
          fld%f3(1,i1,i2,1-i3) = fld%f3(1,i1,i2,nx(3)-i3+1)
        enddo
        do i3 = 1, gcu
          fld%f3(1,i1,i2,nx(3)+i3) = fld%f3(1,i1,i2,i3)
        enddo
      enddo
    enddo
  end select
end subroutine periodic_fill_axion_field

! Fill vector-field periodic guard cells after boundary communication.
subroutine periodic_fill_axion_vecfield( fld, nx, gc_num )
  type(t_vdf), intent(inout) :: fld
  integer, dimension(p_max_dim), intent(in) :: nx
  integer, dimension(2,p_max_dim), intent(in) :: gc_num

  integer :: i1, i2, i3, ic, gcl, gcu
  integer :: l1,u1,l2,u2

  select case ( p_x_dim )
  case (1)
    gcl = gc_num(p_lower,1); gcu = gc_num(p_upper,1)
    do ic = 1, 3
      do i1 = 1, gcl
        fld%f1(ic,1-i1) = fld%f1(ic,nx(1)-i1+1)
      enddo
      do i1 = 1, gcu
        fld%f1(ic,nx(1)+i1) = fld%f1(ic,i1)
      enddo
    enddo
  case (2)
    l1 = lbound(fld%f2,2); u1 = ubound(fld%f2,2)
    gcl = gc_num(p_lower,1); gcu = gc_num(p_upper,1)
    do ic = 1, 3
      do i2 = 1, nx(2)
        do i1 = 1, gcl
          fld%f2(ic,1-i1,i2) = fld%f2(ic,nx(1)-i1+1,i2)
        enddo
        do i1 = 1, gcu
          fld%f2(ic,nx(1)+i1,i2) = fld%f2(ic,i1,i2)
        enddo
      enddo
    enddo
    gcl = gc_num(p_lower,2); gcu = gc_num(p_upper,2)
    do ic = 1, 3
      do i1 = l1, u1
        do i2 = 1, gcl
          fld%f2(ic,i1,1-i2) = fld%f2(ic,i1,nx(2)-i2+1)
        enddo
        do i2 = 1, gcu
          fld%f2(ic,i1,nx(2)+i2) = fld%f2(ic,i1,i2)
        enddo
      enddo
    enddo
  case (3)
    l1 = lbound(fld%f3,2); u1 = ubound(fld%f3,2)
    l2 = lbound(fld%f3,3); u2 = ubound(fld%f3,3)
    gcl = gc_num(p_lower,1); gcu = gc_num(p_upper,1)
    do ic = 1, 3
      do i3 = 1, nx(3)
        do i2 = 1, nx(2)
          do i1 = 1, gcl
            fld%f3(ic,1-i1,i2,i3) = fld%f3(ic,nx(1)-i1+1,i2,i3)
          enddo
          do i1 = 1, gcu
            fld%f3(ic,nx(1)+i1,i2,i3) = fld%f3(ic,i1,i2,i3)
          enddo
        enddo
      enddo
    enddo
    gcl = gc_num(p_lower,2); gcu = gc_num(p_upper,2)
    do ic = 1, 3
      do i3 = 1, nx(3)
        do i1 = l1, u1
          do i2 = 1, gcl
            fld%f3(ic,i1,1-i2,i3) = fld%f3(ic,i1,nx(2)-i2+1,i3)
          enddo
          do i2 = 1, gcu
            fld%f3(ic,i1,nx(2)+i2,i3) = fld%f3(ic,i1,i2,i3)
          enddo
        enddo
      enddo
    enddo
    gcl = gc_num(p_lower,3); gcu = gc_num(p_upper,3)
    do ic = 1, 3
      do i2 = l2, u2
        do i1 = l1, u1
          do i3 = 1, gcl
            fld%f3(ic,i1,i2,1-i3) = fld%f3(ic,i1,i2,nx(3)-i3+1)
          enddo
          do i3 = 1, gcu
            fld%f3(ic,i1,i2,nx(3)+i3) = fld%f3(ic,i1,i2,i3)
          enddo
        enddo
      enddo
    enddo
  end select
end subroutine periodic_fill_axion_vecfield
end module m_axion_field_boundary
