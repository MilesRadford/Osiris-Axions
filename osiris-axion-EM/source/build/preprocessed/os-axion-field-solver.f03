# 1 "axion_field/os-axion-field-solver.f03"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "axion_field/os-axion-field-solver.f03"
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
# 2 "axion_field/os-axion-field-solver.f03" 2
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
# 3 "axion_field/os-axion-field-solver.f03" 2

module m_axion_field_solver_impl

  use m_system
  use m_parameters
  use m_input_file
  use m_vdf_define
  use m_vdf
  use m_vdf_comm
  use m_space
  use m_grid_define
  use m_node_conf
  use m_time_step
  use m_emf_define, only : t_emf
  use m_vdf_report
  use m_diagnostic_utilities
  use m_logprof, only : create_event, begin_event, end_event
  use stringutil
  use m_axion_field_define
  use m_axion_field_boundary, only : apply_axion_bc_scalar, apply_axion_bc_vector

  implicit none
  private

  public :: advance_free_axion_field_impl
  public :: advance_driven_axion_field_impl
  public :: apply_backreaction_axion_field_impl
  public :: reset_perturbation_axion_field_impl
  public :: advance_perturb_b_axion_field_impl
  public :: phi_center_1d, phi_center_2d, phi_center_3d
  public :: phi_t_1d, phi_t_2d, phi_t_3d

contains

! Advance the free massive scalar field by one leapfrog step.
subroutine advance_free_axion_field_impl( this, dt, no_co, send_msg, recv_msg )
  class(t_axion_field), intent(inout) :: this
  real(p_double), intent(in) :: dt
  class(t_node_conf), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  integer :: i1, i2, i3
  real(p_double) :: lap, idx2, idy2, idz2, m2, dt2, phi_nm, phi_np, phi_new
  real(p_double) :: rhs, kx, ky, kz, omega_d, cosfac, twopi, lx, ly, lz, phi0

  if ( .not. this%if_initialized ) return

  this%edotb_max = 0.0_p_double
  this%edotb_sum = 0.0_p_double
  this%src_max = 0.0_p_double
  this%src_sum = 0.0_p_double
  this%src_count = 0

  if ( this%first_advance .and. trim(this%init_type) /= 'zero' .and. &
       trim(this%init_type) /= 'travelling_pulse' .and. trim(this%init_type) /= 'traveling_pulse' .and. &
       trim(this%init_type) /= 'gaussian_pulse' .and. trim(this%init_type) /= 'sin_band' .and. &
       trim(this%init_type) /= 'sinx_band' ) then
    twopi = 2.0_p_double * acos(-1.0_p_double)
    lx = this%xmax_loc(1) - this%xmin_loc(1)
    rhs = this%ma*this%ma
    if ( lx > 0.0_p_double ) then
      kx = twopi * real(this%init_kx,p_double) / lx
      rhs = rhs + ( 2.0_p_double / this%dx(1) * sin( 0.5_p_double * kx * this%dx(1) ) )**2
    endif
    if ( p_x_dim >= 2 ) then
      ly = this%xmax_loc(2) - this%xmin_loc(2)
      if ( ly > 0.0_p_double ) then
        ky = twopi * real(this%init_ky,p_double) / ly
        rhs = rhs + ( 2.0_p_double / this%dx(2) * sin( 0.5_p_double * ky * this%dx(2) ) )**2
      endif
    endif
    if ( p_x_dim >= 3 ) then
      lz = this%xmax_loc(3) - this%xmin_loc(3)
      if ( lz > 0.0_p_double ) then
        kz = twopi * real(this%init_kz,p_double) / lz
        rhs = rhs + ( 2.0_p_double / this%dx(3) * sin( 0.5_p_double * kz * this%dx(3) ) )**2
      endif
    endif
    rhs = max( rhs, 0.0_p_double )
    if ( dt > 0.0_p_double ) then
      omega_d = 2.0_p_double / dt * asin( min( 1.0_p_double, 0.5_p_double * dt * sqrt(rhs) ) )
      cosfac = cos( 0.5_p_double * omega_d * dt )
    else
      cosfac = 1.0_p_double
    endif
    if ( abs(cosfac) < 1.0e-12_p_double ) cosfac = 1.0_p_double

    select case ( p_x_dim )
    case (1)
      do i1 = 1, this%nx(1)
        phi0 = 0.5_p_double * ( real(this%a%f1(1,i1),p_double) + real(this%pi%f1(1,i1),p_double) )
        this%a%f1(1,i1) = real( phi0 * cosfac, p_k_fld )
        this%pi%f1(1,i1) = real( phi0 * cosfac, p_k_fld )
      enddo
    case (2)
      do i2 = 1, this%nx(2)
        do i1 = 1, this%nx(1)
          phi0 = 0.5_p_double * ( real(this%a%f2(1,i1,i2),p_double) + real(this%pi%f2(1,i1,i2),p_double) )
          this%a%f2(1,i1,i2) = real( phi0 * cosfac, p_k_fld )
          this%pi%f2(1,i1,i2) = real( phi0 * cosfac, p_k_fld )
        enddo
      enddo
    case (3)
      do i3 = 1, this%nx(3)
        do i2 = 1, this%nx(2)
          do i1 = 1, this%nx(1)
            phi0 = 0.5_p_double * ( real(this%a%f3(1,i1,i2,i3),p_double) + real(this%pi%f3(1,i1,i2,i3),p_double) )
            this%a%f3(1,i1,i2,i3) = real( phi0 * cosfac, p_k_fld )
            this%pi%f3(1,i1,i2,i3) = real( phi0 * cosfac, p_k_fld )
          enddo
        enddo
      enddo
    end select

    call apply_axion_bc_scalar( this, this%a )
    call apply_axion_bc_scalar( this, this%pi )
    this%first_advance = .false.
  endif

  m2 = this%ma * this%ma
  dt2 = dt * dt
  idx2 = 1.0_p_double / ( this%dx(1) * this%dx(1) )
  if ( p_x_dim >= 2 ) idy2 = 1.0_p_double / ( this%dx(2) * this%dx(2) )
  if ( p_x_dim >= 3 ) idz2 = 1.0_p_double / ( this%dx(3) * this%dx(3) )

  call apply_axion_bc_scalar( this, this%pi )

  this%edotb_max = 0.0_p_double
  this%edotb_sum = 0.0_p_double
  this%src_max = 0.0_p_double
  this%src_sum = 0.0_p_double
  this%src_count = 0

  select case ( p_x_dim )
  case (1)
    do i1 = 1, this%nx(1)
      phi_nm = real(this%a%f1(1,i1),p_double)
      phi_np = real(this%pi%f1(1,i1),p_double)
      lap = ( real(this%pi%f1(1,i1+1),p_double) - 2.0_p_double * phi_np + real(this%pi%f1(1,i1-1),p_double) ) * idx2
      phi_new = 2.0_p_double * phi_np - phi_nm + dt2 * ( lap - m2 * phi_np )
      this%a%f1(1,i1) = real( phi_new, p_k_fld )
    enddo
    do i1 = 1, this%nx(1)
      phi_np = real(this%pi%f1(1,i1),p_double)
      this%pi%f1(1,i1) = this%a%f1(1,i1)
      this%a%f1(1,i1) = real( phi_np, p_k_fld )
    enddo
  case (2)
    do i2 = 1, this%nx(2)
      do i1 = 1, this%nx(1)
        phi_nm = real(this%a%f2(1,i1,i2),p_double)
        phi_np = real(this%pi%f2(1,i1,i2),p_double)
        lap = ( real(this%pi%f2(1,i1+1,i2),p_double) - 2.0_p_double * phi_np + real(this%pi%f2(1,i1-1,i2),p_double) ) * idx2 + &
              ( real(this%pi%f2(1,i1,i2+1),p_double) - 2.0_p_double * phi_np + real(this%pi%f2(1,i1,i2-1),p_double) ) * idy2
        phi_new = 2.0_p_double * phi_np - phi_nm + dt2 * ( lap - m2 * phi_np )
        this%a%f2(1,i1,i2) = real( phi_new, p_k_fld )
      enddo
    enddo
    do i2 = 1, this%nx(2)
      do i1 = 1, this%nx(1)
        phi_np = real(this%pi%f2(1,i1,i2),p_double)
        this%pi%f2(1,i1,i2) = this%a%f2(1,i1,i2)
        this%a%f2(1,i1,i2) = real( phi_np, p_k_fld )
      enddo
    enddo
  case (3)
    do i3 = 1, this%nx(3)
      do i2 = 1, this%nx(2)
        do i1 = 1, this%nx(1)
          phi_nm = real(this%a%f3(1,i1,i2,i3),p_double)
          phi_np = real(this%pi%f3(1,i1,i2,i3),p_double)
          lap = ( real(this%pi%f3(1,i1+1,i2,i3),p_double) - 2.0_p_double * phi_np + real(this%pi%f3(1,i1-1,i2,i3),p_double) ) * idx2 + &
                ( real(this%pi%f3(1,i1,i2+1,i3),p_double) - 2.0_p_double * phi_np + real(this%pi%f3(1,i1,i2-1,i3),p_double) ) * idy2 + &
                ( real(this%pi%f3(1,i1,i2,i3+1),p_double) - 2.0_p_double * phi_np + real(this%pi%f3(1,i1,i2,i3-1),p_double) ) * idz2
          phi_new = 2.0_p_double * phi_np - phi_nm + dt2 * ( lap - m2 * phi_np )
          this%a%f3(1,i1,i2,i3) = real( phi_new, p_k_fld )
        enddo
      enddo
    enddo
    do i3 = 1, this%nx(3)
      do i2 = 1, this%nx(2)
        do i1 = 1, this%nx(1)
          phi_np = real(this%pi%f3(1,i1,i2,i3),p_double)
          this%pi%f3(1,i1,i2,i3) = this%a%f3(1,i1,i2,i3)
          this%a%f3(1,i1,i2,i3) = real( phi_np, p_k_fld )
        enddo
      enddo
    enddo
  end select

  call apply_axion_bc_scalar( this, this%a )
  call apply_axion_bc_scalar( this, this%pi )
end subroutine advance_free_axion_field_impl


! Advance the driven axion field using the total E dot B source term.
subroutine advance_driven_axion_field_impl( this, emf, dt, no_co, send_msg, recv_msg )
  class(t_axion_field), intent(inout) :: this
  class(t_emf), intent(in) :: emf
  real(p_double), intent(in) :: dt
  class(t_node_conf), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  integer :: i1, i2, i3
  real(p_double) :: lap, idx2, idy2, idz2, m2, dt2, phi_nm, phi_np, phi_new
  real(p_double) :: rhs, kx, ky, kz, omega_d, cosfac, twopi, lx, ly, lz, phi0
  real(p_double) :: src, edotb

  if ( .not. this%if_initialized ) return

  if ( this%first_advance .and. trim(this%init_type) /= 'zero' .and. &
       trim(this%init_type) /= 'travelling_pulse' .and. trim(this%init_type) /= 'traveling_pulse' .and. &
       trim(this%init_type) /= 'gaussian_pulse' .and. trim(this%init_type) /= 'sin_band' .and. &
       trim(this%init_type) /= 'sinx_band' ) then
    twopi = 2.0_p_double * acos(-1.0_p_double)
    lx = this%xmax_loc(1) - this%xmin_loc(1)
    rhs = this%ma*this%ma
    if ( lx > 0.0_p_double ) then
      kx = twopi * real(this%init_kx,p_double) / lx
      rhs = rhs + ( 2.0_p_double / this%dx(1) * sin( 0.5_p_double * kx * this%dx(1) ) )**2
    endif
    if ( p_x_dim >= 2 ) then
      ly = this%xmax_loc(2) - this%xmin_loc(2)
      if ( ly > 0.0_p_double ) then
        ky = twopi * real(this%init_ky,p_double) / ly
        rhs = rhs + ( 2.0_p_double / this%dx(2) * sin( 0.5_p_double * ky * this%dx(2) ) )**2
      endif
    endif
    if ( p_x_dim >= 3 ) then
      lz = this%xmax_loc(3) - this%xmin_loc(3)
      if ( lz > 0.0_p_double ) then
        kz = twopi * real(this%init_kz,p_double) / lz
        rhs = rhs + ( 2.0_p_double / this%dx(3) * sin( 0.5_p_double * kz * this%dx(3) ) )**2
      endif
    endif
    rhs = max( rhs, 0.0_p_double )
    if ( dt > 0.0_p_double ) then
      omega_d = 2.0_p_double / dt * asin( min( 1.0_p_double, 0.5_p_double * dt * sqrt(rhs) ) )
      cosfac = cos( 0.5_p_double * omega_d * dt )
    else
      cosfac = 1.0_p_double
    endif
    if ( abs(cosfac) < 1.0e-12_p_double ) cosfac = 1.0_p_double

    select case ( p_x_dim )
    case (1)
      do i1 = 1, this%nx(1)
        phi0 = 0.5_p_double * ( real(this%a%f1(1,i1),p_double) + real(this%pi%f1(1,i1),p_double) )
        this%a%f1(1,i1) = real( phi0 * cosfac, p_k_fld )
        this%pi%f1(1,i1) = real( phi0 * cosfac, p_k_fld )
      enddo
    case (2)
      do i2 = 1, this%nx(2)
        do i1 = 1, this%nx(1)
          phi0 = 0.5_p_double * ( real(this%a%f2(1,i1,i2),p_double) + real(this%pi%f2(1,i1,i2),p_double) )
          this%a%f2(1,i1,i2) = real( phi0 * cosfac, p_k_fld )
          this%pi%f2(1,i1,i2) = real( phi0 * cosfac, p_k_fld )
        enddo
      enddo
    case (3)
      do i3 = 1, this%nx(3)
        do i2 = 1, this%nx(2)
          do i1 = 1, this%nx(1)
            phi0 = 0.5_p_double * ( real(this%a%f3(1,i1,i2,i3),p_double) + real(this%pi%f3(1,i1,i2,i3),p_double) )
            this%a%f3(1,i1,i2,i3) = real( phi0 * cosfac, p_k_fld )
            this%pi%f3(1,i1,i2,i3) = real( phi0 * cosfac, p_k_fld )
          enddo
        enddo
      enddo
    end select

    call apply_axion_bc_scalar( this, this%a )
    call apply_axion_bc_scalar( this, this%pi )
    this%first_advance = .false.
  endif

  m2 = this%ma * this%ma
  dt2 = dt * dt
  idx2 = 1.0_p_double / ( this%dx(1) * this%dx(1) )
  if ( p_x_dim >= 2 ) idy2 = 1.0_p_double / ( this%dx(2) * this%dx(2) )
  if ( p_x_dim >= 3 ) idz2 = 1.0_p_double / ( this%dx(3) * this%dx(3) )

  call apply_axion_bc_scalar( this, this%pi )
  if ( this%pert%if_initialized ) then
    call apply_axion_bc_vector( this, this%pert%e )
    call apply_axion_bc_vector( this, this%pert%b )
  endif


  select case ( p_x_dim )
  case (1)
    do i1 = 1, this%nx(1)
      phi_nm = real(this%a%f1(1,i1),p_double)
      phi_np = real(this%pi%f1(1,i1),p_double)
      lap = ( real(this%pi%f1(1,i1+1),p_double) - 2.0_p_double * phi_np + real(this%pi%f1(1,i1-1),p_double) ) * idx2
      edotb = edotb_total_cell_center_1d( this, emf, i1 )
      this%diag_edotb%f1(1,i1) = real( edotb, p_k_fld )
      src = - this%g_agg * edotb
      this%edotb_max = max( this%edotb_max, abs(edotb) )
      this%edotb_sum = this%edotb_sum + edotb
      this%src_max = max( this%src_max, abs(src) )
      this%src_sum = this%src_sum + src
      this%src_count = this%src_count + 1
      phi_new = 2.0_p_double * phi_np - phi_nm + dt2 * ( lap - m2 * phi_np + src )
      this%a%f1(1,i1) = real( phi_new, p_k_fld )
    enddo
    do i1 = 1, this%nx(1)
      phi_np = real(this%pi%f1(1,i1),p_double)
      this%pi%f1(1,i1) = this%a%f1(1,i1)
      this%a%f1(1,i1) = real( phi_np, p_k_fld )
    enddo
  case (2)
    do i2 = 1, this%nx(2)
      do i1 = 1, this%nx(1)
        phi_nm = real(this%a%f2(1,i1,i2),p_double)
        phi_np = real(this%pi%f2(1,i1,i2),p_double)
        lap = ( real(this%pi%f2(1,i1+1,i2),p_double) - 2.0_p_double * phi_np + real(this%pi%f2(1,i1-1,i2),p_double) ) * idx2 + &
              ( real(this%pi%f2(1,i1,i2+1),p_double) - 2.0_p_double * phi_np + real(this%pi%f2(1,i1,i2-1),p_double) ) * idy2
        edotb = edotb_total_cell_center_2d( this, emf, i1, i2 )
        this%diag_edotb%f2(1,i1,i2) = real( edotb, p_k_fld )
        src = - this%g_agg * edotb
        this%edotb_max = max( this%edotb_max, abs(edotb) )
        this%edotb_sum = this%edotb_sum + edotb
        this%src_max = max( this%src_max, abs(src) )
        this%src_sum = this%src_sum + src
        this%src_count = this%src_count + 1
        phi_new = 2.0_p_double * phi_np - phi_nm + dt2 * ( lap - m2 * phi_np + src )
        this%a%f2(1,i1,i2) = real( phi_new, p_k_fld )
      enddo
    enddo
    do i2 = 1, this%nx(2)
      do i1 = 1, this%nx(1)
        phi_np = real(this%pi%f2(1,i1,i2),p_double)
        this%pi%f2(1,i1,i2) = this%a%f2(1,i1,i2)
        this%a%f2(1,i1,i2) = real( phi_np, p_k_fld )
      enddo
    enddo
  case (3)
    do i3 = 1, this%nx(3)
      do i2 = 1, this%nx(2)
        do i1 = 1, this%nx(1)
          phi_nm = real(this%a%f3(1,i1,i2,i3),p_double)
          phi_np = real(this%pi%f3(1,i1,i2,i3),p_double)
          lap = ( real(this%pi%f3(1,i1+1,i2,i3),p_double) - 2.0_p_double * phi_np + real(this%pi%f3(1,i1-1,i2,i3),p_double) ) * idx2 + &
                ( real(this%pi%f3(1,i1,i2+1,i3),p_double) - 2.0_p_double * phi_np + real(this%pi%f3(1,i1,i2-1,i3),p_double) ) * idy2 + &
                ( real(this%pi%f3(1,i1,i2,i3+1),p_double) - 2.0_p_double * phi_np + real(this%pi%f3(1,i1,i2,i3-1),p_double) ) * idz2
          edotb = edotb_total_cell_center_3d( this, emf, i1, i2, i3 )
          this%diag_edotb%f3(1,i1,i2,i3) = real( edotb, p_k_fld )
          src = - this%g_agg * edotb
          this%edotb_max = max( this%edotb_max, abs(edotb) )
          this%edotb_sum = this%edotb_sum + edotb
          this%src_max = max( this%src_max, abs(src) )
          this%src_sum = this%src_sum + src
          this%src_count = this%src_count + 1
          phi_new = 2.0_p_double * phi_np - phi_nm + dt2 * ( lap - m2 * phi_np + src )
          this%a%f3(1,i1,i2,i3) = real( phi_new, p_k_fld )
        enddo
      enddo
    enddo
    do i3 = 1, this%nx(3)
      do i2 = 1, this%nx(2)
        do i1 = 1, this%nx(1)
          phi_np = real(this%pi%f3(1,i1,i2,i3),p_double)
          this%pi%f3(1,i1,i2,i3) = this%a%f3(1,i1,i2,i3)
          this%a%f3(1,i1,i2,i3) = real( phi_np, p_k_fld )
        enddo
      enddo
    enddo
  end select

  call apply_axion_bc_scalar( this, this%a )
  call apply_axion_bc_scalar( this, this%pi )
  call apply_axion_bc_scalar( this, this%diag_edotb )
end subroutine advance_driven_axion_field_impl


! Evaluate the centered linearised E dot B source in one spatial dimension.
!
! The coherent-field validation uses the field split
!
! E = E0 + dE, B = B0 + dB,
!
! but keeps only terms that are at most first order in the perturbative
! fields:
!
! E dot B -> E0 dot B0 + E0 dot dB + dE dot B0.
!
! The nonlinear perturbative product dE dot dB is deliberately omitted here.
! It can dominate once the perturbative response leaves the linear regime and
! obscures the Beyer-style early-time comparison.
real(p_double) function edotb_total_cell_center_1d( this, emf, i1 ) result(val)
  class(t_axion_field), intent(in) :: this
  class(t_emf), intent(in) :: emf
  integer, intent(in) :: i1
  real(p_double) :: ex0_c, ey0_c, ez0_c, bx0_c, by0_c, bz0_c
  real(p_double) :: exp_c, eyp_c, ezp_c, bxp_c, byp_c, bzp_c

  ex0_c = 0.5_p_double * ( real(emf%e%f1(1,i1),p_double) + real(emf%e%f1(1,i1-1),p_double) )
  ey0_c = real(emf%e%f1(2,i1),p_double)
  ez0_c = real(emf%e%f1(3,i1),p_double)
  bx0_c = real(emf%b%f1(1,i1),p_double)
  by0_c = 0.5_p_double * ( real(emf%b%f1(2,i1),p_double) + real(emf%b%f1(2,i1-1),p_double) )
  bz0_c = 0.5_p_double * ( real(emf%b%f1(3,i1),p_double) + real(emf%b%f1(3,i1-1),p_double) )

  exp_c = 0.0_p_double; eyp_c = 0.0_p_double; ezp_c = 0.0_p_double
  bxp_c = 0.0_p_double; byp_c = 0.0_p_double; bzp_c = 0.0_p_double
  if ( this%pert%if_initialized ) then
    exp_c = 0.5_p_double * ( real(this%pert%e%f1(1,i1),p_double) + real(this%pert%e%f1(1,i1-1),p_double) )
    eyp_c = real(this%pert%e%f1(2,i1),p_double)
    ezp_c = real(this%pert%e%f1(3,i1),p_double)
    bxp_c = real(this%pert%b%f1(1,i1),p_double)
    byp_c = 0.5_p_double * ( real(this%pert%b%f1(2,i1),p_double) + real(this%pert%b%f1(2,i1-1),p_double) )
    bzp_c = 0.5_p_double * ( real(this%pert%b%f1(3,i1),p_double) + real(this%pert%b%f1(3,i1-1),p_double) )
  endif

  val = ex0_c * bx0_c + ey0_c * by0_c + ez0_c * bz0_c + &
        ex0_c * bxp_c + ey0_c * byp_c + ez0_c * bzp_c + &
        exp_c * bx0_c + eyp_c * by0_c + ezp_c * bz0_c
end function edotb_total_cell_center_1d

! Evaluate the centered linearised E dot B source in two spatial dimensions.
real(p_double) function edotb_total_cell_center_2d( this, emf, i1, i2 ) result(val)
  class(t_axion_field), intent(in) :: this
  class(t_emf), intent(in) :: emf
  integer, intent(in) :: i1, i2
  real(p_double) :: ex0_c, ey0_c, ez0_c, bx0_c, by0_c, bz0_c
  real(p_double) :: exp_c, eyp_c, ezp_c, bxp_c, byp_c, bzp_c

  ex0_c = 0.5_p_double * ( real(emf%e%f2(1,i1,i2),p_double) + real(emf%e%f2(1,i1-1,i2),p_double) )
  ey0_c = 0.5_p_double * ( real(emf%e%f2(2,i1,i2),p_double) + real(emf%e%f2(2,i1,i2-1),p_double) )
  ez0_c = real(emf%e%f2(3,i1,i2),p_double)

  bx0_c = 0.5_p_double * ( real(emf%b%f2(1,i1,i2),p_double) + real(emf%b%f2(1,i1,i2-1),p_double) )
  by0_c = 0.5_p_double * ( real(emf%b%f2(2,i1,i2),p_double) + real(emf%b%f2(2,i1-1,i2),p_double) )
  bz0_c = 0.25_p_double * ( real(emf%b%f2(3,i1,i2),p_double) + real(emf%b%f2(3,i1-1,i2),p_double) + &
                             real(emf%b%f2(3,i1,i2-1),p_double) + real(emf%b%f2(3,i1-1,i2-1),p_double) )

  exp_c = 0.0_p_double; eyp_c = 0.0_p_double; ezp_c = 0.0_p_double
  bxp_c = 0.0_p_double; byp_c = 0.0_p_double; bzp_c = 0.0_p_double
  if ( this%pert%if_initialized ) then
    exp_c = 0.5_p_double * ( real(this%pert%e%f2(1,i1,i2),p_double) + real(this%pert%e%f2(1,i1-1,i2),p_double) )
    eyp_c = 0.5_p_double * ( real(this%pert%e%f2(2,i1,i2),p_double) + real(this%pert%e%f2(2,i1,i2-1),p_double) )
    ezp_c = real(this%pert%e%f2(3,i1,i2),p_double)
    bxp_c = 0.5_p_double * ( real(this%pert%b%f2(1,i1,i2),p_double) + real(this%pert%b%f2(1,i1,i2-1),p_double) )
    byp_c = 0.5_p_double * ( real(this%pert%b%f2(2,i1,i2),p_double) + real(this%pert%b%f2(2,i1-1,i2),p_double) )
    bzp_c = 0.25_p_double * ( real(this%pert%b%f2(3,i1,i2),p_double) + real(this%pert%b%f2(3,i1-1,i2),p_double) + &
                               real(this%pert%b%f2(3,i1,i2-1),p_double) + real(this%pert%b%f2(3,i1-1,i2-1),p_double) )
  endif

  val = ex0_c * bx0_c + ey0_c * by0_c + ez0_c * bz0_c + &
        ex0_c * bxp_c + ey0_c * byp_c + ez0_c * bzp_c + &
        exp_c * bx0_c + eyp_c * by0_c + ezp_c * bz0_c
end function edotb_total_cell_center_2d

! Evaluate the centered linearised E dot B source in three spatial dimensions.
real(p_double) function edotb_total_cell_center_3d( this, emf, i1, i2, i3 ) result(val)
  class(t_axion_field), intent(in) :: this
  class(t_emf), intent(in) :: emf
  integer, intent(in) :: i1, i2, i3
  real(p_double) :: ex0_c, ey0_c, ez0_c, bx0_c, by0_c, bz0_c
  real(p_double) :: exp_c, eyp_c, ezp_c, bxp_c, byp_c, bzp_c

  ex0_c = 0.5_p_double * ( real(emf%e%f3(1,i1,i2,i3),p_double) + real(emf%e%f3(1,i1-1,i2,i3),p_double) )
  ey0_c = 0.5_p_double * ( real(emf%e%f3(2,i1,i2,i3),p_double) + real(emf%e%f3(2,i1,i2-1,i3),p_double) )
  ez0_c = 0.5_p_double * ( real(emf%e%f3(3,i1,i2,i3),p_double) + real(emf%e%f3(3,i1,i2,i3-1),p_double) )

  bx0_c = 0.25_p_double * ( real(emf%b%f3(1,i1,i2,i3),p_double) + real(emf%b%f3(1,i1,i2-1,i3),p_double) + &
                             real(emf%b%f3(1,i1,i2,i3-1),p_double) + real(emf%b%f3(1,i1,i2-1,i3-1),p_double) )
  by0_c = 0.25_p_double * ( real(emf%b%f3(2,i1,i2,i3),p_double) + real(emf%b%f3(2,i1-1,i2,i3),p_double) + &
                             real(emf%b%f3(2,i1,i2,i3-1),p_double) + real(emf%b%f3(2,i1-1,i2,i3-1),p_double) )
  bz0_c = 0.25_p_double * ( real(emf%b%f3(3,i1,i2,i3),p_double) + real(emf%b%f3(3,i1-1,i2,i3),p_double) + &
                             real(emf%b%f3(3,i1,i2-1,i3),p_double) + real(emf%b%f3(3,i1-1,i2-1,i3),p_double) )

  exp_c = 0.0_p_double; eyp_c = 0.0_p_double; ezp_c = 0.0_p_double
  bxp_c = 0.0_p_double; byp_c = 0.0_p_double; bzp_c = 0.0_p_double
  if ( this%pert%if_initialized ) then
    exp_c = 0.5_p_double * ( real(this%pert%e%f3(1,i1,i2,i3),p_double) + real(this%pert%e%f3(1,i1-1,i2,i3),p_double) )
    eyp_c = 0.5_p_double * ( real(this%pert%e%f3(2,i1,i2,i3),p_double) + real(this%pert%e%f3(2,i1,i2-1,i3),p_double) )
    ezp_c = 0.5_p_double * ( real(this%pert%e%f3(3,i1,i2,i3),p_double) + real(this%pert%e%f3(3,i1,i2,i3-1),p_double) )
    bxp_c = 0.25_p_double * ( real(this%pert%b%f3(1,i1,i2,i3),p_double) + real(this%pert%b%f3(1,i1,i2-1,i3),p_double) + &
                               real(this%pert%b%f3(1,i1,i2,i3-1),p_double) + real(this%pert%b%f3(1,i1,i2-1,i3-1),p_double) )
    byp_c = 0.25_p_double * ( real(this%pert%b%f3(2,i1,i2,i3),p_double) + real(this%pert%b%f3(2,i1-1,i2,i3),p_double) + &
                               real(this%pert%b%f3(2,i1,i2,i3-1),p_double) + real(this%pert%b%f3(2,i1-1,i2,i3-1),p_double) )
    bzp_c = 0.25_p_double * ( real(this%pert%b%f3(3,i1,i2,i3),p_double) + real(this%pert%b%f3(3,i1-1,i2,i3),p_double) + &
                               real(this%pert%b%f3(3,i1,i2-1,i3),p_double) + real(this%pert%b%f3(3,i1-1,i2-1,i3),p_double) )
  endif

  val = ex0_c * bx0_c + ey0_c * by0_c + ez0_c * bz0_c + &
        ex0_c * bxp_c + ey0_c * byp_c + ez0_c * bzp_c + &
        exp_c * bx0_c + eyp_c * by0_c + ezp_c * bz0_c
end function edotb_total_cell_center_3d


! Accumulate the coherent axion backreaction into the perturbative EM response.
subroutine apply_backreaction_axion_field_impl( this, emf, dt_half, dt_full, jay1 )
  class(t_axion_field), intent(inout) :: this
  class(t_emf), intent(in) :: emf
  real(p_double), intent(in) :: dt_half, dt_full
  class(t_vdf), optional, intent(in) :: jay1

  integer :: i1, i2, i3
  real(p_double) :: curlx, curly, curlz
  real(p_double) :: sx, sy, sz, dmag, jx, jy, jz

  if ( .not. this%if_initialized ) return
  if ( .not. this%pert%if_initialized ) return

  this%de_max = 0.0_p_double
  this%de_sum = 0.0_p_double
  this%de_count = 0

  call apply_axion_bc_scalar( this, this%a )
  call apply_axion_bc_scalar( this, this%pi )
  call apply_axion_bc_vector( this, this%pert%b )

  select case ( p_x_dim )
  case (1)
    do i1 = 1, this%nx(1)
      curlx = 0.0_p_double
      curly = - ( real(this%pert%b%f1(3,i1+1),p_double) - real(this%pert%b%f1(3,i1),p_double) ) / this%dx(1)
      curlz = ( real(this%pert%b%f1(2,i1+1),p_double) - real(this%pert%b%f1(2,i1),p_double) ) / this%dx(1)

      sx = source_ex_1d( this, emf, dt_full, i1 )
      sy = source_ey_1d( this, emf, dt_full, i1 )
      sz = source_ez_1d( this, emf, dt_full, i1 )
      jx = 0.0_p_double; jy = 0.0_p_double; jz = 0.0_p_double
      if ( present(jay1) ) then
        jx = real( jay1%f1(1,i1), p_double )
        jy = real( jay1%f1(2,i1), p_double )
        jz = real( jay1%f1(3,i1), p_double )
      endif

      this%pert%e%f1(1,i1) = this%pert%e%f1(1,i1) + real( dt_half * ( curlx - sx - jx ), p_k_fld )
      this%pert%e%f1(2,i1) = this%pert%e%f1(2,i1) + real( dt_half * ( curly - sy - jy ), p_k_fld )
      this%pert%e%f1(3,i1) = this%pert%e%f1(3,i1) + real( dt_half * ( curlz - sz - jz ), p_k_fld )

      dmag = sqrt( sx*sx + sy*sy + sz*sz )
      this%diag_emsrcmag%f1(1,i1) = real( dmag, p_k_fld )
      this%de_max = max( this%de_max, abs(dt_half) * dmag )
      this%de_sum = this%de_sum + abs(dt_half) * dmag
      this%de_count = this%de_count + 1
    enddo
  case (2)
    do i2 = 1, this%nx(2)
      do i1 = 1, this%nx(1)
        curlx = ( real(this%pert%b%f2(3,i1,i2+1),p_double) - real(this%pert%b%f2(3,i1,i2),p_double) ) / this%dx(2)
        curly = - ( real(this%pert%b%f2(3,i1+1,i2),p_double) - real(this%pert%b%f2(3,i1,i2),p_double) ) / this%dx(1)
        curlz = ( real(this%pert%b%f2(2,i1+1,i2),p_double) - real(this%pert%b%f2(2,i1,i2),p_double) ) / this%dx(1) - &
                ( real(this%pert%b%f2(1,i1,i2+1),p_double) - real(this%pert%b%f2(1,i1,i2),p_double) ) / this%dx(2)

        sx = source_ex_2d( this, emf, dt_full, i1, i2 )
        sy = source_ey_2d( this, emf, dt_full, i1, i2 )
        sz = source_ez_2d( this, emf, dt_full, i1, i2 )
        jx = 0.0_p_double; jy = 0.0_p_double; jz = 0.0_p_double
        if ( present(jay1) ) then
          jx = real( jay1%f2(1,i1,i2), p_double )
          jy = real( jay1%f2(2,i1,i2), p_double )
          jz = real( jay1%f2(3,i1,i2), p_double )
        endif

        this%pert%e%f2(1,i1,i2) = this%pert%e%f2(1,i1,i2) + real( dt_half * ( curlx - sx - jx ), p_k_fld )
        this%pert%e%f2(2,i1,i2) = this%pert%e%f2(2,i1,i2) + real( dt_half * ( curly - sy - jy ), p_k_fld )
        this%pert%e%f2(3,i1,i2) = this%pert%e%f2(3,i1,i2) + real( dt_half * ( curlz - sz - jz ), p_k_fld )

        dmag = sqrt( sx*sx + sy*sy + sz*sz )
        this%diag_emsrcmag%f2(1,i1,i2) = real( dmag, p_k_fld )
        this%de_max = max( this%de_max, abs(dt_half) * dmag )
        this%de_sum = this%de_sum + abs(dt_half) * dmag
        this%de_count = this%de_count + 1
      enddo
    enddo
  case (3)
    do i3 = 1, this%nx(3)
      do i2 = 1, this%nx(2)
        do i1 = 1, this%nx(1)
          curlx = ( real(this%pert%b%f3(3,i1,i2+1,i3),p_double) - real(this%pert%b%f3(3,i1,i2,i3),p_double) ) / this%dx(2) - &
                  ( real(this%pert%b%f3(2,i1,i2,i3+1),p_double) - real(this%pert%b%f3(2,i1,i2,i3),p_double) ) / this%dx(3)
          curly = ( real(this%pert%b%f3(1,i1,i2,i3+1),p_double) - real(this%pert%b%f3(1,i1,i2,i3),p_double) ) / this%dx(3) - &
                  ( real(this%pert%b%f3(3,i1+1,i2,i3),p_double) - real(this%pert%b%f3(3,i1,i2,i3),p_double) ) / this%dx(1)
          curlz = ( real(this%pert%b%f3(2,i1+1,i2,i3),p_double) - real(this%pert%b%f3(2,i1,i2,i3),p_double) ) / this%dx(1) - &
                  ( real(this%pert%b%f3(1,i1,i2+1,i3),p_double) - real(this%pert%b%f3(1,i1,i2,i3),p_double) ) / this%dx(2)

          sx = source_ex_3d( this, emf, dt_full, i1, i2, i3 )
          sy = source_ey_3d( this, emf, dt_full, i1, i2, i3 )
          sz = source_ez_3d( this, emf, dt_full, i1, i2, i3 )
          jx = 0.0_p_double; jy = 0.0_p_double; jz = 0.0_p_double
          if ( present(jay1) ) then
            jx = real( jay1%f3(1,i1,i2,i3), p_double )
            jy = real( jay1%f3(2,i1,i2,i3), p_double )
            jz = real( jay1%f3(3,i1,i2,i3), p_double )
          endif

          this%pert%e%f3(1,i1,i2,i3) = this%pert%e%f3(1,i1,i2,i3) + real( dt_half * ( curlx - sx - jx ), p_k_fld )
          this%pert%e%f3(2,i1,i2,i3) = this%pert%e%f3(2,i1,i2,i3) + real( dt_half * ( curly - sy - jy ), p_k_fld )
          this%pert%e%f3(3,i1,i2,i3) = this%pert%e%f3(3,i1,i2,i3) + real( dt_half * ( curlz - sz - jz ), p_k_fld )

          dmag = sqrt( sx*sx + sy*sy + sz*sz )
          this%diag_emsrcmag%f3(1,i1,i2,i3) = real( dmag, p_k_fld )
          this%de_max = max( this%de_max, abs(dt_half) * dmag )
          this%de_sum = this%de_sum + abs(dt_half) * dmag
          this%de_count = this%de_count + 1
        enddo
      enddo
    enddo
  end select

  call apply_axion_bc_vector( this, this%pert%e )
  call apply_axion_bc_scalar( this, this%diag_emsrcmag )
end subroutine apply_backreaction_axion_field_impl




! Zero the perturbative EM response fields.
subroutine reset_perturbation_axion_field_impl( this )
  class(t_axion_field), intent(inout) :: this

  integer :: i1, i2, i3, ic
  integer :: l1, u1, l2, u2, l3, u3

  if ( .not. this%if_initialized ) return
  if ( .not. this%pert%if_initialized ) return

  this%edotb_max = 0.0_p_double
  this%edotb_sum = 0.0_p_double
  this%src_max = 0.0_p_double
  this%src_sum = 0.0_p_double
  this%src_count = 0
  this%de_max = 0.0_p_double
  this%de_sum = 0.0_p_double
  this%de_count = 0

  select case ( p_x_dim )
  case (1)
    l1 = lbound(this%pert%e%f1,2); u1 = ubound(this%pert%e%f1,2)
    do ic = 1, 3
      do i1 = l1, u1
        this%pert%e%f1(ic,i1) = 0.0_p_k_fld
        this%pert%b%f1(ic,i1) = 0.0_p_k_fld
      enddo
    enddo
  case (2)
    l1 = lbound(this%pert%e%f2,2); u1 = ubound(this%pert%e%f2,2)
    l2 = lbound(this%pert%e%f2,3); u2 = ubound(this%pert%e%f2,3)
    do ic = 1, 3
      do i2 = l2, u2
        do i1 = l1, u1
          this%pert%e%f2(ic,i1,i2) = 0.0_p_k_fld
          this%pert%b%f2(ic,i1,i2) = 0.0_p_k_fld
        enddo
      enddo
    enddo
  case (3)
    l1 = lbound(this%pert%e%f3,2); u1 = ubound(this%pert%e%f3,2)
    l2 = lbound(this%pert%e%f3,3); u2 = ubound(this%pert%e%f3,3)
    l3 = lbound(this%pert%e%f3,4); u3 = ubound(this%pert%e%f3,4)
    do ic = 1, 3
      do i3 = l3, u3
        do i2 = l2, u2
          do i1 = l1, u1
            this%pert%e%f3(ic,i1,i2,i3) = 0.0_p_k_fld
            this%pert%b%f3(ic,i1,i2,i3) = 0.0_p_k_fld
          enddo
        enddo
      enddo
    enddo
  end select

  select case ( p_x_dim )
  case (1)
    l1 = lbound(this%diag_edotb%f1,2); u1 = ubound(this%diag_edotb%f1,2)
    do i1 = l1, u1
      this%diag_edotb%f1(1,i1) = 0.0_p_k_fld
      this%diag_emsrcmag%f1(1,i1) = 0.0_p_k_fld
    enddo
  case (2)
    l1 = lbound(this%diag_edotb%f2,2); u1 = ubound(this%diag_edotb%f2,2)
    l2 = lbound(this%diag_edotb%f2,3); u2 = ubound(this%diag_edotb%f2,3)
    do i2 = l2, u2
      do i1 = l1, u1
        this%diag_edotb%f2(1,i1,i2) = 0.0_p_k_fld
        this%diag_emsrcmag%f2(1,i1,i2) = 0.0_p_k_fld
      enddo
    enddo
  case (3)
    l1 = lbound(this%diag_edotb%f3,2); u1 = ubound(this%diag_edotb%f3,2)
    l2 = lbound(this%diag_edotb%f3,3); u2 = ubound(this%diag_edotb%f3,3)
    l3 = lbound(this%diag_edotb%f3,4); u3 = ubound(this%diag_edotb%f3,4)
    do i3 = l3, u3
      do i2 = l2, u2
        do i1 = l1, u1
          this%diag_edotb%f3(1,i1,i2,i3) = 0.0_p_k_fld
          this%diag_emsrcmag%f3(1,i1,i2,i3) = 0.0_p_k_fld
        enddo
      enddo
    enddo
  end select

  call apply_axion_bc_vector( this, this%pert%e )
  call apply_axion_bc_vector( this, this%pert%b )
  call apply_axion_bc_scalar( this, this%diag_edotb )
  call apply_axion_bc_scalar( this, this%diag_emsrcmag )
end subroutine reset_perturbation_axion_field_impl



! Advance the perturbative magnetic field by a half step.
subroutine advance_perturb_b_axion_field_impl( this, dt_half )
  class(t_axion_field), intent(inout) :: this
  real(p_double), intent(in) :: dt_half

  integer :: i1, i2, i3
  real(p_double) :: curlx, curly, curlz

  if ( .not. this%if_initialized ) return
  if ( .not. this%pert%if_initialized ) return

  call apply_axion_bc_vector( this, this%pert%e )

  select case ( p_x_dim )
  case (1)
    do i1 = 1, this%nx(1)
      ! In 1D, Faraday's law gives
      ! dB2/dt = + dE3/dx, dB3/dt = - dE2/dx.
      ! Keep the y/z curl components paired with the corresponding B update.
      curly = ( real(this%pert%e%f1(3,i1),p_double) - real(this%pert%e%f1(3,i1-1),p_double) ) / this%dx(1)
      curlz = ( real(this%pert%e%f1(2,i1),p_double) - real(this%pert%e%f1(2,i1-1),p_double) ) / this%dx(1)
      this%pert%b%f1(2,i1) = this%pert%b%f1(2,i1) + real( dt_half * curly, p_k_fld )
      this%pert%b%f1(3,i1) = this%pert%b%f1(3,i1) - real( dt_half * curlz, p_k_fld )
    enddo
  case (2)
    do i2 = 1, this%nx(2)
      do i1 = 1, this%nx(1)
        curlx = ( real(this%pert%e%f2(3,i1,i2),p_double) - real(this%pert%e%f2(3,i1,i2-1),p_double) ) / this%dx(2)
        curly = ( real(this%pert%e%f2(3,i1,i2),p_double) - real(this%pert%e%f2(3,i1-1,i2),p_double) ) / this%dx(1)
        curlz = ( real(this%pert%e%f2(2,i1,i2),p_double) - real(this%pert%e%f2(2,i1-1,i2),p_double) ) / this%dx(1) - &
                ( real(this%pert%e%f2(1,i1,i2),p_double) - real(this%pert%e%f2(1,i1,i2-1),p_double) ) / this%dx(2)
        this%pert%b%f2(1,i1,i2) = this%pert%b%f2(1,i1,i2) - real( dt_half * curlx, p_k_fld )
        this%pert%b%f2(2,i1,i2) = this%pert%b%f2(2,i1,i2) + real( dt_half * curly, p_k_fld )
        this%pert%b%f2(3,i1,i2) = this%pert%b%f2(3,i1,i2) - real( dt_half * curlz, p_k_fld )
      enddo
    enddo
  case (3)
    do i3 = 1, this%nx(3)
      do i2 = 1, this%nx(2)
        do i1 = 1, this%nx(1)
          curlx = ( real(this%pert%e%f3(3,i1,i2,i3),p_double) - real(this%pert%e%f3(3,i1,i2-1,i3),p_double) ) / this%dx(2) - &
                  ( real(this%pert%e%f3(2,i1,i2,i3),p_double) - real(this%pert%e%f3(2,i1,i2,i3-1),p_double) ) / this%dx(3)
          curly = ( real(this%pert%e%f3(3,i1,i2,i3),p_double) - real(this%pert%e%f3(3,i1-1,i2,i3),p_double) ) / this%dx(1) - &
                  ( real(this%pert%e%f3(1,i1,i2,i3),p_double) - real(this%pert%e%f3(1,i1,i2,i3-1),p_double) ) / this%dx(3)
          curlz = ( real(this%pert%e%f3(2,i1,i2,i3),p_double) - real(this%pert%e%f3(2,i1-1,i2,i3),p_double) ) / this%dx(1) - &
                  ( real(this%pert%e%f3(1,i1,i2,i3),p_double) - real(this%pert%e%f3(1,i1,i2-1,i3),p_double) ) / this%dx(2)
          this%pert%b%f3(1,i1,i2,i3) = this%pert%b%f3(1,i1,i2,i3) - real( dt_half * curlx, p_k_fld )
          this%pert%b%f3(2,i1,i2,i3) = this%pert%b%f3(2,i1,i2,i3) + real( dt_half * curly, p_k_fld )
          this%pert%b%f3(3,i1,i2,i3) = this%pert%b%f3(3,i1,i2,i3) - real( dt_half * curlz, p_k_fld )
        enddo
      enddo
    enddo
  end select

  call apply_axion_bc_vector( this, this%pert%b )
end subroutine advance_perturb_b_axion_field_impl



! Return the centered axion field value in one spatial dimension.
real(p_double) function phi_center_1d( this, i1 ) result(val)
  class(t_axion_field), intent(in) :: this
  integer, intent(in) :: i1
  val = 0.5_p_double * ( real(this%a%f1(1,i1),p_double) + real(this%pi%f1(1,i1),p_double) )
end function phi_center_1d

! Return the centered axion field value in two spatial dimensions.
real(p_double) function phi_center_2d( this, i1, i2 ) result(val)
  class(t_axion_field), intent(in) :: this
  integer, intent(in) :: i1, i2
  val = 0.5_p_double * ( real(this%a%f2(1,i1,i2),p_double) + real(this%pi%f2(1,i1,i2),p_double) )
end function phi_center_2d

! Return the centered axion field value in three spatial dimensions.
real(p_double) function phi_center_3d( this, i1, i2, i3 ) result(val)
  class(t_axion_field), intent(in) :: this
  integer, intent(in) :: i1, i2, i3
  val = 0.5_p_double * ( real(this%a%f3(1,i1,i2,i3),p_double) + real(this%pi%f3(1,i1,i2,i3),p_double) )
end function phi_center_3d

! Return the centered axion time derivative in one spatial dimension.
real(p_double) function phi_t_1d( this, dt_full, i1 ) result(val)
  class(t_axion_field), intent(in) :: this
  real(p_double), intent(in) :: dt_full
  integer, intent(in) :: i1
  val = ( real(this%pi%f1(1,i1),p_double) - real(this%a%f1(1,i1),p_double) ) / max( dt_full, tiny(1.0_p_double) )
end function phi_t_1d

! Return the centered axion time derivative in two spatial dimensions.
real(p_double) function phi_t_2d( this, dt_full, i1, i2 ) result(val)
  class(t_axion_field), intent(in) :: this
  real(p_double), intent(in) :: dt_full
  integer, intent(in) :: i1, i2
  val = ( real(this%pi%f2(1,i1,i2),p_double) - real(this%a%f2(1,i1,i2),p_double) ) / max( dt_full, tiny(1.0_p_double) )
end function phi_t_2d

! Return the centered axion time derivative in three spatial dimensions.
real(p_double) function phi_t_3d( this, dt_full, i1, i2, i3 ) result(val)
  class(t_axion_field), intent(in) :: this
  real(p_double), intent(in) :: dt_full
  integer, intent(in) :: i1, i2, i3
  val = ( real(this%pi%f3(1,i1,i2,i3),p_double) - real(this%a%f3(1,i1,i2,i3),p_double) ) / max( dt_full, tiny(1.0_p_double) )
end function phi_t_3d

! Evaluate the driven perturbative Ex source in one spatial dimension.
real(p_double) function source_ex_1d( this, emf, dt_full, i1 ) result(val)
  class(t_axion_field), intent(in) :: this
  class(t_emf), intent(in) :: emf
  real(p_double), intent(in) :: dt_full
  integer, intent(in) :: i1
  real(p_double) :: pt, bx
  pt = 0.5_p_double * ( phi_t_1d(this,dt_full,i1) + phi_t_1d(this,dt_full,i1+1) )
  bx = real(emf%b%f1(1,i1),p_double)
  val = this%g_agg * pt * bx
end function source_ex_1d

! Evaluate the driven perturbative Ey source in one spatial dimension.
real(p_double) function source_ey_1d( this, emf, dt_full, i1 ) result(val)
  class(t_axion_field), intent(in) :: this
  class(t_emf), intent(in) :: emf
  real(p_double), intent(in) :: dt_full
  integer, intent(in) :: i1
  real(p_double) :: pt, by, gradx, ez
  pt = phi_t_1d(this,dt_full,i1)
  by = 0.5_p_double * ( real(emf%b%f1(2,i1),p_double) + real(emf%b%f1(2,i1-1),p_double) )
  gradx = ( phi_center_1d(this,i1+1) - phi_center_1d(this,i1-1) ) / ( 2.0_p_double * this%dx(1) )
  ez = real(emf%e%f1(3,i1),p_double)
  val = this%g_agg * ( pt * by - ez * gradx )
end function source_ey_1d

! Evaluate the driven perturbative Ez source in one spatial dimension.
real(p_double) function source_ez_1d( this, emf, dt_full, i1 ) result(val)
  class(t_axion_field), intent(in) :: this
  class(t_emf), intent(in) :: emf
  real(p_double), intent(in) :: dt_full
  integer, intent(in) :: i1
  real(p_double) :: pt, bz, gradx, ey
  pt = phi_t_1d(this,dt_full,i1)
  bz = 0.5_p_double * ( real(emf%b%f1(3,i1),p_double) + real(emf%b%f1(3,i1-1),p_double) )
  gradx = ( phi_center_1d(this,i1+1) - phi_center_1d(this,i1-1) ) / ( 2.0_p_double * this%dx(1) )
  ey = real(emf%e%f1(2,i1),p_double)
  val = this%g_agg * ( pt * bz + ey * gradx )
end function source_ez_1d

! Evaluate the driven perturbative Ex source in two spatial dimensions.
real(p_double) function source_ex_2d( this, emf, dt_full, i1, i2 ) result(val)
  class(t_axion_field), intent(in) :: this
  class(t_emf), intent(in) :: emf
  real(p_double), intent(in) :: dt_full
  integer, intent(in) :: i1, i2
  real(p_double) :: pt, bx, ey, dzphi, ez, dyphi
  pt = 0.5_p_double * ( phi_t_2d(this,dt_full,i1,i2) + phi_t_2d(this,dt_full,i1+1,i2) )
  bx = 0.5_p_double * ( real(emf%b%f2(1,i1,i2),p_double) + real(emf%b%f2(1,i1,i2-1),p_double) )
  ey = 0.5_p_double * ( real(emf%e%f2(2,i1,i2),p_double) + real(emf%e%f2(2,i1+1,i2),p_double) )
  dzphi = 0.0_p_double
  ez = 0.5_p_double * ( real(emf%e%f2(3,i1,i2),p_double) + real(emf%e%f2(3,i1+1,i2),p_double) )
  dyphi = 0.5_p_double * ( ( phi_center_2d(this,i1,i2+1) - phi_center_2d(this,i1,i2-1) ) + &
                            ( phi_center_2d(this,i1+1,i2+1) - phi_center_2d(this,i1+1,i2-1) ) ) / ( 2.0_p_double * this%dx(2) )
  val = this%g_agg * ( pt * bx - ey * dzphi + ez * dyphi )
end function source_ex_2d

! Evaluate the driven perturbative Ey source in two spatial dimensions.
real(p_double) function source_ey_2d( this, emf, dt_full, i1, i2 ) result(val)
  class(t_axion_field), intent(in) :: this
  class(t_emf), intent(in) :: emf
  real(p_double), intent(in) :: dt_full
  integer, intent(in) :: i1, i2
  real(p_double) :: pt, by, ez, dxphi
  pt = 0.5_p_double * ( phi_t_2d(this,dt_full,i1,i2) + phi_t_2d(this,dt_full,i1,i2+1) )
  by = 0.5_p_double * ( real(emf%b%f2(2,i1,i2),p_double) + real(emf%b%f2(2,i1-1,i2),p_double) )
  ez = 0.5_p_double * ( real(emf%e%f2(3,i1,i2),p_double) + real(emf%e%f2(3,i1,i2+1),p_double) )
  dxphi = 0.5_p_double * ( ( phi_center_2d(this,i1+1,i2) - phi_center_2d(this,i1-1,i2) ) + &
                            ( phi_center_2d(this,i1+1,i2+1) - phi_center_2d(this,i1-1,i2+1) ) ) / ( 2.0_p_double * this%dx(1) )
  val = this%g_agg * ( pt * by - ez * dxphi )
end function source_ey_2d

! Evaluate the driven perturbative Ez source in two spatial dimensions.
real(p_double) function source_ez_2d( this, emf, dt_full, i1, i2 ) result(val)
  class(t_axion_field), intent(in) :: this
  class(t_emf), intent(in) :: emf
  real(p_double), intent(in) :: dt_full
  integer, intent(in) :: i1, i2
  real(p_double) :: pt, bz, ex, dyphi, ey, dxphi
  pt = phi_t_2d(this,dt_full,i1,i2)
  bz = 0.25_p_double * ( real(emf%b%f2(3,i1,i2),p_double) + real(emf%b%f2(3,i1-1,i2),p_double) + &
                         real(emf%b%f2(3,i1,i2-1),p_double) + real(emf%b%f2(3,i1-1,i2-1),p_double) )
  ex = 0.5_p_double * ( real(emf%e%f2(1,i1,i2),p_double) + real(emf%e%f2(1,i1,i2-1),p_double) )
  ey = 0.5_p_double * ( real(emf%e%f2(2,i1,i2),p_double) + real(emf%e%f2(2,i1-1,i2),p_double) )
  dyphi = ( phi_center_2d(this,i1,i2+1) - phi_center_2d(this,i1,i2-1) ) / ( 2.0_p_double * this%dx(2) )
  dxphi = ( phi_center_2d(this,i1+1,i2) - phi_center_2d(this,i1-1,i2) ) / ( 2.0_p_double * this%dx(1) )
  val = this%g_agg * ( pt * bz - ex * dyphi + ey * dxphi )
end function source_ez_2d

! Evaluate the driven perturbative Ex source in three spatial dimensions.
real(p_double) function source_ex_3d( this, emf, dt_full, i1, i2, i3 ) result(val)
  class(t_axion_field), intent(in) :: this
  class(t_emf), intent(in) :: emf
  real(p_double), intent(in) :: dt_full
  integer, intent(in) :: i1, i2, i3
  real(p_double) :: pt, bx, term1, term2
  pt = 0.5_p_double * ( phi_t_3d(this,dt_full,i1,i2,i3) + phi_t_3d(this,dt_full,i1+1,i2,i3) )
  bx = 0.25_p_double * ( real(emf%b%f3(1,i1,i2,i3),p_double) + real(emf%b%f3(1,i1,i2-1,i3),p_double) + &
                         real(emf%b%f3(1,i1,i2,i3-1),p_double) + real(emf%b%f3(1,i1,i2-1,i3-1),p_double) )
  term1 = 0.25_p_double * ( real(emf%e%f3(2,i1,i2,i3),p_double) + real(emf%e%f3(2,i1,i2-1,i3),p_double) ) * &
          ( phi_center_3d(this,i1,i2,i3+1) - phi_center_3d(this,i1,i2,i3-1) ) / ( 2.0_p_double * this%dx(3) ) + &
          0.25_p_double * ( real(emf%e%f3(2,i1+1,i2,i3),p_double) + real(emf%e%f3(2,i1+1,i2-1,i3),p_double) ) * &
          ( phi_center_3d(this,i1+1,i2,i3+1) - phi_center_3d(this,i1+1,i2,i3-1) ) / ( 2.0_p_double * this%dx(3) )
  term2 = 0.25_p_double * ( real(emf%e%f3(3,i1,i2,i3),p_double) + real(emf%e%f3(3,i1,i2,i3-1),p_double) ) * &
          ( phi_center_3d(this,i1,i2+1,i3) - phi_center_3d(this,i1,i2-1,i3) ) / ( 2.0_p_double * this%dx(2) ) + &
          0.25_p_double * ( real(emf%e%f3(3,i1+1,i2,i3),p_double) + real(emf%e%f3(3,i1+1,i2,i3-1),p_double) ) * &
          ( phi_center_3d(this,i1+1,i2+1,i3) - phi_center_3d(this,i1+1,i2-1,i3) ) / ( 2.0_p_double * this%dx(2) )
  val = this%g_agg * ( pt * bx - term1 + term2 )
end function source_ex_3d

! Evaluate the driven perturbative Ey source in three spatial dimensions.
real(p_double) function source_ey_3d( this, emf, dt_full, i1, i2, i3 ) result(val)
  class(t_axion_field), intent(in) :: this
  class(t_emf), intent(in) :: emf
  real(p_double), intent(in) :: dt_full
  integer, intent(in) :: i1, i2, i3
  real(p_double) :: pt, by, term1, term2
  pt = 0.5_p_double * ( phi_t_3d(this,dt_full,i1,i2,i3) + phi_t_3d(this,dt_full,i1,i2+1,i3) )
  by = 0.25_p_double * ( real(emf%b%f3(2,i1,i2,i3),p_double) + real(emf%b%f3(2,i1-1,i2,i3),p_double) + &
                         real(emf%b%f3(2,i1,i2,i3-1),p_double) + real(emf%b%f3(2,i1-1,i2,i3-1),p_double) )
  term1 = 0.25_p_double * ( real(emf%e%f3(3,i1,i2,i3),p_double) + real(emf%e%f3(3,i1-1,i2,i3),p_double) ) * &
          ( phi_center_3d(this,i1+1,i2,i3) - phi_center_3d(this,i1-1,i2,i3) ) / ( 2.0_p_double * this%dx(1) ) + &
          0.25_p_double * ( real(emf%e%f3(3,i1,i2+1,i3),p_double) + real(emf%e%f3(3,i1-1,i2+1,i3),p_double) ) * &
          ( phi_center_3d(this,i1+1,i2+1,i3) - phi_center_3d(this,i1-1,i2+1,i3) ) / ( 2.0_p_double * this%dx(1) )
  term2 = 0.25_p_double * ( real(emf%e%f3(1,i1,i2,i3),p_double) + real(emf%e%f3(1,i1,i2,i3-1),p_double) ) * &
          ( phi_center_3d(this,i1,i2,i3+1) - phi_center_3d(this,i1,i2,i3-1) ) / ( 2.0_p_double * this%dx(3) ) + &
          0.25_p_double * ( real(emf%e%f3(1,i1,i2+1,i3),p_double) + real(emf%e%f3(1,i1,i2+1,i3-1),p_double) ) * &
          ( phi_center_3d(this,i1,i2+1,i3+1) - phi_center_3d(this,i1,i2+1,i3-1) ) / ( 2.0_p_double * this%dx(3) )
  val = this%g_agg * ( pt * by - term1 + term2 )
end function source_ey_3d

! Evaluate the driven perturbative Ez source in three spatial dimensions.
real(p_double) function source_ez_3d( this, emf, dt_full, i1, i2, i3 ) result(val)
  class(t_axion_field), intent(in) :: this
  class(t_emf), intent(in) :: emf
  real(p_double), intent(in) :: dt_full
  integer, intent(in) :: i1, i2, i3
  real(p_double) :: pt, bz, term1, term2
  pt = 0.5_p_double * ( phi_t_3d(this,dt_full,i1,i2,i3) + phi_t_3d(this,dt_full,i1,i2,i3+1) )
  bz = 0.25_p_double * ( real(emf%b%f3(3,i1,i2,i3),p_double) + real(emf%b%f3(3,i1-1,i2,i3),p_double) + &
                         real(emf%b%f3(3,i1,i2-1,i3),p_double) + real(emf%b%f3(3,i1-1,i2-1,i3),p_double) )
  term1 = 0.25_p_double * ( real(emf%e%f3(1,i1,i2,i3),p_double) + real(emf%e%f3(1,i1,i2-1,i3),p_double) ) * &
          ( phi_center_3d(this,i1,i2+1,i3) - phi_center_3d(this,i1,i2-1,i3) ) / ( 2.0_p_double * this%dx(2) ) + &
          0.25_p_double * ( real(emf%e%f3(1,i1,i2,i3+1),p_double) + real(emf%e%f3(1,i1,i2-1,i3+1),p_double) ) * &
          ( phi_center_3d(this,i1,i2+1,i3+1) - phi_center_3d(this,i1,i2-1,i3+1) ) / ( 2.0_p_double * this%dx(2) )
  term2 = 0.25_p_double * ( real(emf%e%f3(2,i1,i2,i3),p_double) + real(emf%e%f3(2,i1-1,i2,i3),p_double) ) * &
          ( phi_center_3d(this,i1+1,i2,i3) - phi_center_3d(this,i1-1,i2,i3) ) / ( 2.0_p_double * this%dx(1) ) + &
          0.25_p_double * ( real(emf%e%f3(2,i1,i2,i3+1),p_double) + real(emf%e%f3(2,i1-1,i2,i3+1),p_double) ) * &
          ( phi_center_3d(this,i1+1,i2,i3+1) - phi_center_3d(this,i1-1,i2,i3+1) ) / ( 2.0_p_double * this%dx(1) )
  val = this%g_agg * ( pt * bz - term1 + term2 )
end function source_ez_3d


end module m_axion_field_solver_impl

! Forward the free axion-field update to the compiled implementation module.
subroutine advance_free_axion_field( this, dt, no_co, send_msg, recv_msg )
  use m_parameters, only : p_double
  use m_axion_field_define, only : t_axion_field
  use m_node_conf, only : t_node_conf
  use m_vdf_comm, only : t_vdf_msg
  use m_axion_field_solver_impl, only : advance_free_axion_field_impl
  implicit none
  class(t_axion_field), intent(inout) :: this
  real(p_double), intent(in) :: dt
  class(t_node_conf), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg
  call advance_free_axion_field_impl( this, dt, no_co, send_msg, recv_msg )
end subroutine advance_free_axion_field

! Forward the driven axion-field update to the compiled implementation module.
subroutine advance_driven_axion_field( this, emf, dt, no_co, send_msg, recv_msg )
  use m_parameters, only : p_double
  use m_axion_field_define, only : t_axion_field
  use m_emf_define, only : t_emf
  use m_node_conf, only : t_node_conf
  use m_vdf_comm, only : t_vdf_msg
  use m_axion_field_solver_impl, only : advance_driven_axion_field_impl
  implicit none
  class(t_axion_field), intent(inout) :: this
  class(t_emf), intent(in) :: emf
  real(p_double), intent(in) :: dt
  class(t_node_conf), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg
  call advance_driven_axion_field_impl( this, emf, dt, no_co, send_msg, recv_msg )
end subroutine advance_driven_axion_field

! Forward the axion backreaction update to the compiled implementation module.
subroutine apply_backreaction_axion_field( this, emf, dt_half, dt_full, jay1 )
  use m_parameters, only : p_double
  use m_axion_field_define, only : t_axion_field
  use m_emf_define, only : t_emf
  use m_vdf_define, only : t_vdf
  use m_axion_field_solver_impl, only : apply_backreaction_axion_field_impl
  implicit none
  class(t_axion_field), intent(inout) :: this
  class(t_emf), intent(in) :: emf
  real(p_double), intent(in) :: dt_half, dt_full
  class(t_vdf), optional, intent(in) :: jay1

  if ( present(jay1) ) then
    call apply_backreaction_axion_field_impl( this, emf, dt_half, dt_full, jay1 )
  else
    call apply_backreaction_axion_field_impl( this, emf, dt_half, dt_full )
  endif
end subroutine apply_backreaction_axion_field

! Forward the perturbation reset to the compiled implementation module.
subroutine reset_perturbation_axion_field( this )
  use m_axion_field_define, only : t_axion_field
  use m_axion_field_solver_impl, only : reset_perturbation_axion_field_impl
  implicit none
  class(t_axion_field), intent(inout) :: this
  call reset_perturbation_axion_field_impl( this )
end subroutine reset_perturbation_axion_field

! Forward the perturbative magnetic-field half-step to the compiled implementation module.
subroutine advance_perturb_b_axion_field( this, dt_half )
  use m_parameters, only : p_double
  use m_axion_field_define, only : t_axion_field
  use m_axion_field_solver_impl, only : advance_perturb_b_axion_field_impl
  implicit none
  class(t_axion_field), intent(inout) :: this
  real(p_double), intent(in) :: dt_half
  call advance_perturb_b_axion_field_impl( this, dt_half )
end subroutine advance_perturb_b_axion_field
