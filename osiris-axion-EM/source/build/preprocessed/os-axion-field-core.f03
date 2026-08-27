# 1 "axion_field/os-axion-field-core.f03"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "axion_field/os-axion-field-core.f03"
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
# 2 "axion_field/os-axion-field-core.f03" 2
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
# 3 "axion_field/os-axion-field-core.f03" 2

module m_axion_field_core_impl

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
  use m_axion_field_boundary, only : canonical_axion_bc, apply_axion_bc_scalar, apply_axion_bc_vector

  implicit none
  private

  public :: read_input_axion_field_impl
  public :: init_axion_field_impl
  public :: cleanup_axion_field_impl
  public :: move_window_axion_field_impl
  public :: update_boundary_axion_field_impl
  public :: list_algorithm_axion_field_impl

contains

! Read the axion-field namelist and store the user configuration.
subroutine read_input_axion_field_impl( this, input_file )
  class(t_axion_field), intent(inout) :: this
  class(t_input_file), intent(inout) :: input_file

  real(p_double) :: ma, g_agg, init_ampl, init_lambda, init_phase, init_x0, init_y0, init_z0
  real(p_double) :: init_dir(3), init_w_long, init_w_trans1, init_w_trans2, init_vphase
  real(p_double) :: init_phase_step
  character(len=32) :: init_type
  character(len=32) :: bc_x_lower, bc_x_upper, bc_y_lower, bc_y_upper, bc_z_lower, bc_z_upper
  integer :: init_kx, init_ky, init_kz, init_kx_min, init_kx_max
  logical :: init_band_norm
  integer :: ndump_fac, ndump_fac_ave, ndump_fac_lineout, prec, n_tavg
  integer, dimension(p_x_dim) :: n_ave
  integer, dimension(p_n_report_type) :: ndump_fac_all
  character( len = p_max_reports_len ), dimension( p_max_reports ) :: reports
  integer :: ierr

  namelist /nl_axion_field/ ma, g_agg, init_type, init_ampl, init_kx, init_ky, init_kz, &
                            init_kx_min, init_kx_max, init_phase_step, init_band_norm, &
                            init_lambda, init_phase, init_x0, init_y0, init_z0, init_dir, &
                            init_w_long, init_w_trans1, init_w_trans2, init_vphase, &
                            bc_x_lower, bc_x_upper, bc_y_lower, bc_y_upper, bc_z_lower, bc_z_upper, &
                            ndump_fac, ndump_fac_ave, ndump_fac_lineout, prec, n_ave, n_tavg, reports

  ma = 0.0_p_double
  g_agg = 0.0_p_double
  init_type = 'zero'
  init_ampl = 0.0_p_double
  init_kx = 1
  init_ky = 0
  init_kz = 0
  init_kx_min = 1
  init_kx_max = 1
  init_phase_step = 0.0_p_double
  init_band_norm = .false.
  init_lambda = 0.0_p_double
  init_phase = 0.0_p_double
  init_x0 = 0.0_p_double
  init_y0 = 0.0_p_double
  init_z0 = 0.0_p_double
  init_dir = (/ 1.0_p_double, 0.0_p_double, 0.0_p_double /)
  init_w_long = 0.0_p_double
  init_w_trans1 = 0.0_p_double
  init_w_trans2 = 0.0_p_double
  init_vphase = 1.0_p_double
  bc_x_lower = 'periodic'
  bc_x_upper = 'periodic'
  bc_y_lower = 'periodic'
  bc_y_upper = 'periodic'
  bc_z_lower = 'periodic'
  bc_z_upper = 'periodic'
  ndump_fac = 0
  ndump_fac_ave = 0
  ndump_fac_lineout = 0
  prec = p_single
  n_ave = -1
  n_tavg = -1
  reports = '-'

  call get_namelist( input_file, 'nl_axion_field', ierr )

  if ( ierr == 0 ) then
    read( input_file%nml_text, nml = nl_axion_field, iostat = ierr )
    if ( ierr /= 0 ) then
      write(err_buf__,*) 'Error reading axion_field parameters';call err__("axion_field/os-axion-field-core.f03",101)
      call abort_program( p_err_invalid )
    endif
  else if ( ierr < 0 ) then
    write(err_buf__,*) 'Error processing axion_field parameters';call err__("axion_field/os-axion-field-core.f03",105)
    call abort_program( p_err_invalid )
  endif

  this%ma = ma
  this%g_agg = g_agg
  this%init_type = trim(init_type)
  this%init_ampl = init_ampl
  if ( ( trim(init_type) == 'sin_band' .or. trim(init_type) == 'sinx_band' ) .and. &
       init_kx_max == 1 .and. init_kx > 1 .and. init_kx_min == 1 ) then
    init_kx_max = init_kx
  endif

  this%init_kx = init_kx
  this%init_ky = init_ky
  this%init_kz = init_kz
  this%init_kx_min = init_kx_min
  this%init_kx_max = init_kx_max
  this%init_phase_step = init_phase_step
  this%init_band_norm = init_band_norm
  this%init_lambda = init_lambda
  this%init_phase = init_phase
  this%init_x0 = init_x0
  this%init_y0 = init_y0
  this%init_z0 = init_z0
  this%init_dir = init_dir
  this%init_w_long = init_w_long
  this%init_w_trans1 = init_w_trans1
  this%init_w_trans2 = init_w_trans2
  this%init_vphase = init_vphase
  this%bc_type(p_lower,1) = canonical_axion_bc( bc_x_lower )
  this%bc_type(p_upper,1) = canonical_axion_bc( bc_x_upper )
  this%bc_type(p_lower,2) = canonical_axion_bc( bc_y_lower )
  this%bc_type(p_upper,2) = canonical_axion_bc( bc_y_upper )
  this%bc_type(p_lower,3) = canonical_axion_bc( bc_z_lower )
  this%bc_type(p_upper,3) = canonical_axion_bc( bc_z_upper )

  ndump_fac_all(p_full) = ndump_fac
  ndump_fac_all(p_savg) = ndump_fac_ave
  ndump_fac_all(p_senv) = ndump_fac_ave
  ndump_fac_all(p_line) = ndump_fac_lineout
  ndump_fac_all(p_slice) = ndump_fac_lineout

  call this%diag%init_report_quants()
  call new( this%diag%reports, reports, this%diag%report_quants, ndump_fac_all, n_ave, n_tavg, prec, p_x_dim, ierr )
  if ( ierr /= 0 ) then
    print *, '(*error*) Invalid axion diagnostic report'
    stop
  endif
end subroutine read_input_axion_field_impl

! Allocate and initialize the axion-field state, diagnostics, and perturbation fields.
subroutine init_axion_field_impl( this, g_space, grid, no_co, gc_num_in, dt_init )
  class(t_axion_field), intent(inout) :: this
  type(t_space), intent(in) :: g_space
  class(t_grid), intent(in) :: grid
  class(t_node_conf), intent(in) :: no_co
  integer, dimension(2,p_x_dim), intent(in) :: gc_num_in
  real(p_double), intent(in) :: dt_init

  real(p_double), dimension(p_x_dim) :: ldx
  integer :: i1, i2, i3
  integer :: gx, gy, gz
  real(p_double) :: x, y, z, twopi, phase, amp
  real(p_double) :: lx_band, kx_band, rhs_band, omega_band, cosfac_band, mode_amp, phi_sum
  real(p_double) :: dir(3), dir_norm, k0, time_old, time_new, phi_old, phi_new
  real(p_double) :: s, q1, q2, eta, gauss_arg, w_long, w_t1, w_t2
  real(p_double) :: rx, ry, rz, e1(3), e2(3), ref(3), ref_norm
  real(p_double) :: omega_d, vphase_d, vgroup_d, v_env, ds_eff, arg_disp, denom
  integer :: idim, imode, nmodes

  twopi = 2.0_p_double * acos(-1.0_p_double)

  this%x_dim = p_x_dim
  this%nx(1:p_x_dim) = grid%my_nx(3,1:p_x_dim)
  this%gc_num(:,1:p_x_dim) = gc_num_in(:,1:p_x_dim)

  ldx = ( xmax(g_space) - xmin(g_space) ) / grid%g_nx(1:p_x_dim)
  this%dx(1:p_x_dim) = ldx(1:p_x_dim)
  this%xmin_loc(1:p_x_dim) = xmin(g_space)
  this%xmax_loc(1:p_x_dim) = xmax(g_space)

  call this%a%new( p_x_dim, 1, this%nx(1:p_x_dim), this%gc_num(:,1:p_x_dim), ldx, .true. )
  call this%pi%new( p_x_dim, 1, this%nx(1:p_x_dim), this%gc_num(:,1:p_x_dim), ldx, .true. )
  call this%pert%e%new( p_x_dim, 3, this%nx(1:p_x_dim), this%gc_num(:,1:p_x_dim), ldx, .true. )
  call this%pert%b%new( p_x_dim, 3, this%nx(1:p_x_dim), this%gc_num(:,1:p_x_dim), ldx, .true. )
  call this%diag_edotb%new( p_x_dim, 1, this%nx(1:p_x_dim), this%gc_num(:,1:p_x_dim), ldx, .true. )
  call this%diag_emsrcmag%new( p_x_dim, 1, this%nx(1:p_x_dim), this%gc_num(:,1:p_x_dim), ldx, .true. )
  this%pert%if_initialized = .true.

  amp = this%init_ampl

  select case ( trim(this%init_type) )
  case ( 'zero' )
    ! already zeroed
  case ( 'sin_mode', 'sinx' )
    select case ( p_x_dim )
    case (1)
      do i1 = 1, this%nx(1)
        gx = grid%my_nx(p_lower,1) + i1 - 1
        x = xmin(g_space,1) + ( gx - 0.5_p_double ) * ldx(1)
        phase = twopi * real(this%init_kx,p_double) * ( x - xmin(g_space,1) ) / ( xmax(g_space,1) - xmin(g_space,1) )
        this%a%f1(1,i1) = real( amp * sin( phase ), p_k_fld )
        this%pi%f1(1,i1) = this%a%f1(1,i1)
      enddo
    case (2)
      do i2 = 1, this%nx(2)
        gy = grid%my_nx(p_lower,2) + i2 - 1
        y = xmin(g_space,2) + ( gy - 0.5_p_double ) * ldx(2)
        do i1 = 1, this%nx(1)
          gx = grid%my_nx(p_lower,1) + i1 - 1
          x = xmin(g_space,1) + ( gx - 0.5_p_double ) * ldx(1)
          phase = twopi * real(this%init_kx,p_double) * ( x - xmin(g_space,1) ) / ( xmax(g_space,1) - xmin(g_space,1) ) + &
                  twopi * real(this%init_ky,p_double) * ( y - xmin(g_space,2) ) / ( xmax(g_space,2) - xmin(g_space,2) )
          this%a%f2(1,i1,i2) = real( amp * sin( phase ), p_k_fld )
          this%pi%f2(1,i1,i2) = this%a%f2(1,i1,i2)
        enddo
      enddo
    case (3)
      do i3 = 1, this%nx(3)
        gz = grid%my_nx(p_lower,3) + i3 - 1
        z = xmin(g_space,3) + ( gz - 0.5_p_double ) * ldx(3)
        do i2 = 1, this%nx(2)
          gy = grid%my_nx(p_lower,2) + i2 - 1
          y = xmin(g_space,2) + ( gy - 0.5_p_double ) * ldx(2)
          do i1 = 1, this%nx(1)
            gx = grid%my_nx(p_lower,1) + i1 - 1
            x = xmin(g_space,1) + ( gx - 0.5_p_double ) * ldx(1)
            phase = twopi * real(this%init_kx,p_double) * ( x - xmin(g_space,1) ) / ( xmax(g_space,1) - xmin(g_space,1) ) + &
                    twopi * real(this%init_ky,p_double) * ( y - xmin(g_space,2) ) / ( xmax(g_space,2) - xmin(g_space,2) ) + &
                    twopi * real(this%init_kz,p_double) * ( z - xmin(g_space,3) ) / ( xmax(g_space,3) - xmin(g_space,3) )
            this%a%f3(1,i1,i2,i3) = real( amp * sin( phase ), p_k_fld )
            this%pi%f3(1,i1,i2,i3) = this%a%f3(1,i1,i2,i3)
          enddo
        enddo
      enddo
    end select
  case ( 'sin_band', 'sinx_band' )
    ! Seed a longitudinal Fourier band along x1. Each mode is initialized with
    ! its own discrete half-step correction so the solver must not later apply
    ! the single-mode first-advance correction to this init_type.
    if ( this%init_kx_min < 1 ) then
      write(err_buf__,*) 'sin_band requires init_kx_min >= 1';call err__("axion_field/os-axion-field-core.f03",247)
      call abort_program( p_err_invalid )
    endif
    if ( this%init_kx_max < this%init_kx_min ) then
      write(err_buf__,*) 'sin_band requires init_kx_max >= init_kx_min';call err__("axion_field/os-axion-field-core.f03",251)
      call abort_program( p_err_invalid )
    endif
    lx_band = xmax(g_space,1) - xmin(g_space,1)
    if ( lx_band <= 0.0_p_double ) then
      write(err_buf__,*) 'sin_band requires positive x1 domain length';call err__("axion_field/os-axion-field-core.f03",256)
      call abort_program( p_err_invalid )
    endif
    nmodes = this%init_kx_max - this%init_kx_min + 1
    mode_amp = amp
    if ( this%init_band_norm .and. nmodes > 0 ) then
      mode_amp = amp / sqrt( real(nmodes,p_double) )
    endif

    select case ( p_x_dim )
    case (1)
      do i1 = 1, this%nx(1)
        gx = grid%my_nx(p_lower,1) + i1 - 1
        x = xmin(g_space,1) + ( gx - 0.5_p_double ) * ldx(1)
        phi_sum = 0.0_p_double
        do imode = this%init_kx_min, this%init_kx_max
          kx_band = twopi * real(imode,p_double) / lx_band
          rhs_band = this%ma*this%ma + &
                     ( 2.0_p_double / this%dx(1) * sin( 0.5_p_double * kx_band * this%dx(1) ) )**2
          if ( dt_init > 0.0_p_double ) then
            omega_band = 2.0_p_double / dt_init * &
                         asin( min( 1.0_p_double, 0.5_p_double * dt_init * sqrt( max(rhs_band,0.0_p_double) ) ) )
            cosfac_band = cos( 0.5_p_double * omega_band * dt_init )
          else
            cosfac_band = 1.0_p_double
          endif
          if ( abs(cosfac_band) < 1.0e-12_p_double ) cosfac_band = 1.0_p_double
          phase = kx_band * ( x - xmin(g_space,1) ) + this%init_phase + &
                  real(imode - this%init_kx_min,p_double) * this%init_phase_step
          phi_sum = phi_sum + mode_amp * cosfac_band * sin( phase )
        enddo
        this%a%f1(1,i1) = real( phi_sum, p_k_fld )
        this%pi%f1(1,i1) = real( phi_sum, p_k_fld )
      enddo
    case (2)
      do i2 = 1, this%nx(2)
        do i1 = 1, this%nx(1)
          gx = grid%my_nx(p_lower,1) + i1 - 1
          x = xmin(g_space,1) + ( gx - 0.5_p_double ) * ldx(1)
          phi_sum = 0.0_p_double
          do imode = this%init_kx_min, this%init_kx_max
            kx_band = twopi * real(imode,p_double) / lx_band
            rhs_band = this%ma*this%ma + &
                       ( 2.0_p_double / this%dx(1) * sin( 0.5_p_double * kx_band * this%dx(1) ) )**2
            if ( dt_init > 0.0_p_double ) then
              omega_band = 2.0_p_double / dt_init * &
                           asin( min( 1.0_p_double, 0.5_p_double * dt_init * sqrt( max(rhs_band,0.0_p_double) ) ) )
              cosfac_band = cos( 0.5_p_double * omega_band * dt_init )
            else
              cosfac_band = 1.0_p_double
            endif
            if ( abs(cosfac_band) < 1.0e-12_p_double ) cosfac_band = 1.0_p_double
            phase = kx_band * ( x - xmin(g_space,1) ) + this%init_phase + &
                    real(imode - this%init_kx_min,p_double) * this%init_phase_step
            phi_sum = phi_sum + mode_amp * cosfac_band * sin( phase )
          enddo
          this%a%f2(1,i1,i2) = real( phi_sum, p_k_fld )
          this%pi%f2(1,i1,i2) = real( phi_sum, p_k_fld )
        enddo
      enddo
    case (3)
      do i3 = 1, this%nx(3)
        do i2 = 1, this%nx(2)
          do i1 = 1, this%nx(1)
            gx = grid%my_nx(p_lower,1) + i1 - 1
            x = xmin(g_space,1) + ( gx - 0.5_p_double ) * ldx(1)
            phi_sum = 0.0_p_double
            do imode = this%init_kx_min, this%init_kx_max
              kx_band = twopi * real(imode,p_double) / lx_band
              rhs_band = this%ma*this%ma + &
                         ( 2.0_p_double / this%dx(1) * sin( 0.5_p_double * kx_band * this%dx(1) ) )**2
              if ( dt_init > 0.0_p_double ) then
                omega_band = 2.0_p_double / dt_init * &
                             asin( min( 1.0_p_double, 0.5_p_double * dt_init * sqrt( max(rhs_band,0.0_p_double) ) ) )
                cosfac_band = cos( 0.5_p_double * omega_band * dt_init )
              else
                cosfac_band = 1.0_p_double
              endif
              if ( abs(cosfac_band) < 1.0e-12_p_double ) cosfac_band = 1.0_p_double
              phase = kx_band * ( x - xmin(g_space,1) ) + this%init_phase + &
                      real(imode - this%init_kx_min,p_double) * this%init_phase_step
              phi_sum = phi_sum + mode_amp * cosfac_band * sin( phase )
            enddo
            this%a%f3(1,i1,i2,i3) = real( phi_sum, p_k_fld )
            this%pi%f3(1,i1,i2,i3) = real( phi_sum, p_k_fld )
          enddo
        enddo
      enddo
    end select

  case ( 'travelling_pulse', 'traveling_pulse', 'gaussian_pulse' )
    if ( this%init_lambda <= 0.0_p_double ) then
      write(err_buf__,*) 'travelling_pulse requires init_lambda > 0';call err__("axion_field/os-axion-field-core.f03",348)
      call abort_program( p_err_invalid )
    endif
    dir = this%init_dir
    if ( p_x_dim < 3 ) dir(p_x_dim+1:3) = 0.0_p_double
    dir_norm = sqrt( sum( dir**2 ) )
    if ( dir_norm <= 0.0_p_double ) then
      dir = (/ 1.0_p_double, 0.0_p_double, 0.0_p_double /)
      dir_norm = 1.0_p_double
    endif
    dir = dir / dir_norm
    w_long = this%init_w_long
    if ( w_long <= 0.0_p_double ) w_long = this%init_lambda
    w_t1 = this%init_w_trans1
    if ( w_t1 <= 0.0_p_double ) w_t1 = w_long
    w_t2 = this%init_w_trans2
    if ( w_t2 <= 0.0_p_double ) w_t2 = w_t1
    k0 = twopi / this%init_lambda
    time_old = -0.5_p_double * dt_init
    time_new = 0.5_p_double * dt_init

    ds_eff = 0.0_p_double
    do idim = 1, p_x_dim
      ds_eff = ds_eff + ( dir(idim) * ldx(idim) )**2
    enddo
    ds_eff = sqrt( max( ds_eff, tiny(1.0_p_double) ) )
    if ( this%init_lambda / ds_eff < 8.0_p_double ) then
      print *, '(*error*) travelling_pulse / gaussian_pulse under-resolved on this grid'
      print *, '          lambda / ds_eff = ', this%init_lambda / ds_eff, '  (require >= 8)'
      call abort_program( p_err_invalid )
    endif

    arg_disp = 0.25_p_double * dt_init * dt_init * this%ma * this%ma
    do idim = 1, p_x_dim
      arg_disp = arg_disp + ( dt_init / ldx(idim) )**2 * sin( 0.5_p_double * k0 * dir(idim) * ldx(idim) )**2
    enddo
    arg_disp = min( max( arg_disp, 0.0_p_double ), 1.0_p_double - tiny(1.0_p_double) )
    omega_d = 2.0_p_double * asin( sqrt( arg_disp ) ) / dt_init
    if ( k0 > tiny(1.0_p_double) ) then
      vphase_d = omega_d / k0
    else
      vphase_d = 0.0_p_double
    endif
    denom = sin( omega_d * dt_init )
    if ( abs(denom) > tiny(1.0_p_double) ) then
      vgroup_d = 0.0_p_double
      do idim = 1, p_x_dim
        vgroup_d = vgroup_d + dir(idim) * sin( k0 * dir(idim) * ldx(idim) ) / ldx(idim)
      enddo
      vgroup_d = dt_init * vgroup_d / denom
    else
      vgroup_d = vphase_d
    endif
    v_env = this%init_vphase
    if ( v_env <= 0.0_p_double ) v_env = vgroup_d

    if ( p_x_dim == 1 ) then
      e1 = (/ 0.0_p_double, 1.0_p_double, 0.0_p_double /)
      e2 = (/ 0.0_p_double, 0.0_p_double, 1.0_p_double /)
    else if ( p_x_dim == 2 ) then
      e1 = (/ -dir(2), dir(1), 0.0_p_double /)
      ref_norm = sqrt( e1(1)**2 + e1(2)**2 + e1(3)**2 )
      if ( ref_norm <= 0.0_p_double ) e1 = (/ 0.0_p_double, 1.0_p_double, 0.0_p_double /)
      e1 = e1 / sqrt( max( sum(e1**2), tiny(1.0_p_double) ) )
      e2 = (/ 0.0_p_double, 0.0_p_double, 1.0_p_double /)
    else
      ref = (/ 0.0_p_double, 0.0_p_double, 1.0_p_double /)
      if ( abs(dir(3)) > 0.9_p_double ) ref = (/ 0.0_p_double, 1.0_p_double, 0.0_p_double /)
      e1(1) = dir(2)*ref(3) - dir(3)*ref(2)
      e1(2) = dir(3)*ref(1) - dir(1)*ref(3)
      e1(3) = dir(1)*ref(2) - dir(2)*ref(1)
      e1 = e1 / sqrt( max( sum(e1**2), tiny(1.0_p_double) ) )
      e2(1) = dir(2)*e1(3) - dir(3)*e1(2)
      e2(2) = dir(3)*e1(1) - dir(1)*e1(3)
      e2(3) = dir(1)*e1(2) - dir(2)*e1(1)
      e2 = e2 / sqrt( max( sum(e2**2), tiny(1.0_p_double) ) )
    endif

    select case ( p_x_dim )
    case (1)
      do i1 = 1, this%nx(1)
        gx = grid%my_nx(p_lower,1) + i1 - 1
        x = xmin(g_space,1) + ( gx - 0.5_p_double ) * ldx(1)
        rx = x - this%init_x0; ry = 0.0_p_double; rz = 0.0_p_double
        s = rx*dir(1)
        q1 = 0.0_p_double; q2 = 0.0_p_double
        eta = s - v_env*time_old
        gauss_arg = (eta/w_long)**2
        phi_old = amp * exp( -gauss_arg ) * sin( k0*s - omega_d*time_old + this%init_phase )
        eta = s - v_env*time_new
        gauss_arg = (eta/w_long)**2
        phi_new = amp * exp( -gauss_arg ) * sin( k0*s - omega_d*time_new + this%init_phase )
        this%a%f1(1,i1) = real( phi_old, p_k_fld )
        this%pi%f1(1,i1) = real( phi_new, p_k_fld )
      enddo
    case (2)
      do i2 = 1, this%nx(2)
        gy = grid%my_nx(p_lower,2) + i2 - 1
        y = xmin(g_space,2) + ( gy - 0.5_p_double ) * ldx(2)
        do i1 = 1, this%nx(1)
          gx = grid%my_nx(p_lower,1) + i1 - 1
          x = xmin(g_space,1) + ( gx - 0.5_p_double ) * ldx(1)
          rx = x - this%init_x0; ry = y - this%init_y0; rz = 0.0_p_double
          s = rx*dir(1) + ry*dir(2)
          q1 = rx*e1(1) + ry*e1(2)
          q2 = 0.0_p_double
          eta = s - v_env*time_old
          gauss_arg = (eta/w_long)**2 + (q1/w_t1)**2
          phi_old = amp * exp( -gauss_arg ) * sin( k0*s - omega_d*time_old + this%init_phase )
          eta = s - v_env*time_new
          gauss_arg = (eta/w_long)**2 + (q1/w_t1)**2
          phi_new = amp * exp( -gauss_arg ) * sin( k0*s - omega_d*time_new + this%init_phase )
          this%a%f2(1,i1,i2) = real( phi_old, p_k_fld )
          this%pi%f2(1,i1,i2) = real( phi_new, p_k_fld )
        enddo
      enddo
    case (3)
      do i3 = 1, this%nx(3)
        gz = grid%my_nx(p_lower,3) + i3 - 1
        z = xmin(g_space,3) + ( gz - 0.5_p_double ) * ldx(3)
        do i2 = 1, this%nx(2)
          gy = grid%my_nx(p_lower,2) + i2 - 1
          y = xmin(g_space,2) + ( gy - 0.5_p_double ) * ldx(2)
          do i1 = 1, this%nx(1)
            gx = grid%my_nx(p_lower,1) + i1 - 1
            x = xmin(g_space,1) + ( gx - 0.5_p_double ) * ldx(1)
            rx = x - this%init_x0; ry = y - this%init_y0; rz = z - this%init_z0
            s = rx*dir(1) + ry*dir(2) + rz*dir(3)
            q1 = rx*e1(1) + ry*e1(2) + rz*e1(3)
            q2 = rx*e2(1) + ry*e2(2) + rz*e2(3)
            eta = s - v_env*time_old
            gauss_arg = (eta/w_long)**2 + (q1/w_t1)**2 + (q2/w_t2)**2
            phi_old = amp * exp( -gauss_arg ) * sin( k0*s - omega_d*time_old + this%init_phase )
            eta = s - v_env*time_new
            gauss_arg = (eta/w_long)**2 + (q1/w_t1)**2 + (q2/w_t2)**2
            phi_new = amp * exp( -gauss_arg ) * sin( k0*s - omega_d*time_new + this%init_phase )
            this%a%f3(1,i1,i2,i3) = real( phi_old, p_k_fld )
            this%pi%f3(1,i1,i2,i3) = real( phi_new, p_k_fld )
          enddo
        enddo
      enddo
    end select
  case default
    write(err_buf__,*) 'Unsupported axion_field init_type = ', trim(this%init_type);call err__("axion_field/os-axion-field-core.f03",491)
    call abort_program( p_err_invalid )
  end select

  call apply_axion_bc_scalar( this, this%a )
  call apply_axion_bc_scalar( this, this%pi )
  call apply_axion_bc_vector( this, this%pert%e )
  call apply_axion_bc_vector( this, this%pert%b )

  this%if_initialized = .true.
  this%first_advance = .true.
  call this%diag%init()
end subroutine init_axion_field_impl

! Release axion-field storage and diagnostic resources.
subroutine cleanup_axion_field_impl( this )
  class(t_axion_field), intent(inout) :: this
  call this%a%cleanup()
  call this%pi%cleanup()
  call this%diag_edotb%cleanup()
  call this%diag_emsrcmag%cleanup()
  if ( this%pert%if_initialized ) then
    call this%pert%e%cleanup()
    call this%pert%b%cleanup()
    this%pert%if_initialized = .false.
  endif
  call this%diag%cleanup()
  this%if_initialized = .false.
  this%first_advance = .true.
end subroutine cleanup_axion_field_impl

! Shift the axion-field state when the moving window advances.
subroutine move_window_axion_field_impl( this, g_space )
  class(t_axion_field), intent(inout) :: this
  type(t_space), intent(in) :: g_space
  if ( .not. this%if_initialized ) return
  call move_window( this%a, g_space )
  call move_window( this%pi, g_space )
  call move_window( this%diag_edotb, g_space )
  call move_window( this%diag_emsrcmag, g_space )
  if ( this%pert%if_initialized ) then
    call move_window( this%pert%e, g_space )
    call move_window( this%pert%b, g_space )
  endif
end subroutine move_window_axion_field_impl

! Exchange guard cells and apply physical boundary conditions for the axion field.
subroutine update_boundary_axion_field_impl( this, no_co, send_msg, recv_msg )
  class(t_axion_field), intent(inout) :: this
  class(t_node_conf), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg
  if ( .not. this%if_initialized ) return
  call apply_axion_bc_scalar( this, this%a )
  call apply_axion_bc_scalar( this, this%pi )
  call apply_axion_bc_scalar( this, this%diag_edotb )
  call apply_axion_bc_scalar( this, this%diag_emsrcmag )
  if ( this%pert%if_initialized ) then
    call apply_axion_bc_vector( this, this%pert%e )
    call apply_axion_bc_vector( this, this%pert%b )
  endif
end subroutine update_boundary_axion_field_impl

! Print a summary of the coherent axion-field algorithm and settings.
subroutine list_algorithm_axion_field_impl( this )
  class(t_axion_field), intent(in) :: this
  print *, 'Axion field: coherent KG field with separate perturbative EM container (pass 2)'
  print *, '  g_agg is used as an internal simulation-unit coupling in the normalized axion source terms'
  print *, '  ma is expected in normalized frequency units (typically ma/omega_ref)'
  print *, '  for OSIRIS laser drivers, E0(code) ~ omega0 * a0 so A0_est ~ max|E0_perp| / omega0'
  print *, '  init_type = ', trim(this%init_type)
  if ( trim(this%init_type) == 'sin_mode' .or. trim(this%init_type) == 'sinx' ) then
    print *, '  init_kx,ky,kz = ', this%init_kx, this%init_ky, this%init_kz
    print *, '  init_ampl = ', this%init_ampl
    print *, '  init_phase = ', this%init_phase
  endif
  if ( trim(this%init_type) == 'sin_band' .or. trim(this%init_type) == 'sinx_band' ) then
    print *, '  init_kx_min,max = ', this%init_kx_min, this%init_kx_max
    print *, '  init_ampl per mode = ', this%init_ampl
    print *, '  init_phase = ', this%init_phase
    print *, '  init_phase_step = ', this%init_phase_step
    print *, '  init_band_norm = ', this%init_band_norm
  endif
  print *, '  bc x(lower,upper) = ', trim(this%bc_type(p_lower,1)), trim(this%bc_type(p_upper,1))
  if ( p_x_dim >= 2 ) print *, '  bc y(lower,upper) = ', trim(this%bc_type(p_lower,2)), trim(this%bc_type(p_upper,2))
  if ( p_x_dim >= 3 ) print *, '  bc z(lower,upper) = ', trim(this%bc_type(p_lower,3)), trim(this%bc_type(p_upper,3))
  if ( trim(this%init_type) == 'travelling_pulse' .or. trim(this%init_type) == 'traveling_pulse' .or. trim(this%init_type) == 'gaussian_pulse' ) then
    print *, '  init_lambda = ', this%init_lambda
    print *, '  init_phase  = ', this%init_phase
    print *, '  init_x0,y0,z0 = ', this%init_x0, this%init_y0, this%init_z0
    print *, '  init_dir    = ', this%init_dir
    print *, '  init_w_long = ', this%init_w_long
    print *, '  init_w_trans1 = ', this%init_w_trans1
    print *, '  init_w_trans2 = ', this%init_w_trans2
    print *, '  init_vphase = ', this%init_vphase
    if ( this%init_vphase <= 0.0_p_double ) print *, '  init_vphase <= 0 -> auto discrete group velocity for envelope'
  endif
end subroutine list_algorithm_axion_field_impl

end module m_axion_field_core_impl

! Forward the axion-field input read to the compiled implementation module.
subroutine read_input_axion_field( this, input_file )
  use m_axion_field_define, only : t_axion_field
  use m_input_file, only : t_input_file
  use m_axion_field_core_impl, only : read_input_axion_field_impl
  implicit none
  class(t_axion_field), intent(inout) :: this
  class(t_input_file), intent(inout) :: input_file
  call read_input_axion_field_impl( this, input_file )
end subroutine read_input_axion_field

! Forward axion-field initialization to the compiled implementation module.
subroutine init_axion_field( this, g_space, grid, no_co, gc_num_in, dt_init )
  use m_parameters, only : p_x_dim, p_double
  use m_axion_field_define, only : t_axion_field
  use m_space, only : t_space
  use m_grid_define, only : t_grid
  use m_node_conf, only : t_node_conf
  use m_axion_field_core_impl, only : init_axion_field_impl
  implicit none
  class(t_axion_field), intent(inout) :: this
  type(t_space), intent(in) :: g_space
  class(t_grid), intent(in) :: grid
  class(t_node_conf), intent(in) :: no_co
  integer, dimension(2,p_x_dim), intent(in) :: gc_num_in
  real(p_double), intent(in) :: dt_init
  call init_axion_field_impl( this, g_space, grid, no_co, gc_num_in, dt_init )
end subroutine init_axion_field

! Forward axion-field cleanup to the compiled implementation module.
subroutine cleanup_axion_field( this )
  use m_axion_field_define, only : t_axion_field
  use m_axion_field_core_impl, only : cleanup_axion_field_impl
  implicit none
  class(t_axion_field), intent(inout) :: this
  call cleanup_axion_field_impl( this )
end subroutine cleanup_axion_field

! Forward the axion-field moving-window update to the compiled implementation module.
subroutine move_window_axion_field( this, g_space )
  use m_axion_field_define, only : t_axion_field
  use m_space, only : t_space
  use m_axion_field_core_impl, only : move_window_axion_field_impl
  implicit none
  class(t_axion_field), intent(inout) :: this
  type(t_space), intent(in) :: g_space
  call move_window_axion_field_impl( this, g_space )
end subroutine move_window_axion_field

! Forward the axion-field boundary update to the compiled implementation module.
subroutine update_boundary_axion_field( this, no_co, send_msg, recv_msg )
  use m_axion_field_define, only : t_axion_field
  use m_node_conf, only : t_node_conf
  use m_vdf_comm, only : t_vdf_msg
  use m_axion_field_core_impl, only : update_boundary_axion_field_impl
  implicit none
  class(t_axion_field), intent(inout) :: this
  class(t_node_conf), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg
  call update_boundary_axion_field_impl( this, no_co, send_msg, recv_msg )
end subroutine update_boundary_axion_field

! Forward the axion-field algorithm listing to the compiled implementation module.
subroutine list_algorithm_axion_field( this )
  use m_axion_field_define, only : t_axion_field
  use m_axion_field_core_impl, only : list_algorithm_axion_field_impl
  implicit none
  class(t_axion_field), intent(in) :: this
  call list_algorithm_axion_field_impl( this )
end subroutine list_algorithm_axion_field
