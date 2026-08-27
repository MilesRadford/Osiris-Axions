# 1 "axion_field/os-axion-field-diagnostics.f03"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "axion_field/os-axion-field-diagnostics.f03"
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
# 2 "axion_field/os-axion-field-diagnostics.f03" 2
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
# 3 "axion_field/os-axion-field-diagnostics.f03" 2

module m_axion_field_diagnostics_impl

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
  use m_axion_field_diag_define
  use m_axion_field_solver_impl, only : phi_center_1d, phi_center_2d, phi_center_3d, &
                                        phi_t_1d, phi_t_2d, phi_t_3d

  implicit none
  private

  public :: avail_report_quants_axion_impl
  public :: init_report_quants_axion_impl
  public :: read_input_diag_axion_impl
  public :: init_diag_axion_impl
  public :: cleanup_diag_axion_impl
  public :: report_axion_field_impl
  public :: report_energy_axion_field_impl

contains

! Return the number of supported axion-field diagnostic quantities.
function avail_report_quants_axion_impl( this ) result(nq)
  class(t_axion_diag), intent(in) :: this
  integer :: nq
  nq = size( p_axion_report_quants )
end function avail_report_quants_axion_impl

! Populate the list of valid axion-field diagnostic quantity names.
subroutine init_report_quants_axion_impl( this )
  class(t_axion_diag), intent(inout) :: this
  if ( .not. associated( this%report_quants ) ) then
    allocate( this%report_quants( this%avail_report_quants() ) )
  endif
  this%report_quants(1:size(p_axion_report_quants)) = p_axion_report_quants
end subroutine init_report_quants_axion_impl

! Read the axion-field diagnostics namelist.
subroutine read_input_diag_axion_impl( this, input_file )
  class(t_axion_diag), intent(inout) :: this
  class(t_input_file), intent(inout) :: input_file
  integer :: ndump_fac, ndump_fac_ave, ndump_fac_lineout, prec, n_tavg, ierr
  integer, dimension(p_x_dim) :: n_ave
  integer, dimension(p_n_report_type) :: ndump_fac_all
  character( len = p_max_reports_len ), dimension( p_max_reports ) :: reports
  namelist /nl_diag_axion_field/ ndump_fac, ndump_fac_ave, ndump_fac_lineout, prec, n_ave, n_tavg, reports

  ndump_fac = 0
  ndump_fac_ave = 0
  ndump_fac_lineout = 0
  prec = p_single
  n_ave = -1
  n_tavg = -1
  reports = '-'

  call get_namelist( input_file, 'nl_diag_axion_field', ierr )
  if ( ierr == 0 ) then
    read (input_file%nml_text, nml = nl_diag_axion_field, iostat = ierr)
    if ( ierr /= 0 ) then
      print *, 'Error reading diag_axion_field parameters'
      stop
    endif
  endif

  ndump_fac_all(p_full) = ndump_fac
  ndump_fac_all(p_savg) = ndump_fac_ave
  ndump_fac_all(p_senv) = ndump_fac_ave
  ndump_fac_all(p_line) = ndump_fac_lineout
  ndump_fac_all(p_slice) = ndump_fac_lineout

  call this%init_report_quants()
  call new( this%reports, reports, this%report_quants, ndump_fac_all, n_ave, n_tavg, prec, p_x_dim, ierr )
  if ( ierr /= 0 ) then
    print *, '(*error*) Invalid axion diagnostic report'
    stop
  endif
end subroutine read_input_diag_axion_impl

! Build the axion-field diagnostic report descriptors and event handle.
subroutine init_diag_axion_impl( this )
  class(t_axion_diag), intent(inout) :: this
  type(t_vdf_report), pointer :: report

  report => this%reports
  do
    if ( .not. associated( report ) ) exit
    report%xname = (/'x1', 'x2', 'x3'/)
    report%xlabel = (/'x_1', 'x_2', 'x_3'/)
    report%xunits = (/'c / \omega_p', 'c / \omega_p', 'c / \omega_p'/)
    report%time_units = '1 / \omega_p'
    report%dt = 1.0_p_double
    report%fileLabel = ''
    report%basePath = trim(path_mass) // 'FLD' // p_dir_sep

    select case ( report%quant )
    case ( p_phi )
      report%label = '\phi'
      report%units = 'axion field'
      report%offset_t = 0.0_p_double
      report%offset_x = (/0.5_p_double, 0.5_p_double, 0.5_p_double/)
    case ( p_phidot )
      report%label = '\partial_t \phi'
      report%units = 'axion field / (1/\omega_p)'
      report%offset_t = 0.0_p_double
      report%offset_x = (/0.5_p_double, 0.5_p_double, 0.5_p_double/)
    case ( p_aold )
      report%label = '\phi^{n-1/2}'
      report%units = 'axion field'
      report%offset_t = -0.5_p_double
      report%offset_x = (/0.5_p_double, 0.5_p_double, 0.5_p_double/)
    case ( p_anew )
      report%label = '\phi^{n+1/2}'
      report%units = 'axion field'
      report%offset_t = 0.5_p_double
      report%offset_x = (/0.5_p_double, 0.5_p_double, 0.5_p_double/)
    case ( p_ae1, p_ae2, p_ae3 )
      report%label = 'E_'//char(ichar('0') + 1 + report%quant - p_ae1)//'^{ax}'
      report%units = 'm_e c \omega_p e^{-1}'
      report%offset_t = 0.0_p_double
      report%offset_x = (/0.0_p_double, 0.0_p_double, 0.0_p_double/)
      report%offset_x(report%quant-p_ae1+1) = 0.5_p_double
    case ( p_ab1, p_ab2, p_ab3 )
      report%label = 'B_'//char(ichar('0') + 1 + report%quant - p_ab1)//'^{ax}'
      report%units = 'm_e c \omega_p e^{-1}'
      report%offset_t = 0.0_p_double
      report%offset_x = (/0.5_p_double, 0.5_p_double, 0.5_p_double/)
      report%offset_x(report%quant-p_ab1+1) = 0.0_p_double
    case ( p_edotb )
      report%label = 'E \cdot B'
      report%units = 'field^2'
      report%offset_t = 0.0_p_double
      report%offset_x = (/0.5_p_double, 0.5_p_double, 0.5_p_double/)
    case ( p_emsrcmag )
      report%label = '|S_{E1}|'
      report%units = 'EM source'
      report%offset_t = 0.0_p_double
      report%offset_x = (/0.5_p_double, 0.5_p_double, 0.5_p_double/)
    end select

    report => report%next
  enddo

  if ( diag_axion_ev == 0 ) diag_axion_ev = create_event('axion field diagnostics')
end subroutine init_diag_axion_impl

! Release axion-field diagnostic report storage.
subroutine cleanup_diag_axion_impl( this )
  class(t_axion_diag), intent(inout) :: this
  call cleanup( this%reports )
end subroutine cleanup_diag_axion_impl

! Write the configured axion-field diagnostic datasets for the current dump step.
subroutine report_axion_field_impl( this, g_space, grid, no_co, tstep, t )
  class(t_axion_field), intent(inout) :: this
  type(t_space), intent(in) :: g_space
  class(t_grid), intent(in) :: grid
  class(t_node_conf), intent(in) :: no_co
  type(t_time_step), intent(in) :: tstep
  real(p_double), intent(in) :: t
  type(t_vdf_report), pointer :: rep
  type(t_vdf) :: tmp
  integer :: i1,i2,i3
  real(p_double) :: idt

  if ( .not. this%if_initialized ) return
  call begin_event( diag_axion_ev )
  rep => this%diag%reports
  do
    if ( .not. associated(rep) ) exit
    if ( if_report( rep, tstep ) ) then
      select case( rep%quant )
      case ( p_phi )
        call tmp%new( this%a, f_dim = 1 )
        select case ( p_x_dim )
        case (1)
          do i1=1,this%nx(1)
            tmp%f1(1,i1)=real(0.5_p_double*(real(this%a%f1(1,i1),p_double)+real(this%pi%f1(1,i1),p_double)),p_k_fld)
          enddo
        case (2)
          do i2=1,this%nx(2); do i1=1,this%nx(1)
            tmp%f2(1,i1,i2)=real(0.5_p_double*(real(this%a%f2(1,i1,i2),p_double)+real(this%pi%f2(1,i1,i2),p_double)),p_k_fld)
          enddo; enddo
        case (3)
          do i3=1,this%nx(3); do i2=1,this%nx(2); do i1=1,this%nx(1)
            tmp%f3(1,i1,i2,i3)=real(0.5_p_double*(real(this%a%f3(1,i1,i2,i3),p_double)+real(this%pi%f3(1,i1,i2,i3),p_double)),p_k_fld)
          enddo; enddo; enddo
        end select
        call report_vdf( rep, tmp, 1, g_space, grid, no_co, tstep, t )
        call tmp%cleanup()
      case ( p_phidot )
        idt = 1.0_p_double / dt(tstep)
        call tmp%new( this%a, f_dim = 1 )
        select case ( p_x_dim )
        case (1)
          do i1=1,this%nx(1)
            tmp%f1(1,i1)=real((real(this%pi%f1(1,i1),p_double)-real(this%a%f1(1,i1),p_double))*idt,p_k_fld)
          enddo
        case (2)
          do i2=1,this%nx(2); do i1=1,this%nx(1)
            tmp%f2(1,i1,i2)=real((real(this%pi%f2(1,i1,i2),p_double)-real(this%a%f2(1,i1,i2),p_double))*idt,p_k_fld)
          enddo; enddo
        case (3)
          do i3=1,this%nx(3); do i2=1,this%nx(2); do i1=1,this%nx(1)
            tmp%f3(1,i1,i2,i3)=real((real(this%pi%f3(1,i1,i2,i3),p_double)-real(this%a%f3(1,i1,i2,i3),p_double))*idt,p_k_fld)
          enddo; enddo; enddo
        end select
        call report_vdf( rep, tmp, 1, g_space, grid, no_co, tstep, t )
        call tmp%cleanup()
      case ( p_aold )
        call report_vdf( rep, this%a, 1, g_space, grid, no_co, tstep, t )
      case ( p_anew )
        call report_vdf( rep, this%pi, 1, g_space, grid, no_co, tstep, t )
      case ( p_ae1, p_ae2, p_ae3 )
        call report_vdf( rep, this%pert%e, rep%quant - p_ae1 + 1, g_space, grid, no_co, tstep, t )
      case ( p_ab1, p_ab2, p_ab3 )
        call report_vdf( rep, this%pert%b, rep%quant - p_ab1 + 1, g_space, grid, no_co, tstep, t )
      case ( p_edotb )
        call report_vdf( rep, this%diag_edotb, 1, g_space, grid, no_co, tstep, t )
      case ( p_emsrcmag )
        call report_vdf( rep, this%diag_emsrcmag, 1, g_space, grid, no_co, tstep, t )
      end select
    endif
    rep => rep%next
  enddo
  call end_event( diag_axion_ev )
end subroutine report_axion_field_impl

! Compute and print integrated axion and perturbation energy diagnostics.
subroutine report_energy_axion_field_impl( this, no_co, tstep, emf, omega0_ref )
  class(t_axion_field), intent(in) :: this
  class(t_node_conf), intent(in) :: no_co
  type(t_time_step), intent(in) :: tstep
  class(t_emf), optional, intent(in) :: emf
  real(p_double), optional, intent(in) :: omega0_ref

  integer :: i1, i2, i3
  real(p_double) :: ekin, egrad, emass, vol, phi_c, phi_t, gradx, grady, gradz
  real(p_double) :: en1, bn1
  real(p_double) :: pert_e_local, pert_b_local, amp_local, edotb_avg, src_avg, de_avg
  real(p_double) :: pert_e_max, pert_b_max
  real(p_double) :: e0_max, b0_max, e0_perp_max, a0_est, gA0_est, omega0_used, ma_over_omega0, omega0_over_ma

  if ( .not. this%if_initialized ) return
  if ( mod( n(tstep), max(1, ndump(tstep)) ) /= 0 ) return

  ekin = 0.0_p_double
  egrad = 0.0_p_double
  emass = 0.0_p_double
  pert_e_local = 0.0_p_double
  pert_b_local = 0.0_p_double
  pert_e_max = 0.0_p_double
  pert_b_max = 0.0_p_double
  amp_local = 0.0_p_double
  e0_max = 0.0_p_double
  b0_max = 0.0_p_double
  e0_perp_max = 0.0_p_double
  a0_est = 0.0_p_double
  gA0_est = 0.0_p_double
  omega0_used = 0.0_p_double
  ma_over_omega0 = 0.0_p_double
  omega0_over_ma = 0.0_p_double

  vol = 1.0_p_double
  do i1 = 1, p_x_dim
    vol = vol * this%dx(i1)
  enddo

  select case ( p_x_dim )
  case(1)
    do i1 = 1, this%nx(1)
      phi_c = phi_center_1d(this,i1)
      phi_t = phi_t_1d(this,dt(tstep),i1)
      gradx = ( phi_center_1d(this,i1+1) - phi_center_1d(this,i1-1) ) / ( 2.0_p_double * this%dx(1) )
      ekin = ekin + 0.5_p_double * phi_t**2
      egrad = egrad + 0.5_p_double * gradx**2
      emass = emass + 0.5_p_double * this%ma**2 * phi_c**2
      amp_local = max( amp_local, abs(phi_c) )
      en1 = real(this%pert%e%f1(1,i1),p_double)**2 + real(this%pert%e%f1(2,i1),p_double)**2 + &
            real(this%pert%e%f1(3,i1),p_double)**2
      bn1 = real(this%pert%b%f1(1,i1),p_double)**2 + real(this%pert%b%f1(2,i1),p_double)**2 + &
            real(this%pert%b%f1(3,i1),p_double)**2
      pert_e_local = pert_e_local + 0.5_p_double * en1
      pert_b_local = pert_b_local + 0.5_p_double * bn1
      pert_e_max = max( pert_e_max, sqrt(en1) )
      pert_b_max = max( pert_b_max, sqrt(bn1) )
    enddo
  case(2)
    do i2 = 1, this%nx(2)
      do i1 = 1, this%nx(1)
        phi_c = phi_center_2d(this,i1,i2)
        phi_t = phi_t_2d(this,dt(tstep),i1,i2)
        gradx = ( phi_center_2d(this,i1+1,i2) - phi_center_2d(this,i1-1,i2) ) / ( 2.0_p_double * this%dx(1) )
        grady = ( phi_center_2d(this,i1,i2+1) - phi_center_2d(this,i1,i2-1) ) / ( 2.0_p_double * this%dx(2) )
        ekin = ekin + 0.5_p_double * phi_t**2
        egrad = egrad + 0.5_p_double * ( gradx**2 + grady**2 )
        emass = emass + 0.5_p_double * this%ma**2 * phi_c**2
        amp_local = max( amp_local, abs(phi_c) )
        en1 = real(this%pert%e%f2(1,i1,i2),p_double)**2 + real(this%pert%e%f2(2,i1,i2),p_double)**2 + &
              real(this%pert%e%f2(3,i1,i2),p_double)**2
        bn1 = real(this%pert%b%f2(1,i1,i2),p_double)**2 + real(this%pert%b%f2(2,i1,i2),p_double)**2 + &
              real(this%pert%b%f2(3,i1,i2),p_double)**2
        pert_e_local = pert_e_local + 0.5_p_double * en1
        pert_b_local = pert_b_local + 0.5_p_double * bn1
        pert_e_max = max( pert_e_max, sqrt(en1) )
        pert_b_max = max( pert_b_max, sqrt(bn1) )
      enddo
    enddo
  case(3)
    do i3 = 1, this%nx(3)
      do i2 = 1, this%nx(2)
        do i1 = 1, this%nx(1)
          phi_c = phi_center_3d(this,i1,i2,i3)
          phi_t = phi_t_3d(this,dt(tstep),i1,i2,i3)
          gradx = ( phi_center_3d(this,i1+1,i2,i3) - phi_center_3d(this,i1-1,i2,i3) ) / ( 2.0_p_double * this%dx(1) )
          grady = ( phi_center_3d(this,i1,i2+1,i3) - phi_center_3d(this,i1,i2-1,i3) ) / ( 2.0_p_double * this%dx(2) )
          gradz = ( phi_center_3d(this,i1,i2,i3+1) - phi_center_3d(this,i1,i2,i3-1) ) / ( 2.0_p_double * this%dx(3) )
          ekin = ekin + 0.5_p_double * phi_t**2
          egrad = egrad + 0.5_p_double * ( gradx**2 + grady**2 + gradz**2 )
          emass = emass + 0.5_p_double * this%ma**2 * phi_c**2
          amp_local = max( amp_local, abs(phi_c) )
          en1 = real(this%pert%e%f3(1,i1,i2,i3),p_double)**2 + real(this%pert%e%f3(2,i1,i2,i3),p_double)**2 + &
                real(this%pert%e%f3(3,i1,i2,i3),p_double)**2
          bn1 = real(this%pert%b%f3(1,i1,i2,i3),p_double)**2 + real(this%pert%b%f3(2,i1,i2,i3),p_double)**2 + &
                real(this%pert%b%f3(3,i1,i2,i3),p_double)**2
          pert_e_local = pert_e_local + 0.5_p_double * en1
          pert_b_local = pert_b_local + 0.5_p_double * bn1
          pert_e_max = max( pert_e_max, sqrt(en1) )
          pert_b_max = max( pert_b_max, sqrt(bn1) )
        enddo
      enddo
    enddo
  end select

  if ( present(emf) ) then
    call background_field_peaks_axion( this, emf, e0_max, b0_max, e0_perp_max )
  endif
  if ( present(omega0_ref) ) omega0_used = abs(omega0_ref)
  if ( omega0_used > tiny(1.0_p_double) ) then
    a0_est = e0_perp_max / omega0_used
    gA0_est = abs(this%g_agg) * a0_est
    ma_over_omega0 = abs(this%ma) / omega0_used
  endif
  if ( abs(this%ma) > tiny(1.0_p_double) ) then
    omega0_over_ma = omega0_used / abs(this%ma)
  endif

  edotb_avg = 0.0_p_double
  src_avg = 0.0_p_double
  de_avg = 0.0_p_double
  if ( this%src_count > 0 ) then
    edotb_avg = this%edotb_sum / real(this%src_count,p_double)
    src_avg = this%src_sum / real(this%src_count,p_double)
  endif
  if ( this%de_count > 0 ) de_avg = this%de_sum / real(this%de_count,p_double)

  if ( root(no_co) ) then
    write(*,'(A,1PE14.6,A,1PE14.6,A,1PE14.6,A,1PE14.6)') &
         'axion_field: ekin=', ekin*vol, ' egrad=', egrad*vol, &
         ' emass=', emass*vol, ' max|phi|=', amp_local
    write(*,'(A,1PE14.6,A,1PE14.6,A,1PE14.6,A,1PE14.6)') &
         '             pertE=', pert_e_local*vol, ' pertB=', pert_b_local*vol, &
         ' max|E1|=', pert_e_max, ' max|B1|=', pert_b_max
    write(*,'(A,1PE14.6,A,1PE14.6,A,1PE14.6,A,1PE14.6)') &
         '             max|E.B|=', this%edotb_max, ' avg(E.B)=', edotb_avg, &
         ' max|src|=', this%src_max, ' avg(src)=', src_avg
    write(*,'(A,1PE14.6,A,1PE14.6)') '             max|dE1|=', this%de_max, ' avg|dE1|=', de_avg
    if ( present(emf) ) then
      write(*,'(A,1PE14.6,A,1PE14.6,A,1PE14.6)') &
           '             g_agg(code)=', this%g_agg, ' ma(code)=', this%ma, ' max|E0|=', e0_max
      write(*,'(A,1PE14.6,A,1PE14.6)') '             max|B0|=', b0_max, ' max|E0_perp|=', e0_perp_max
      if ( omega0_used > tiny(1.0_p_double) ) then
        write(*,'(A,1PE14.6,A,1PE14.6,A,1PE14.6)') &
             '             omega0(code)=', omega0_used, ' A0_est=max|E0_perp|/omega0=', a0_est, ' |gA0|_est=', gA0_est
        write(*,'(A,1PE14.6,A,1PE14.6)') '             ma/omega0=', ma_over_omega0, ' omega0/ma=', omega0_over_ma
      else
        write(*,'(A)') '             omega0(code)=unavailable -> A0_est not printed'
      endif
    endif
  endif
end subroutine report_energy_axion_field_impl



! Estimate peak background field amplitudes for normalization diagnostics.
subroutine background_field_peaks_axion( this, emf, e0_max, b0_max, e0_perp_max )
  class(t_axion_field), intent(in) :: this
  class(t_emf), intent(in) :: emf
  real(p_double), intent(out) :: e0_max, b0_max, e0_perp_max

  integer :: i1, i2, i3
  real(p_double) :: ex, ey, ez, bx, by, bz

  e0_max = 0.0_p_double
  b0_max = 0.0_p_double
  e0_perp_max = 0.0_p_double

  select case ( p_x_dim )
  case (1)
    do i1 = 1, this%nx(1)
      ex = real(emf%e%f1(1,i1),p_double)
      ey = real(emf%e%f1(2,i1),p_double)
      ez = real(emf%e%f1(3,i1),p_double)
      bx = real(emf%b%f1(1,i1),p_double)
      by = real(emf%b%f1(2,i1),p_double)
      bz = real(emf%b%f1(3,i1),p_double)
      e0_max = max( e0_max, sqrt(ex*ex + ey*ey + ez*ez) )
      b0_max = max( b0_max, sqrt(bx*bx + by*by + bz*bz) )
      e0_perp_max = max( e0_perp_max, sqrt(ey*ey + ez*ez) )
    enddo
  case (2)
    do i2 = 1, this%nx(2)
      do i1 = 1, this%nx(1)
        ex = real(emf%e%f2(1,i1,i2),p_double)
        ey = real(emf%e%f2(2,i1,i2),p_double)
        ez = real(emf%e%f2(3,i1,i2),p_double)
        bx = real(emf%b%f2(1,i1,i2),p_double)
        by = real(emf%b%f2(2,i1,i2),p_double)
        bz = real(emf%b%f2(3,i1,i2),p_double)
        e0_max = max( e0_max, sqrt(ex*ex + ey*ey + ez*ez) )
        b0_max = max( b0_max, sqrt(bx*bx + by*by + bz*bz) )
        e0_perp_max = max( e0_perp_max, sqrt(ey*ey + ez*ez) )
      enddo
    enddo
  case (3)
    do i3 = 1, this%nx(3)
      do i2 = 1, this%nx(2)
        do i1 = 1, this%nx(1)
          ex = real(emf%e%f3(1,i1,i2,i3),p_double)
          ey = real(emf%e%f3(2,i1,i2,i3),p_double)
          ez = real(emf%e%f3(3,i1,i2,i3),p_double)
          bx = real(emf%b%f3(1,i1,i2,i3),p_double)
          by = real(emf%b%f3(2,i1,i2,i3),p_double)
          bz = real(emf%b%f3(3,i1,i2,i3),p_double)
          e0_max = max( e0_max, sqrt(ex*ex + ey*ey + ez*ez) )
          b0_max = max( b0_max, sqrt(bx*bx + by*by + bz*bz) )
          e0_perp_max = max( e0_perp_max, sqrt(ey*ey + ez*ez) )
        enddo
      enddo
    enddo
  end select
end subroutine background_field_peaks_axion

end module m_axion_field_diagnostics_impl

! Forward the diagnostic-quantity count query to the compiled implementation module.
function avail_report_quants_axion( this ) result(nq)
  use m_axion_field_diag_define, only : t_axion_diag
  use m_axion_field_diagnostics_impl, only : avail_report_quants_axion_impl
  implicit none
  class(t_axion_diag), intent(in) :: this
  integer :: nq
  nq = avail_report_quants_axion_impl( this )
end function avail_report_quants_axion

! Forward diagnostic-quantity initialization to the compiled implementation module.
subroutine init_report_quants_axion( this )
  use m_axion_field_diag_define, only : t_axion_diag
  use m_axion_field_diagnostics_impl, only : init_report_quants_axion_impl
  implicit none
  class(t_axion_diag), intent(inout) :: this
  call init_report_quants_axion_impl( this )
end subroutine init_report_quants_axion

! Forward axion diagnostic input parsing to the compiled implementation module.
subroutine read_input_diag_axion( this, input_file )
  use m_axion_field_diag_define, only : t_axion_diag
  use m_input_file, only : t_input_file
  use m_axion_field_diagnostics_impl, only : read_input_diag_axion_impl
  implicit none
  class(t_axion_diag), intent(inout) :: this
  class(t_input_file), intent(inout) :: input_file
  call read_input_diag_axion_impl( this, input_file )
end subroutine read_input_diag_axion

! Forward axion diagnostic initialization to the compiled implementation module.
subroutine init_diag_axion( this )
  use m_axion_field_diag_define, only : t_axion_diag
  use m_axion_field_diagnostics_impl, only : init_diag_axion_impl
  implicit none
  class(t_axion_diag), intent(inout) :: this
  call init_diag_axion_impl( this )
end subroutine init_diag_axion

! Forward axion diagnostic cleanup to the compiled implementation module.
subroutine cleanup_diag_axion( this )
  use m_axion_field_diag_define, only : t_axion_diag
  use m_axion_field_diagnostics_impl, only : cleanup_diag_axion_impl
  implicit none
  class(t_axion_diag), intent(inout) :: this
  call cleanup_diag_axion_impl( this )
end subroutine cleanup_diag_axion

! Forward axion-field reporting to the compiled implementation module.
subroutine report_axion_field( this, g_space, grid, no_co, tstep, t )
  use m_parameters, only : p_double
  use m_axion_field_define, only : t_axion_field
  use m_space, only : t_space
  use m_grid_define, only : t_grid
  use m_node_conf, only : t_node_conf
  use m_time_step, only : t_time_step
  use m_axion_field_diagnostics_impl, only : report_axion_field_impl
  implicit none
  class(t_axion_field), intent(inout) :: this
  type(t_space), intent(in) :: g_space
  class(t_grid), intent(in) :: grid
  class(t_node_conf), intent(in) :: no_co
  type(t_time_step), intent(in) :: tstep
  real(p_double), intent(in) :: t
  call report_axion_field_impl( this, g_space, grid, no_co, tstep, t )
end subroutine report_axion_field

! Forward axion energy reporting to the compiled implementation module.
subroutine report_energy_axion_field( this, no_co, tstep, emf, omega0_ref )
  use m_parameters, only : p_double
  use m_axion_field_define, only : t_axion_field
  use m_node_conf, only : t_node_conf
  use m_time_step, only : t_time_step
  use m_emf_define, only : t_emf
  use m_axion_field_diagnostics_impl, only : report_energy_axion_field_impl
  implicit none
  class(t_axion_field), intent(in) :: this
  class(t_node_conf), intent(in) :: no_co
  type(t_time_step), intent(in) :: tstep
  class(t_emf), optional, intent(in) :: emf
  real(p_double), optional, intent(in) :: omega0_ref
  call report_energy_axion_field_impl( this, no_co, tstep, emf, omega0_ref )
end subroutine report_energy_axion_field
