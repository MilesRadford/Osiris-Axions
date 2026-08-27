#include "os-preprocess.fpp"
#include "os-config.h"

module m_axion_bremsstrahlung

#include "memory/memory.h"
  !use m_random_class

  use m_parameters
  use m_system
  use m_node_conf
  use m_random
  use m_vdf_define
  use m_vdf
  use m_species_define

implicit none

private

! -----
! CLASS DEFINITIONS
! -----

type :: t_axion_bremsstrahlung_info

  logical :: if_bremsstrahlung
  integer :: Z_ion
  type( t_vdf ) :: rho_ion

  ! Multiplicative factor applied on the (macroparticle) emission mean.
  ! Legacy knob: if w_axion <= 0, emitted axion weight is |q_e|/proba_mult.
  real( p_double ) :: proba_mult

  ! Axion model selector
  !  0 : photon-clone (uses Martinez photon Coulomb brems as template)
  !  1 : photon-clone with simple g_ae^2 scaling (placeholder)
  !  2 : thermal plasma axion bremsstrahlung (solar-style, Option A)
  integer :: ax_model

  ! Physics knobs (used mainly for ax_model = 1/2)
  real( p_double ) :: g_ae
  real( p_double ) :: cs_scale

  ! Thermal plasma parameters (ax_model = 2)
  real( p_double ) :: Te_keV
  logical :: include_ee
  logical :: assume_neutral
  real( p_double ) :: ne_factor
  logical :: use_screening

  ! Inverse process (absorption) toggle (used mainly for ax_model = 2)
  ! If enabled, axions may be absorbed via inverse bremsstrahlung using
  ! detailed balance with the same local emissivity model used for emission.
  logical :: if_inverse

  ! If true, use Bose-Einstein equilibrium occupation in detailed-balance
  ! absorption coefficient. If false, uses the Maxwell-Boltzmann limit.
  logical :: use_bose

  ! Composition/ion-mass assumptions (ax_model = 2)
  real( p_double ) :: X_H
  real( p_double ) :: A_ion

  ! Safety cap on Poisson mean (variance control). If <= 0, caller may set a default.
  real( p_double ) :: lambda_cap

  ! Variance-control knob: if > 0, emitted axions are created with this weight
  ! (i.e. each macroparticle represents w_axion real axions).
  real( p_double ) :: w_axion

  ! -------------------------------------------------------------------
  ! Plasma feedback / evolving temperature (Layer 0/2/3 scaffolding)
  !
  ! If apply_feedback is enabled, emission/absorption deposit energy and
  ! momentum into per-cell accumulators (depE/depP). A separate operator
  ! (implemented elsewhere) can then apply these deposits to the plasma.
  !
  ! If evolve_Te is enabled, the local temperature field Te_cell (keV) is
  ! used instead of the constant Te_keV when computing rates/spectra.
  ! -------------------------------------------------------------------
  logical :: apply_feedback
  real( p_double ) :: fe_electrons
  logical :: evolve_Te
  real( p_double ) :: Te_floor_keV
  real( p_double ) :: Te_ceil_keV
  real( p_double ) :: Te_relax
  logical :: feedback_verbose

  ! Per-cell deposits (units: depE in mec^2, depP in mec)
  type( t_vdf ) :: depE
  type( t_vdf ) :: depP

  ! Local electron temperature field (keV)
  type( t_vdf ) :: Te_cell

  ! Simple rank-local counters (useful while validating)
  real( p_double ) :: fb_dE_step
  real( p_double ), dimension(3) :: fb_dP_step
  integer :: fb_n_emit_step

  ! Kept for compatibility with the photon-brems template code path
  integer :: cs_model
  logical :: energy_damp
  real( p_double ) :: p_emit_cutoff

  ! Pointer to the ion species for bremsstrahlung
  class( t_species ), pointer :: ion

  integer, dimension(:), allocatable :: i_ion

contains

  procedure :: read_input => read_input_axion_bremsstrahlung
  procedure :: ensure_feedback => ensure_feedback_axion_bremsstrahlung
  procedure :: clear_feedback  => clear_feedback_axion_bremsstrahlung

end type

public :: t_axion_bremsstrahlung_info

contains

  !-----------------------------------------------------------------------------------------
  ! Reads info for axion bremsstrahlung
  !-----------------------------------------------------------------------------------------
  subroutine read_input_axion_bremsstrahlung( this, input_file, num_qed )

    use m_input_file, only : t_input_file, get_namelist

    implicit none

    class( t_axion_bremsstrahlung_info ), intent(inout) :: this
    type( t_input_file ), intent(inout) :: input_file
    integer, intent(in) :: num_qed

    ! local variables (namelist)
    logical :: if_bremsstrahlung, energy_damp
    integer :: ierr, Z_ion, cs_model, ax_model
    real(p_double) :: proba_mult, g_ae, cs_scale, Te_keV, ne_factor, w_axion, X_H, A_ion, lambda_cap
    logical :: include_ee, assume_neutral, use_screening
    logical :: if_inverse, use_bose
    logical :: apply_feedback, evolve_Te, feedback_verbose
    real(p_double) :: fe_electrons, Te_floor_keV, Te_ceil_keV, Te_relax
    integer, dimension(:), allocatable :: i_ion

    namelist /nl_qed_axion_bremsstrahlung/ &
      if_bremsstrahlung, Z_ion, proba_mult, cs_model, energy_damp, i_ion, &
      ax_model, g_ae, cs_scale, &
      Te_keV, include_ee, assume_neutral, ne_factor, use_screening, &
      X_H, A_ion, lambda_cap, &
      if_inverse, use_bose, &
      w_axion, &
      apply_feedback, fe_electrons, evolve_Te, Te_floor_keV, Te_ceil_keV, Te_relax, feedback_verbose

    ! defaults
    if_bremsstrahlung = .false.
    Z_ion = -1
    proba_mult = 1.0_p_double
    cs_model = -1
    energy_damp = .true.

    ax_model = 0
    g_ae     = 0.0_p_double
    cs_scale = 1.0_p_double

    Te_keV = 0.0_p_double
    include_ee = .false.
    assume_neutral = .true.
    ne_factor = 1.0_p_double
    use_screening = .true.

    if_inverse = .false.
    use_bose   = .true.

    X_H = 1.0_p_double
    A_ion = 1.0_p_double
    lambda_cap = 0.0_p_double

    w_axion = 0.0_p_double

    apply_feedback  = .false.
    fe_electrons    = 1.0_p_double
    evolve_Te        = .false.
    Te_floor_keV     = 1.0d-6
    Te_ceil_keV      = 1.0d6
    Te_relax         = 1.0_p_double
    feedback_verbose = .false.

    allocate(i_ion(num_qed))
    i_ion = 0

    call get_namelist( input_file, "nl_qed_axion_bremsstrahlung", ierr )

    if ( ierr == 0 ) then
      read (input_file%nml_text, nml = nl_qed_axion_bremsstrahlung, iostat = ierr)
      if (ierr /= 0) then
        print *, ""
        print *, "   Error reading qed_axion_bremsstrahlung parameters "
        print *, "   aborting..."
        stop
      endif
    else
      SCR_ROOT("   - no particle axion bremsstrahlung specified")
    endif

    this%if_bremsstrahlung = if_bremsstrahlung
    this%Z_ion = Z_ion
    this%proba_mult = proba_mult
    this%cs_model = cs_model
    this%energy_damp = energy_damp
    this%i_ion = i_ion

    this%ax_model = ax_model
    this%g_ae = g_ae
    this%cs_scale = cs_scale

    this%Te_keV = Te_keV
    this%include_ee = include_ee
    this%assume_neutral = assume_neutral
    this%ne_factor = ne_factor
    this%use_screening = use_screening

    this%if_inverse = if_inverse
    this%use_bose   = use_bose

    this%X_H = X_H
    this%A_ion = A_ion
    this%lambda_cap = lambda_cap

    this%w_axion = w_axion

    ! feedback / evolving Te
    this%apply_feedback  = apply_feedback
    this%fe_electrons    = max(0.0_p_double, min(1.0_p_double, fe_electrons))
    this%evolve_Te       = evolve_Te
    this%Te_floor_keV    = Te_floor_keV
    this%Te_ceil_keV     = Te_ceil_keV
    this%Te_relax        = Te_relax
    this%feedback_verbose = feedback_verbose

    ! If Te is evolving, we will need the feedback storage.
    if (this%evolve_Te) this%apply_feedback = .true.

    this%fb_dE_step = 0.0_p_double
    this%fb_dP_step = 0.0_p_double
    this%fb_n_emit_step = 0

  end subroutine read_input_axion_bremsstrahlung

  !-----------------------------------------------------------------------------------------
  ! Ensure feedback VDFs are allocated and initialized.
  ! We allocate depE (1 comp), depP (3 comps), and Te_cell (1 comp) with the same
  ! local mesh extents as rho_ion.
  !-----------------------------------------------------------------------------------------
  subroutine ensure_feedback_axion_bremsstrahlung( this )

    implicit none

    class( t_axion_bremsstrahlung_info ), intent(inout) :: this

    integer :: nx, ny, nz
    integer :: lx, ux, ly, uy, lz, uz

    if (.not. this%apply_feedback) return

    select case ( p_x_dim )
    case ( 1 )
      if (.not. associated(this%rho_ion%f1)) return
      if (.not. associated(this%depE%f1)) then
lx = lbound(this%rho_ion%f1, 2)
ux = ubound(this%rho_ion%f1, 2)
allocate( this%depE%f1(1, lx:ux) )
allocate( this%depP%f1(3, lx:ux) )
allocate( this%Te_cell%f1(1, lx:ux) )
this%depE%f1 = 0.0_p_double
this%depP%f1 = 0.0_p_double
this%Te_cell%f1 = this%Te_keV
      end if

    case ( 2 )
      if (.not. associated(this%rho_ion%f2)) return
      if (.not. associated(this%depE%f2)) then
lx = lbound(this%rho_ion%f2, 2)
ux = ubound(this%rho_ion%f2, 2)
ly = lbound(this%rho_ion%f2, 3)
uy = ubound(this%rho_ion%f2, 3)
allocate( this%depE%f2(1, lx:ux, ly:uy) )
allocate( this%depP%f2(3, lx:ux, ly:uy) )
allocate( this%Te_cell%f2(1, lx:ux, ly:uy) )
this%depE%f2 = 0.0_p_double
this%depP%f2 = 0.0_p_double
this%Te_cell%f2 = this%Te_keV
      end if

    case ( 3 )
      if (.not. associated(this%rho_ion%f3)) return
      if (.not. associated(this%depE%f3)) then
lx = lbound(this%rho_ion%f3, 2)
ux = ubound(this%rho_ion%f3, 2)
ly = lbound(this%rho_ion%f3, 3)
uy = ubound(this%rho_ion%f3, 3)
lz = lbound(this%rho_ion%f3, 4)
uz = ubound(this%rho_ion%f3, 4)
allocate( this%depE%f3(1, lx:ux, ly:uy, lz:uz) )
allocate( this%depP%f3(3, lx:ux, ly:uy, lz:uz) )
allocate( this%Te_cell%f3(1, lx:ux, ly:uy, lz:uz) )
this%depE%f3 = 0.0_p_double
this%depP%f3 = 0.0_p_double
this%Te_cell%f3 = this%Te_keV
      end if
    end select

  end subroutine ensure_feedback_axion_bremsstrahlung

  !-----------------------------------------------------------------------------------------
  ! Clear per-step deposits and counters.
  !-----------------------------------------------------------------------------------------
  subroutine clear_feedback_axion_bremsstrahlung( this )

    implicit none

    class( t_axion_bremsstrahlung_info ), intent(inout) :: this

    if (.not. this%apply_feedback) return

    select case ( p_x_dim )
    case ( 1 )
      if (associated(this%depE%f1)) this%depE%f1 = 0.0_p_double
      if (associated(this%depP%f1)) this%depP%f1 = 0.0_p_double
    case ( 2 )
      if (associated(this%depE%f2)) this%depE%f2 = 0.0_p_double
      if (associated(this%depP%f2)) this%depP%f2 = 0.0_p_double
    case ( 3 )
      if (associated(this%depE%f3)) this%depE%f3 = 0.0_p_double
      if (associated(this%depP%f3)) this%depP%f3 = 0.0_p_double
    end select

    this%fb_dE_step = 0.0_p_double
    this%fb_dP_step = 0.0_p_double
    this%fb_n_emit_step = 0

  end subroutine clear_feedback_axion_bremsstrahlung

end module m_axion_bremsstrahlung
