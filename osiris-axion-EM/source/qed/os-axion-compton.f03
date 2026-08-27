#include "os-preprocess.fpp"
#include "os-config.h"

module m_axion_compton

#include "memory/memory.h"

  use m_parameters
  use m_system
  use m_node_conf
  use m_random
  use m_vdf_define
  use m_vdf
  use m_species_define

implicit none

private

!---------------------------------------------------------------------------------------------------
! Thermal-bath Compton-like axion production (gamma + e -> e + a) information container
!
! This is designed to mirror t_axion_bremsstrahlung_info (ax_model=2) so that:
!  - rates use local ion density proxy rho_ion (Option A)
!  - emission uses Poisson criterion with unbiased variance-control weights
!  - inverse (absorption) can be implemented using detailed balance
!  - optional plasma feedback via per-cell depE/depP and evolving Te_cell
!---------------------------------------------------------------------------------------------------

type :: t_axion_compton_info

  logical :: if_compton
  integer :: Z_ion
  type( t_vdf ) :: rho_ion

  ! Legacy knob: if w_axion <= 0, emitted axion weight is |q_e|/proba_mult.
  real( p_double ) :: proba_mult

  ! Coupling and overall scaling
  real( p_double ) :: g_ae
  real( p_double ) :: cs_scale

  ! Thermal bath parameters
  real( p_double ) :: Te_keV
  logical :: assume_neutral
  real( p_double ) :: ne_factor

  ! Inverse process (absorption) toggle
  logical :: if_inverse
  logical :: use_bose

  ! Composition proxies (enter Raffelt emissivity)
  real( p_double ) :: X_H
  real( p_double ) :: A_ion

  ! Safety cap on Poisson mean (variance control). If <= 0, caller may set a default.
  real( p_double ) :: lambda_cap

  ! Variance-control knob: if > 0, emitted axions are created with this weight.
  real( p_double ) :: w_axion

  ! Plasma feedback / evolving temperature scaffolding (mirrors axion-brems)
  logical :: apply_feedback
  real( p_double ) :: fe_electrons
  logical :: evolve_Te
  real( p_double ) :: Te_floor_keV
  real( p_double ) :: Te_ceil_keV
  real( p_double ) :: Te_relax
  logical :: feedback_verbose

  type( t_vdf ) :: depE
  type( t_vdf ) :: depP
  type( t_vdf ) :: Te_cell

  real( p_double ) :: fb_dE_step
  real( p_double ), dimension(3) :: fb_dP_step
  integer :: fb_n_emit_step

  ! Pointer to the ion species
  class( t_species ), pointer :: ion
  integer, dimension(:), allocatable :: i_ion

contains

  procedure :: read_input       => read_input_axion_compton
  procedure :: ensure_feedback  => ensure_feedback_axion_compton
  procedure :: clear_feedback   => clear_feedback_axion_compton

end type t_axion_compton_info

public :: t_axion_compton_info

contains

  !-----------------------------------------------------------------------------------------
  ! Reads info for axion Compton-like production (thermal bath)
  !-----------------------------------------------------------------------------------------
  subroutine read_input_axion_compton( this, input_file, num_qed )

    use m_input_file, only : t_input_file, get_namelist

    implicit none

    class( t_axion_compton_info ), intent(inout) :: this
    type( t_input_file ), intent(inout) :: input_file
    integer, intent(in) :: num_qed

    integer :: ierr, Z_ion
    logical :: if_compton
    real(p_double) :: proba_mult, g_ae, cs_scale, Te_keV, ne_factor, w_axion, X_H, A_ion, lambda_cap
    logical :: assume_neutral, if_inverse, use_bose
    logical :: apply_feedback, evolve_Te, feedback_verbose
    real(p_double) :: fe_electrons, Te_floor_keV, Te_ceil_keV, Te_relax
    integer, dimension(:), allocatable :: i_ion

    character(len=256) :: msg
    namelist /nl_qed_axion_compton/ &
      if_compton, Z_ion, proba_mult, i_ion, &
      g_ae, cs_scale, &
      Te_keV, assume_neutral, ne_factor, &
      X_H, A_ion, lambda_cap, &
      if_inverse, use_bose, &
      w_axion, &
      apply_feedback, fe_electrons, evolve_Te, Te_floor_keV, Te_ceil_keV, Te_relax, feedback_verbose

    ! defaults
    if_compton = .false.
    Z_ion = -1
    proba_mult = 1.0_p_double

    g_ae     = 0.0_p_double
    cs_scale = 1.0_p_double

    Te_keV = 0.0_p_double
    assume_neutral = .true.
    ne_factor = 1.0_p_double

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

    ! NOTE on OSIRIS input decks:
    ! - In most builds, a deck block `qed_axion_compton{...}` is converted into a Fortran
    !   namelist group named `nl_qed_axion_compton`.
    ! - In some branches/merges, the block may be stored under the raw block name instead.
    !   (e.g. `qed_axion_compton`).
    !
    ! To match the behaviour of the working axion-brems and axion-primakoff readers while
    ! staying robust across branches, we try the canonical `nl_` name first and then fall
    ! back to the raw block name.
    call get_namelist( input_file, "nl_qed_axion_compton", ierr )
    if ( ierr /= 0 ) call get_namelist( input_file, "nl_qed_compton", ierr )
    if ( ierr /= 0 ) call get_namelist( input_file, "qed_axion_compton", ierr )
    if ( ierr /= 0 ) call get_namelist( input_file, "qed_compton", ierr )

    if ( ierr == 0 ) then
      read (input_file%nml_text, nml = nl_qed_axion_compton, iostat = ierr)
      if (ierr /= 0) then
        print *, ""
        print *, "   Error reading qed_axion_compton parameters "
        print *, "   aborting..."
        stop
      endif
    else
      SCR_ROOT("   - no particle axion compton specified")
    endif

    this%if_compton = if_compton
    this%Z_ion = Z_ion
    this%proba_mult = proba_mult
    this%i_ion = i_ion

    this%g_ae = g_ae
    this%cs_scale = cs_scale

    this%Te_keV = Te_keV
    this%assume_neutral = assume_neutral
    this%ne_factor = ne_factor

    this%if_inverse = if_inverse
    this%use_bose   = use_bose

    this%X_H = X_H
    this%A_ion = A_ion
    this%lambda_cap = lambda_cap

    this%w_axion = w_axion

    this%apply_feedback   = apply_feedback
    this%fe_electrons     = max(0.0_p_double, min(1.0_p_double, fe_electrons))
    this%evolve_Te        = evolve_Te
    this%Te_floor_keV     = Te_floor_keV
    this%Te_ceil_keV      = Te_ceil_keV
    this%Te_relax         = Te_relax
    this%feedback_verbose = feedback_verbose

    if (this%evolve_Te) this%apply_feedback = .true.

    this%fb_dE_step = 0.0_p_double
    this%fb_dP_step = 0.0_p_double
    this%fb_n_emit_step = 0

    if (this%if_compton .or. this%if_inverse) then
      SCR_ROOT(" ")
      SCR_ROOT(" [AXION_COMPTON DEBUG] values after namelist read_input_axion_compton:")

      write(msg,'(A,L1)')  "   if_compton      = ", this%if_compton      ; SCR_ROOT(trim(msg))
      write(msg,'(A,L1)')  "   if_inverse      = ", this%if_inverse      ; SCR_ROOT(trim(msg))
      write(msg,'(A,L1)')  "   use_bose        = ", this%use_bose        ; SCR_ROOT(trim(msg))
      write(msg,'(A,L1)')  "   assume_neutral  = ", this%assume_neutral  ; SCR_ROOT(trim(msg))

      write(msg,'(A,I0)')  "   Z_ion           = ", this%Z_ion           ; SCR_ROOT(trim(msg))
      write(msg,'(A,ES24.16)') "   ne_factor       = ", this%ne_factor     ; SCR_ROOT(trim(msg))

      write(msg,'(A,ES24.16)') "   Te_keV          = ", this%Te_keV        ; SCR_ROOT(trim(msg))
      write(msg,'(A,ES24.16)') "   g_ae            = ", this%g_ae          ; SCR_ROOT(trim(msg))
      write(msg,'(A,ES24.16)') "   cs_scale        = ", this%cs_scale      ; SCR_ROOT(trim(msg))

      write(msg,'(A,ES24.16)') "   w_axion         = ", this%w_axion       ; SCR_ROOT(trim(msg))
      write(msg,'(A,ES24.16)') "   proba_mult      = ", this%proba_mult    ; SCR_ROOT(trim(msg))
      write(msg,'(A,ES24.16)') "   lambda_cap      = ", this%lambda_cap    ; SCR_ROOT(trim(msg))

      write(msg,'(A,ES24.16)') "   X_H             = ", this%X_H           ; SCR_ROOT(trim(msg))
      write(msg,'(A,ES24.16)') "   A_ion           = ", this%A_ion         ; SCR_ROOT(trim(msg))

      write(msg,'(A,L1)')  "   apply_feedback  = ", this%apply_feedback  ; SCR_ROOT(trim(msg))
      write(msg,'(A,L1)')  "   evolve_Te       = ", this%evolve_Te       ; SCR_ROOT(trim(msg))
      write(msg,'(A,ES24.16)') "   fe_electrons    = ", this%fe_electrons  ; SCR_ROOT(trim(msg))
      write(msg,'(A,ES24.16)') "   Te_floor_keV    = ", this%Te_floor_keV  ; SCR_ROOT(trim(msg))
      write(msg,'(A,ES24.16)') "   Te_ceil_keV     = ", this%Te_ceil_keV   ; SCR_ROOT(trim(msg))
      write(msg,'(A,ES24.16)') "   Te_relax        = ", this%Te_relax      ; SCR_ROOT(trim(msg))

      if (allocated(this%i_ion)) then
        SCR_ROOT("   i_ion (first up to 8) = ")
        if (size(this%i_ion) >= 1) then
          write(msg,'(8(I0,1X))') this%i_ion(1:min(8,size(this%i_ion))) ; SCR_ROOT(trim(msg))
        end if
      else
        SCR_ROOT("   i_ion not allocated")
      end if

      SCR_ROOT(" ")
    end if


  end subroutine read_input_axion_compton

  !-----------------------------------------------------------------------------------------
  ! Ensure feedback VDFs are allocated and initialized (mirrors axion-brems).
  !-----------------------------------------------------------------------------------------
  subroutine ensure_feedback_axion_compton( this )

    implicit none

    class( t_axion_compton_info ), intent(inout) :: this
    integer :: lx, ux, ly, uy, lz, uz

    if (.not. this%apply_feedback) return

    select case ( p_x_dim )
    case ( 1 )
      if (.not. associated(this%rho_ion%f1)) return
      if (.not. associated(this%depE%f1)) then
        lx = lbound(this%rho_ion%f1, 2); ux = ubound(this%rho_ion%f1, 2)
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
        lx = lbound(this%rho_ion%f2, 2); ux = ubound(this%rho_ion%f2, 2)
        ly = lbound(this%rho_ion%f2, 3); uy = ubound(this%rho_ion%f2, 3)
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
        lx = lbound(this%rho_ion%f3, 2); ux = ubound(this%rho_ion%f3, 2)
        ly = lbound(this%rho_ion%f3, 3); uy = ubound(this%rho_ion%f3, 3)
        lz = lbound(this%rho_ion%f3, 4); uz = ubound(this%rho_ion%f3, 4)
        allocate( this%depE%f3(1, lx:ux, ly:uy, lz:uz) )
        allocate( this%depP%f3(3, lx:ux, ly:uy, lz:uz) )
        allocate( this%Te_cell%f3(1, lx:ux, ly:uy, lz:uz) )
        this%depE%f3 = 0.0_p_double
        this%depP%f3 = 0.0_p_double
        this%Te_cell%f3 = this%Te_keV
      end if
    end select

  end subroutine ensure_feedback_axion_compton

  !-----------------------------------------------------------------------------------------
  ! Clear per-step deposits and counters.
  !-----------------------------------------------------------------------------------------
  subroutine clear_feedback_axion_compton( this )

    implicit none

    class( t_axion_compton_info ), intent(inout) :: this

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

  end subroutine clear_feedback_axion_compton

end module m_axion_compton
