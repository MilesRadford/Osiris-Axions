module m_axion_field_define

#include "os-config.h"
#include "os-preprocess.fpp"

  use m_system
  use m_parameters
  use m_input_file, only : t_input_file
  use m_vdf_define
  use m_vdf
  use m_vdf_comm
  use m_space
  use m_grid_define
  use m_node_conf
  use m_time_step
  use m_emf_define, only : t_emf
  use m_axion_field_diag_define, only : t_axion_diag

  implicit none
  private

  type :: t_axion_perturb_em
    type(t_vdf) :: e
    type(t_vdf) :: b
    logical :: if_initialized = .false.
  end type t_axion_perturb_em

  type, public :: t_axion_field
    type(t_vdf) :: a
    type(t_vdf) :: pi

    real(p_double) :: ma = 0.0_p_double
    real(p_double) :: g_agg = 0.0_p_double

    character(len=32) :: init_type = 'zero'
    real(p_double)    :: init_ampl = 0.0_p_double
    integer           :: init_kx = 1
    integer           :: init_ky = 0
    integer           :: init_kz = 0
    integer           :: init_kx_min = 1
    integer           :: init_kx_max = 1
    real(p_double)    :: init_phase_step = 0.0_p_double
    logical           :: init_band_norm = .false.
    real(p_double)    :: init_lambda = 0.0_p_double
    real(p_double)    :: init_phase = 0.0_p_double
    real(p_double)    :: init_x0 = 0.0_p_double
    real(p_double)    :: init_y0 = 0.0_p_double
    real(p_double)    :: init_z0 = 0.0_p_double
    real(p_double)    :: init_dir(3) = (/ 1.0_p_double, 0.0_p_double, 0.0_p_double /)
    real(p_double)    :: init_w_long = 0.0_p_double
    real(p_double)    :: init_w_trans1 = 0.0_p_double
    real(p_double)    :: init_w_trans2 = 0.0_p_double
    real(p_double)    :: init_vphase = 1.0_p_double
    character(len=32) :: bc_type(2,p_max_dim) = reshape( (/ &
      'periodic                        ', 'periodic                        ', &
      'periodic                        ', 'periodic                        ', &
      'periodic                        ', 'periodic                        ' /), (/2,p_max_dim/) )

    logical :: if_initialized = .false.
    logical :: first_advance = .true.
    integer :: x_dim = -1
    integer, dimension(p_max_dim) :: nx = 0
    integer, dimension(2,p_max_dim) :: gc_num = 0
    real(p_double), dimension(p_max_dim) :: dx = 0.0_p_double
    real(p_double), dimension(p_max_dim) :: xmin_loc = 0.0_p_double
    real(p_double), dimension(p_max_dim) :: xmax_loc = 0.0_p_double
    real(p_double) :: edotb_max = 0.0_p_double
    real(p_double) :: edotb_sum = 0.0_p_double
    real(p_double) :: src_max   = 0.0_p_double
    real(p_double) :: src_sum   = 0.0_p_double
    integer        :: src_count = 0
    real(p_double) :: de_max    = 0.0_p_double
    real(p_double) :: de_sum    = 0.0_p_double
    integer        :: de_count  = 0
    type(t_vdf)     :: diag_edotb
    type(t_vdf)     :: diag_emsrcmag
    type(t_axion_perturb_em) :: pert
    type(t_axion_diag)       :: diag
  contains
    procedure :: read_input         => read_input_axion_field
    procedure :: init               => init_axion_field
    procedure :: cleanup            => cleanup_axion_field
    procedure :: move_window        => move_window_axion_field
    procedure :: update_boundary    => update_boundary_axion_field
    procedure :: advance_free       => advance_free_axion_field
    procedure :: advance_driven     => advance_driven_axion_field
    procedure :: apply_backreaction => apply_backreaction_axion_field
    procedure :: reset_perturbation => reset_perturbation_axion_field
    procedure :: advance_perturb_b  => advance_perturb_b_axion_field
    procedure :: report_energy      => report_energy_axion_field
    procedure :: report             => report_axion_field
    procedure :: list_algorithm     => list_algorithm_axion_field
  end type t_axion_field

  interface
    ! Forward the axion-field input read to the compiled implementation module.
    subroutine read_input_axion_field( this, input_file )
      import t_axion_field, t_input_file
      class(t_axion_field), intent(inout) :: this
      class(t_input_file),  intent(inout) :: input_file
    end subroutine read_input_axion_field
  end interface

  interface
    ! Forward axion-field initialization to the compiled implementation module.
    subroutine init_axion_field( this, g_space, grid, no_co, gc_num_in, dt_init )
      import t_axion_field, t_space, t_grid, t_node_conf, p_x_dim, p_double
      class(t_axion_field), intent(inout) :: this
      type(t_space),        intent(in)    :: g_space
      class(t_grid),        intent(in)    :: grid
      class(t_node_conf),   intent(in)    :: no_co
      integer, dimension(2,p_x_dim), intent(in) :: gc_num_in
      real(p_double), intent(in) :: dt_init
    end subroutine init_axion_field
  end interface

  interface
    ! Forward axion-field cleanup to the compiled implementation module.
    subroutine cleanup_axion_field( this )
      import t_axion_field
      class(t_axion_field), intent(inout) :: this
    end subroutine cleanup_axion_field
  end interface

  interface
    ! Forward the axion-field moving-window update to the compiled implementation module.
    subroutine move_window_axion_field( this, g_space )
      import t_axion_field, t_space
      class(t_axion_field), intent(inout) :: this
      type(t_space),        intent(in)    :: g_space
    end subroutine move_window_axion_field
  end interface

  interface
    ! Forward the axion-field boundary update to the compiled implementation module.
    subroutine update_boundary_axion_field( this, no_co, send_msg, recv_msg )
      import t_axion_field, t_node_conf, t_vdf_msg
      class(t_axion_field), intent(inout) :: this
      class(t_node_conf),   intent(in)    :: no_co
      type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg
    end subroutine update_boundary_axion_field
  end interface

  interface
    ! Forward the free axion-field update to the compiled implementation module.
    subroutine advance_free_axion_field( this, dt, no_co, send_msg, recv_msg )
      import t_axion_field, t_node_conf, t_vdf_msg, p_double
      class(t_axion_field), intent(inout) :: this
      real(p_double),       intent(in)    :: dt
      class(t_node_conf),   intent(in)    :: no_co
      type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg
    end subroutine advance_free_axion_field
  end interface

  interface
    ! Forward the driven axion-field update to the compiled implementation module.
    subroutine advance_driven_axion_field( this, emf, dt, no_co, send_msg, recv_msg )
      import t_axion_field, t_emf, t_node_conf, t_vdf_msg, p_double
      class(t_axion_field), intent(inout) :: this
      class(t_emf),         intent(in)    :: emf
      real(p_double),       intent(in)    :: dt
      class(t_node_conf),   intent(in)    :: no_co
      type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg
    end subroutine advance_driven_axion_field
  end interface

  interface
    ! Forward the axion backreaction update to the compiled implementation module.
    subroutine apply_backreaction_axion_field( this, emf, dt_half, dt_full, jay1 )
      import t_axion_field, t_emf, t_vdf, p_double
      class(t_axion_field), intent(inout) :: this
      class(t_emf),         intent(in)    :: emf
      real(p_double),       intent(in)    :: dt_half, dt_full
      class(t_vdf), optional, intent(in)  :: jay1
    end subroutine apply_backreaction_axion_field
  end interface

  interface
    ! Forward the perturbation reset to the compiled implementation module.
    subroutine reset_perturbation_axion_field( this )
      import t_axion_field
      class(t_axion_field), intent(inout) :: this
    end subroutine reset_perturbation_axion_field
  end interface

  interface
    ! Forward the perturbative magnetic-field half-step to the compiled implementation module.
    subroutine advance_perturb_b_axion_field( this, dt_half )
      import t_axion_field, p_double
      class(t_axion_field), intent(inout) :: this
      real(p_double),       intent(in)    :: dt_half
    end subroutine advance_perturb_b_axion_field
  end interface

  interface
    ! Forward axion energy reporting to the compiled implementation module.
    subroutine report_energy_axion_field( this, no_co, tstep, emf, omega0_ref )
      import t_axion_field, t_node_conf, t_time_step, t_emf, p_double
      class(t_axion_field), intent(in) :: this
      class(t_node_conf),   intent(in) :: no_co
      type(t_time_step),    intent(in) :: tstep
      class(t_emf), optional, intent(in) :: emf
      real(p_double), optional, intent(in) :: omega0_ref
    end subroutine report_energy_axion_field
  end interface

  interface
    ! Forward axion-field reporting to the compiled implementation module.
    subroutine report_axion_field( this, g_space, grid, no_co, tstep, t )
      import t_axion_field, t_space, t_grid, t_node_conf, t_time_step, p_double
      class(t_axion_field), intent(inout) :: this
      type(t_space), intent(in) :: g_space
      class(t_grid), intent(in) :: grid
      class(t_node_conf), intent(in) :: no_co
      type(t_time_step), intent(in) :: tstep
      real(p_double), intent(in) :: t
    end subroutine report_axion_field
  end interface

  interface
    ! Forward the axion-field algorithm listing to the compiled implementation module.
    subroutine list_algorithm_axion_field( this )
      import t_axion_field
      class(t_axion_field), intent(in) :: this
    end subroutine list_algorithm_axion_field
  end interface

end module m_axion_field_define
