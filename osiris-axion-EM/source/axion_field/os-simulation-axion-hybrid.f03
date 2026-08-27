module m_simulation_axion_hybrid

#include "os-config.h"
#include "os-preprocess.fpp"

  use m_parameters, only : p_x_dim
  use m_system
  use m_simulation_qed
  use m_axion_field_define
  use m_current_define, only : t_current
  use m_axion_fps_context, only : set_axion_fps_context, clear_axion_fps_context
  use m_input_file
  use m_zpulse, only : t_zpulse_list
  use m_zpulse_std, only : t_zpulse

  implicit none
  private

  type, extends( t_simulation_qed ) :: t_simulation_axion_hybrid
    class(t_axion_field), pointer :: axf => null()
    class(t_current), pointer     :: jay1 => null()
  contains
    procedure :: allocate_objs    => allocate_objs_axion_hybrid
    procedure :: read_input       => read_input_axion_hybrid
    procedure :: init             => init_axion_hybrid
    procedure :: cleanup          => cleanup_axion_hybrid
    procedure :: iter             => iter_axion_hybrid
    procedure :: list_algorithm   => list_algorithm_axion_hybrid
  end type t_simulation_axion_hybrid

  public :: t_simulation_axion_hybrid

contains

! Allocate the extra objects needed by the hybrid axion simulation.
subroutine allocate_objs_axion_hybrid( sim )
  class(t_simulation_axion_hybrid), intent(inout) :: sim

  if ( .not. associated(sim%axf) ) allocate( sim%axf )
  if ( .not. associated(sim%jay1) ) allocate( sim%jay1 )
  call sim%t_simulation_qed%allocate_objs()
end subroutine allocate_objs_axion_hybrid

! Read the inherited QED input and then the coherent axion-field input.
subroutine read_input_axion_hybrid( sim, input_file )
  class(t_simulation_axion_hybrid), intent(inout) :: sim
  class(t_input_file), intent(inout) :: input_file

  call sim%t_simulation_qed%read_input( input_file )
  if ( associated(sim%axf) ) call sim%axf%read_input( input_file )
end subroutine read_input_axion_hybrid

! Initialize the inherited simulation state and then construct the coherent axion field.
subroutine init_axion_hybrid( sim, grid_gc_min )
  use m_particles, only : get_min_gc, interpolation
  use m_time_step, only : dt
  use m_space, only : if_move
  use m_species_define, only : t_species
  use m_restart, only : t_restart_handle

  class(t_simulation_axion_hybrid), intent(inout) :: sim
  integer, dimension(2,p_x_dim), intent(in) :: grid_gc_min
  integer, dimension(2,p_x_dim) :: ax_gc
  type(t_restart_handle) :: restart_handle_dummy
  class(t_species), pointer :: species

  call sim%t_simulation_qed%init( grid_gc_min )

  if ( associated(sim%jay1) ) then
    call sim%jay1%allocate_objs()
    sim%jay1%curr_smooth = sim%jay%curr_smooth
    call sim%jay1%init( sim%grid, sim%jay%get_min_gc(), sim%jay%dx(), &
                        sim%no_co, if_move(sim%g_space), sim%emf%bc_type(), &
                        interpolation(sim%part), .false., restart_handle_dummy, sim%options )
    call sim%jay1%pf(1)%zero()
  endif

  species => sim%part%species
  do
    if ( .not. associated(species) ) exit
    call species%enable_fps_state()
    species => species%next
  end do

  if ( associated(sim%axf) ) then
    ax_gc = get_min_gc( sim%part )
    call sim%axf%init( sim%g_space, sim%grid, sim%no_co, ax_gc, dt(sim%tstep) )
    call sim%axf%update_boundary( sim%no_co, sim%bnd%send_vdf, sim%bnd%recv_vdf )
  endif
end subroutine init_axion_hybrid

! Clean up the extra hybrid objects before delegating to the base simulation cleanup.
subroutine cleanup_axion_hybrid( sim )
  class(t_simulation_axion_hybrid), intent(inout) :: sim

  if ( associated(sim%axf) ) then
    call sim%axf%cleanup()
    deallocate( sim%axf )
    nullify( sim%axf )
  endif

  if ( associated(sim%jay1) ) then
    call sim%jay1%cleanup()
    deallocate( sim%jay1 )
    nullify( sim%jay1 )
  endif

  call sim%t_simulation_qed%cleanup()
end subroutine cleanup_axion_hybrid

! Advance one timestep of the hybrid EM-particle-coherent-axion algorithm.
subroutine iter_axion_hybrid( sim )

  use m_emf
  use m_particles
  use m_logprof
  use m_particles_charge
  use m_space
  use m_vdf_define
  use m_time_step, only : dt, n, ndump, advance
  use m_time, only : t, update
  use m_current_define
  use m_particles_qed
  use m_zpulse
  use m_zpulse_std, only : t_zpulse
  use m_antenna_array
  use m_simulation, only : loopev
  use m_node_conf, only : root, n_threads
  use m_parameters, only : p_lower, p_upper, p_cylindrical_b, p_double

  implicit none

  class ( t_simulation_axion_hybrid ), intent(inout) :: sim
  real(p_double) :: pump_omega0, pump_a0

  call strongest_zpulse_params( sim%zpulse_list, pump_omega0, pump_a0 )

  call begin_event(loopev)

  call advance( sim%tstep, 1 )
  call update( sim%time, dt(sim%tstep), n(sim%tstep) )

  if ( root(sim%no_co) ) then
     print '(A,I10,A,F0.3)', 'n = ', n(sim%tstep), ', t = ', t(sim%time)
  endif

  call sim%jay%report(sim%g_space, sim%grid, sim%no_co, sim%tstep, t(sim%time) )

  call sim%zpulse_list%launch( sim%emf, sim%g_space, &
               sim%grid%my_nx(p_lower, :), sim%grid%g_nx, &
               t(sim%time), dt(sim%tstep), sim%no_co )

  call move_boundary( sim%g_space, sim % emf % dx(), dt(sim%tstep) )

  call sim%emf % move_window( sim%g_space, sim%grid%my_nx( p_lower, : ), &
                    sim%no_co%on_edge( 1, p_upper ) )

  call sim%part % move_window( sim%g_space, sim%grid, &
                     sim%jay, sim%no_co, sim%bnd )

  if ( associated(sim%axf) ) then
    call sim%axf%move_window( sim%g_space )
  endif

  if (antenna_on(sim%antenna_array)) then
     call sim % antenna_array % launch( sim%emf, dt(sim%tstep), t(sim%time), &
                                        sim%g_space, sim%grid%my_nx( p_lower, :) )
  end if

  call sim%emf%update_boundary( sim%g_space, sim%no_co, sim%bnd%send_vdf, sim%bnd%recv_vdf )

  call sim%emf%update_particle_fld( n(sim%tstep), t(sim%time), dt(sim%tstep), sim%g_space, &
                            sim%grid%my_nx( p_lower, : ))


if ( associated(sim%axf) ) then
  call sim%axf%update_boundary( sim%no_co, sim%bnd%send_vdf, sim%bnd%recv_vdf )

  ! Paper-style staggered coherent update, first half-step.
  call sim%axf%apply_backreaction( sim%emf, 0.5_p_double * dt(sim%tstep), dt(sim%tstep) )
  call sim%axf%advance_perturb_b( 0.5_p_double * dt(sim%tstep) )
  call sim%axf%advance_driven( sim%emf, dt(sim%tstep), sim%no_co, sim%bnd%send_vdf, sim%bnd%recv_vdf )
endif

if ( sim % emf % if_charge_cons( sim%tstep ) ) then
    call update_charge( sim%part, sim%g_space, sim%grid, sim%no_co, sim%bnd%send_vdf, sim%bnd%recv_vdf )
    ERROR( "Charge conservation diagnostic is temporarily disabled" )
    call abort_program()
  endif

  call sim % emf % report( sim%g_space, sim%grid,  sim%no_co, sim%tstep, t(sim%time), sim%bnd%send_vdf, sim%bnd%recv_vdf )
  if ( associated(sim%axf) ) then
    call sim%axf%report( sim%g_space, sim%grid, sim%no_co, sim%tstep, t(sim%time) )
  endif
  call sim % part % report( sim%emf, sim%g_space, sim%grid, sim%no_co, &
                    sim%tstep, t(sim%time), sim%bnd%send_vdf, sim%bnd%recv_vdf )

  if ( associated(sim%axf) ) then
    call set_axion_fps_context( sim%axf%pert%e, sim%axf%pert%b, sim%emf%part_grid_center )
    call sim % part % advance_deposit( sim%emf, sim%jay, sim%tstep, t(sim%time), &
                                       sim%no_co, sim%options, sim%jay1 )
    call clear_axion_fps_context()
  else
    call sim % part % advance_deposit( sim%emf, sim%jay, sim%tstep, t(sim%time), &
                                       sim%no_co, sim%options )
  endif

  call sim % part % update_boundary( sim%jay, sim%g_space, sim%no_co, dt(sim%tstep), sim%bnd )

  call sim % emf % report_energy( sim%no_co, sim%tstep, t(sim%time))
  call sim % part % report_energy( sim%no_co, sim%tstep, t(sim%time))

  call cathode(sim%part, sim%jay%pf(1), t(sim%time), dt(sim%tstep), &
               sim%no_co, sim%grid%coordinates)

  call ionize_neutral(sim%part, sim%emf, dt(sim%tstep), sim%grid%coordinates)

  call sort_collide( sim%part, n(sim%tstep), t(sim%time), n_threads(sim%no_co) )

  if ( sim%grid%coordinates == p_cylindrical_b ) then
    call sim%jay%normalize_cyl( )
  endif

  call sim%jay%update_boundary( sim%no_co, sim%bnd%send_vdf, sim%bnd%recv_vdf )
  call sim%jay%smooth()
  if ( associated(sim%jay1) ) then
    call sim%jay1%update_boundary( sim%no_co, sim%bnd%send_vdf, sim%bnd%recv_vdf )
    call sim%jay1%smooth()
  endif

  call sim % emf % advance( sim%jay, dt(sim%tstep), sim%bnd%send_vdf, sim%bnd%recv_vdf )

  if ( associated(sim%axf) ) then
    call sim%axf%advance_perturb_b( 0.5_p_double * dt(sim%tstep) )
    call sim%axf%apply_backreaction( sim%emf, 0.5_p_double * dt(sim%tstep), dt(sim%tstep) )
    call sim%axf%update_boundary( sim%no_co, sim%bnd%send_vdf, sim%bnd%recv_vdf )
    call sim%axf%report_energy( sim%no_co, sim%tstep, sim%emf, pump_omega0 )
  endif

  call end_event(loopev)

end subroutine iter_axion_hybrid


! Extract the strongest laser pulse parameters for diagnostic normalization.
subroutine strongest_zpulse_params( zpulse_list, omega0, a0 )
  type(t_zpulse_list), intent(in) :: zpulse_list
  real(p_double), intent(out) :: omega0, a0

  class(t_zpulse), pointer :: pulse
  real(p_double) :: score, best_score

  omega0 = 0.0_p_double
  a0 = 0.0_p_double
  best_score = -1.0_p_double

  pulse => zpulse_list%list
  do
    if ( .not. associated(pulse) ) exit
    score = abs( pulse%omega0 * pulse%a0 )
    if ( score > best_score ) then
      best_score = score
      omega0 = pulse%omega0
      a0 = pulse%a0
    endif
    pulse => pulse%next
  enddo
end subroutine strongest_zpulse_params

! Print the hybrid simulation algorithm summary, including the axion-field component.
subroutine list_algorithm_axion_hybrid( sim )
  class(t_simulation_axion_hybrid), intent(in) :: sim
  call sim%t_simulation_qed%list_algorithm()
  if ( associated(sim%axf) ) call sim%axf%list_algorithm()
end subroutine list_algorithm_axion_hybrid

end module m_simulation_axion_hybrid
