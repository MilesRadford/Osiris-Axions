
module m_walk_data_structures
	use m_system
	use m_parameters
	use m_time_step

	use m_grid_define
	use m_emf_define
	use m_particles_define
	use m_node_conf
	use m_grid
	use m_space
	use m_time
	use m_emf
	use m_particles
	use m_zpulse
	use m_antenna_array
	use m_current
	use m_vdf_define
	use m_simulation_cuda_define

	!use m_current_diag
	use iso_c_binding

	implicit none

	private

	interface setup_cuda
	  module procedure setup_cuda
	end interface

	interface

		!******************************************************************************************************************
		!**************** Interfaces to 'C' functions that setup Oriris's state in C accesiable structures ****************
		!-----------------------------------------------------------------------------------------
		subroutine cuda_basic_configuration(node_id, timestep_at_initalization, filename, filename_length) bind(C)
		  use iso_c_binding
		  integer(c_int), value,    intent(in)  :: node_id
		  integer(c_int), value,    intent(in)  :: timestep_at_initalization
		  character(kind=c_char),   intent(in)  :: filename(*)
		  integer(c_size_t), value, intent(in)  :: filename_length
		end subroutine

		subroutine init_cuda_nodeconf(no_num, nx, my_aid, my_ngp, neighbor, ifpr, comm) bind(C)
		  use iso_c_binding
		  integer(c_int), value,    intent(in)  :: no_num
		  type(c_ptr), value,    intent(in)  :: nx
		  integer(c_int), value,    intent(in)  :: my_aid
		  type(c_ptr), value,    intent(in)  :: my_ngp
		  type(c_ptr), value,    intent(in)  :: neighbor
		  type(c_ptr), value,    intent(in)  :: ifpr
		  integer(c_int), value,    intent(in)  :: comm
		end subroutine

		subroutine init_cuda_grid(nnodes, g_nx, min_nx, max_nx, my_nx, nx_node) bind(C)
		  use iso_c_binding
		  type(c_ptr), value,    intent(in)  :: nnodes
		  type(c_ptr), value,    intent(in)  :: g_nx
		  type(c_ptr), value,    intent(in)  :: min_nx
		  type(c_ptr), value,    intent(in)  :: max_nx
		  type(c_ptr), value,    intent(in)  :: my_nx
		  type(c_ptr), value,    intent(in)  :: nx_node
		end subroutine


		subroutine init_cuda_space(x_bnd, x_bnd_old,move_num, xmoved, total_i_moved, if_move, mov_vel, x_bnd_initial) bind(C)
		  use iso_c_binding
		  type(c_ptr), value,    intent(in)  :: x_bnd
		  type(c_ptr), value,    intent(in)  :: x_bnd_old
		  type(c_ptr), value,    intent(in)  :: move_num
		  type(c_ptr), value,    intent(in)  :: xmoved
		  type(c_ptr), value,    intent(in)  :: total_i_moved
		  type(c_ptr), value,    intent(in)  :: if_move
		  real(c_double), value, intent(in)  :: mov_vel
		  type(c_ptr), value,    intent(in)  :: x_bnd_initial
		end subroutine

		subroutine init_cuda_particles(num_species, num_cathode) bind(C)
		  use iso_c_binding
		  integer(c_int),     	intent(in)  :: num_species
		  integer(c_int), 			intent(in)  :: num_cathode
		end subroutine
!extern "C" void init_cuda_species(	int * species_index, int* particles_per_cell, double* g_box, double *dx, int *my_nx_p, int* coordinates,
!									int* num_par_max, int * num_par, int *num_created, float__fld *rqm, float__fld *rq_real, float__fld *q_real, float__fld *m_real,
!									int* interpolation, int* ix, float__fld *x, float__fld *p, float__fld *q, double* total_kinetic_4energy, double *total_number_particles_in_system_initially,
!									int* boundry_conditions__type, int* boundry_conditions__thermal_type,
!									float__part* boundry_conditions__pth_bnd, float__part* boundry_conditions__pfl_bnd, double* total_xmoved )

	subroutine init_cuda_species(species_index, particles_per_cell, g_box, dx, my_nx_p, coordinates, &
													num_par_max, num_par, num_created, rqm, rq_real, q_real, m_real, &
										interpolation, ix, x, p, q, total_kinetic_4energy, total_number_particles_in_system_initially, &
										boundry_conditions__type, boundry_conditions__thermal_type, &
										boundry_conditions__pth_bnd, boundry_conditions__pfl_bnd, total_xmoved) bind(C)
	  use iso_c_binding
	  use  m_parameters

	  integer(c_int),     intent(in)  :: species_index
	  type(c_ptr), value,    intent(in)  :: particles_per_cell
	  type(c_ptr), value,    intent(in)  :: g_box
	  type(c_ptr), value,    intent(in)  :: dx
	  type(c_ptr), value,    intent(in)  :: my_nx_p
	  integer(c_int),     intent(in)  :: coordinates
	  integer(c_int),    intent(in)  :: num_par_max
	  integer(c_int),     intent(in)  :: num_par
	  integer(c_int),     intent(in)  :: num_created
	  real(p_k_part),			intent(in) :: rqm
	  real(p_k_part),			intent(in) :: rq_real
	  real(p_k_part),			intent(in) :: q_real
	  real(p_k_part),			intent(in) :: m_real

	  integer(c_int),     intent(in)  :: interpolation
	  type(c_ptr), value,    intent(in)  :: ix
	  type(c_ptr), value,    intent(in)  :: x
	  type(c_ptr), value,    intent(in)  :: p
	  type(c_ptr), value,    intent(in)  :: q
	  type(c_ptr), value,    intent(in)  :: total_kinetic_4energy
	  real(c_double),			intent(in) :: total_number_particles_in_system_initially

	  type(c_ptr), value,    intent(in)  :: boundry_conditions__type
	  type(c_ptr), value,    intent(in)  :: boundry_conditions__thermal_type
	  type(c_ptr), value,    intent(in)  :: boundry_conditions__pth_bnd
	  type(c_ptr), value,    intent(in)  :: boundry_conditions__pfl_bnd
	  type(c_ptr), value,    intent(in)  :: total_xmoved

	end subroutine

	subroutine init_species_diagnostics(ndump_fac_ene, ndump_fac_temp, &
                                ndump_fac_raw, ndump_fac_tracks, &
                                phase_space__ndump_fac, phase_space__ndump_fac_tavg, &
                                phase_space__ndump_n_tavg, species_index)  bind(C)
		use iso_c_binding
		integer(c_int),     intent(in)  :: ndump_fac_ene
		integer(c_int),     intent(in)  :: ndump_fac_temp
		integer(c_int),     intent(in)  :: ndump_fac_raw
		integer(c_int),     intent(in)  :: ndump_fac_tracks
		integer(c_int),     intent(in)  :: phase_space__ndump_fac
		integer(c_int),     intent(in)  :: phase_space__ndump_fac_tavg
		integer(c_int),     intent(in)  :: phase_space__ndump_n_tavg
		integer(c_int),     intent(in)  :: species_index
	end subroutine

	subroutine init_cuda_current(coordinates, gix_pos, bc_type, interpolation, &
				type, order, max_order,scoef, digital_A,digital_fc ) bind(c)
		use iso_c_binding
		integer(c_int),     intent(in)  :: coordinates
		type(c_ptr), value,    intent(in)  :: gix_pos
		type(c_ptr), value,    intent(in)  :: bc_type
		integer(c_int),     intent(in)  :: interpolation

		type(c_ptr), value,    intent(in)  :: type
		type(c_ptr), value,    intent(in)  :: order
		integer(c_int),     intent(in)  :: max_order
		type(c_ptr), value,    intent(in)  :: scoef
		type(c_ptr), value,    intent(in)  :: digital_A
		type(c_ptr), value,    intent(in)  :: digital_fc
	end subroutine

	subroutine init_cuda_current_vdf(	jay_xdim, jay_fdim, jay_f1,  jay_f2,  jay_f3, jay_nx, jay_gc_num, jay_dx) bind(C)
																!ifsmooth, smooth_level, nsm, omega0, dir, swfj, scoef) bind(c)
		use iso_c_binding
		integer(c_int),     intent(in)  :: jay_xdim
		integer(c_int),     intent(in)  :: jay_fdim
		type(c_ptr), value,    intent(in)  :: jay_f1
		type(c_ptr), value,    intent(in)  :: jay_f2
		type(c_ptr), value,    intent(in)  :: jay_f3
		type(c_ptr), value,    intent(in)  :: jay_nx
		type(c_ptr), value,    intent(in)  :: jay_gc_num
		type(c_ptr), value,    intent(in)  :: jay_dx

		!type(c_ptr), value,     intent(in)  :: ifsmooth
		!type(c_ptr), value,     intent(in)  :: smooth_level
		!type(c_ptr), value,     intent(in)  :: nsm
		!real(c_double),     intent(in)  :: omega0
		!real(c_int),     intent(in)  :: dir
		!type(c_ptr), value,     intent(in)  :: swfj
		!type(c_ptr), value,     intent(in)  :: scoef
	end subroutine

	subroutine init_cuda_emf(	e_xdim, e_fdim, e_f1,  e_f2,  e_f3, e_nx, e_gc_num, e_dx, &
										b_xdim, b_fdim,  b_f1,  b_f2,  b_f3, b_nx, b_gc_num, b_dx, &
										epart_xdim, epart_fdim,  epart_f1,  epart_f2,  epart_f3, epart_nx, epart_gc_num, epart_dx, &
										bpart_xdim, bpart_fdim,  bpart_f1,  bpart_f2,  bpart_f3, bpart_nx, bpart_gc_num, bpart_dx, &
										exte_xdim, exte_fdim,  exte_f1,  exte_f2,  exte_f3, exte_nx, exte_gc_num, exte_dx, &
										extb_xdim, extb_fdim,  extb_f1,  extb_f2,  extb_f3, extb_nx, extb_gc_num, extb_dx, &
										solver, coordinates, & !e_scratch_f1, e_scratch_f2, e_scratch_f3, b_scratch_f1,b_scratch_f2,  b_scratch_f3,
										ifextfld, emf_bound_types, &
										smooth_type, smooth_niter, smooth_nmax, &
										type, order, max_order,scoef, digital_A,digital_fc) bind(C)
										!ifsmooth, smooth_level, nsm, &
										!omega0, dir, swfj,&
										!scoef) bind(C)
		use iso_c_binding
		integer(c_int),     intent(in)  :: e_xdim
		integer(c_int),     intent(in)  :: e_fdim
		type(c_ptr), value,    intent(in)  :: e_f1
		type(c_ptr), value,    intent(in)  :: e_f2
		type(c_ptr), value,    intent(in)  :: e_f3
		type(c_ptr), value,    intent(in)  :: e_nx
		type(c_ptr), value,    intent(in)  :: e_gc_num
		type(c_ptr), value,    intent(in)  :: e_dx

		integer(c_int),     intent(in)  :: b_xdim
		integer(c_int),     intent(in)  :: b_fdim
		type(c_ptr), value,    intent(in)  :: b_f1
		type(c_ptr), value,    intent(in)  :: b_f2
		type(c_ptr), value,    intent(in)  :: b_f3
		type(c_ptr), value,    intent(in)  :: b_nx
		type(c_ptr), value,    intent(in)  :: b_gc_num
		type(c_ptr), value,    intent(in)  :: b_dx

		integer(c_int),     intent(in)  :: epart_xdim
		integer(c_int),     intent(in)  :: epart_fdim
		type(c_ptr), value,    intent(in)  :: epart_f1
		type(c_ptr), value,    intent(in)  :: epart_f2
		type(c_ptr), value,    intent(in)  :: epart_f3
		type(c_ptr), value,    intent(in)  :: epart_nx
		type(c_ptr), value,    intent(in)  :: epart_gc_num
		type(c_ptr), value,    intent(in)  :: epart_dx

		integer(c_int),     intent(in)  :: bpart_xdim
		integer(c_int),     intent(in)  :: bpart_fdim
		type(c_ptr), value,    intent(in)  :: bpart_f1
		type(c_ptr), value,    intent(in)  :: bpart_f2
		type(c_ptr), value,    intent(in)  :: bpart_f3
		type(c_ptr), value,    intent(in)  :: bpart_nx
		type(c_ptr), value,    intent(in)  :: bpart_gc_num
		type(c_ptr), value,    intent(in)  :: bpart_dx

		integer(c_int),     intent(in)  :: exte_xdim
		integer(c_int),     intent(in)  :: exte_fdim
		type(c_ptr), value,    intent(in)  :: exte_f1
		type(c_ptr), value,    intent(in)  :: exte_f2
		type(c_ptr), value,    intent(in)  :: exte_f3
		type(c_ptr), value,    intent(in)  :: exte_nx
		type(c_ptr), value,    intent(in)  :: exte_gc_num
		type(c_ptr), value,    intent(in)  :: exte_dx

		integer(c_int),     intent(in)  :: extb_xdim
		integer(c_int),     intent(in)  :: extb_fdim
		type(c_ptr), value,    intent(in)  :: extb_f1
		type(c_ptr), value,    intent(in)  :: extb_f2
		type(c_ptr), value,    intent(in)  :: extb_f3
		type(c_ptr), value,    intent(in)  :: extb_nx
		type(c_ptr), value,    intent(in)  :: extb_gc_num
		type(c_ptr), value,    intent(in)  :: extb_dx

		integer(c_int),     intent(in)  :: solver
		integer(c_int),     intent(in)  :: coordinates
		integer(c_int),     intent(in)  :: ifextfld
		type(c_ptr), value,    intent(in)  :: emf_bound_types
		integer(c_int),     intent(in)  :: smooth_type
		integer(c_int),     intent(in)  :: smooth_niter
		integer(c_int),     intent(in)  :: smooth_nmax

		type(c_ptr), value,    intent(in)  :: type
		type(c_ptr), value,    intent(in)  :: order
		integer(c_int),     intent(in)  :: max_order
		type(c_ptr), value,    intent(in)  :: scoef
		type(c_ptr), value,    intent(in)  :: digital_A
		type(c_ptr), value,    intent(in)  :: digital_fc

		!type(c_ptr), value,     intent(in)  :: ifsmooth
		!type(c_ptr), value,     intent(in)  :: smooth_level
		!type(c_ptr), value,     intent(in)  :: nsm
		!real(c_double),     intent(in)  :: omega0
		!real(c_int),     intent(in)  :: dir
		!type(c_ptr), value,     intent(in)  :: swfj
		!type(c_ptr), value,,     intent(in)  :: scoef
		end subroutine

		subroutine check_c_state() bind(C)
		end subroutine

		subroutine osiris_cuda__quickening() bind(C)
			use iso_c_binding
		end subroutine

	end interface
	public :: setup_cuda

	contains

	subroutine setup_cuda( sim  )

		implicit none
		class(t_simulation_cuda) :: sim

		integer :: initial_time_step, i, j
		real(c_double) 								::  total_par
		type(c_funptr) :: cproc
		integer, dimension( 2, p_max_dim ) :: neighbor_array
		type(t_vdf) :: temp_vdf
		!real(p_double), dimension(2,p_max_dim) :: x_bnd_old__temp
		!real(p_double), dimension(p_max_dim) :: x_bnd_old__temp

		initial_time_step = n( sim%tstep )
		call cuda_basic_configuration(	sim%no_co%my_aid(), &
																		n(sim%tstep), &
																		trim( sim%options%input_file ), &
																		len(trim( sim%options%input_file ), &
																		kind=c_size_t) )

		do i=1,p_max_dim
		  neighbor_array( 1, i) = neighbor( sim%no_co, 1, i)
		  neighbor_array( 2, i) = neighbor( sim%no_co, 2, i)
		end do

		call init_cuda_nodeconf (	no_num(sim%no_co), &
															c_loc( nx(sim%no_co) ), &
															sim%no_co%my_aid(), &
															c_loc( my_ngp(sim%no_co) ), &
		                          c_loc(neighbor_array), &
		                          c_loc(periodic( sim%no_co )), &
		                          comm( sim%no_co) )

		call init_cuda_grid(			c_loc( sim%grid%nnodes), &
															c_loc(sim%grid%g_nx), &
															c_loc(sim%grid%min_nx), &
		                          c_loc(sim%grid%max_nx), &
		                          c_loc(sim%grid%my_nx), &
		                          c_loc( sim%grid%nx_node) )



		call init_cuda_space(			c_loc(sim%g_space ), &
															c_loc(sim%g_space%x_bnd_old), &
															c_loc(sim%g_space%move_num), &
															c_loc(sim%g_space%xmoved), &
															c_loc(sim%g_space%total_i_moved), &
															c_loc(sim%g_space%if_move), &
															sim%g_space%move_vel, &
															c_loc(sim%g_space%x_bnd_initial) )

		call init_cuda_particles( sim%part%num_species, sim%part%num_cathode )

		do i=1, sim%part%num_species
			total_par = real( sim%part%species(i)%num_par )
			!call sum_up( sim%no_co, total_par )

	    call init_cuda_species(sim%part%species(i)%sp_id, c_loc(sim%part%species(i)%num_par_x), c_loc(sim%part%species(i)%g_box), c_loc(sim%part%species(i)%dx), &
								    c_loc(sim%part%species(i)%my_nx_p), sim%part%species(i)%coordinates, &
								    sim%part%species(i)%num_par_max,sim%part%species(i)%num_par,sim%part%species(i)%num_created, &
                                    sim%part%species(i)%rqm,sim%part%species(i)%rq_real,sim%part%species(i)%q_real,sim%part%species(i)%m_real, &
								    sim%part%species(i)%interpolation,c_loc(sim%part%species(i)%ix),&
									c_loc(sim%part%species(i)%x),c_loc(sim%part%species(i)%p),c_loc(sim%part%species(i)%q), &
									c_loc(sim%part%species(i)%energy), total_par, &
									c_loc(sim%part%species(i)%bnd_con%type), c_loc(sim%part%species(i)%bnd_con%thermal_type), &
									c_loc(sim%part%species(i)%bnd_con%uth_bnd), c_loc(sim%part%species(i)%bnd_con%ufl_bnd), &
									c_loc(sim%part%species(i)%total_xmoved ))

	    call init_species_diagnostics(  sim%part%species(i)%diag%ndump_fac_ene, &
	                                    sim%part%species(i)%diag%ndump_fac_temp, &
	                                    sim%part%species(i)%diag%ndump_fac_raw, &
	                                    sim%part%species(i)%diag%ndump_fac_tracks, &
	                                    sim%part%species(i)%diag%phasespaces%ndump_fac, &
	                                    sim%part%species(i)%diag%phasespaces%ndump_fac_tavg, &
	                                    sim%part%species(i)%diag%phasespaces%n_tavg, &
	                                    i)

		end do


	   call init_cuda_emf(  &
	   	sim%emf%e%x_dim, sim%emf%e%f_dim, c_loc(sim%emf%e%f1), c_loc(sim%emf%e%f2), c_loc(sim%emf%e%f3), c_loc(sim%emf%e%nx), c_loc(sim%emf%e%gc_num), c_loc(sim%emf%e%dx), &
	    sim%emf%b%x_dim, sim%emf%b%f_dim, c_loc(sim%emf%b%f1), c_loc(sim%emf%b%f2), c_loc(sim%emf%b%f3), c_loc(sim%emf%b%nx), c_loc(sim%emf%b%gc_num), c_loc(sim%emf%b%dx), &
			sim%emf%e_part%x_dim, sim%emf%e_part%f_dim, c_loc(sim%emf%e_part%f1), c_loc(sim%emf%e_part%f2), c_loc(sim%emf%e_part%f3), c_loc(sim%emf%e_part%nx), c_loc(sim%emf%e_part%gc_num), c_loc(sim%emf%e_part%dx), &
			sim%emf%b_part%x_dim, sim%emf%b_part%f_dim, c_loc(sim%emf%b_part%f1), c_loc(sim%emf%b_part%f2), c_loc(sim%emf%b_part%f3), c_loc(sim%emf%b_part%nx), c_loc(sim%emf%b_part%gc_num), c_loc(sim%emf%b_part%dx), &
			sim%emf%ext_e%x_dim, sim%emf%ext_e%f_dim, c_loc(sim%emf%ext_e%f1), c_loc(sim%emf%ext_e%f2), c_loc(sim%emf%ext_e%f3), c_loc(sim%emf%ext_e%nx), c_loc(sim%emf%ext_e%gc_num), c_loc(sim%emf%ext_e%dx), &
			sim%emf%ext_b%x_dim, sim%emf%ext_b%f_dim, c_loc(sim%emf%ext_b%f1), c_loc(sim%emf%ext_b%f2), c_loc(sim%emf%ext_b%f3), c_loc(sim%emf%ext_b%nx), c_loc(sim%emf%ext_b%gc_num), c_loc(sim%emf%ext_b%dx), &
			sim%emf%solver, sim%emf%coordinates, & !sim%emf%e_scratch%f1, sim%emf%e_scratch%f2, sim%emf%e_scratch%f3 ,sim%emf%b_scratch%f1, sim%emf%b_scratch%f2, sim%emf%b_scratch%f3,&
			sim%emf%ext_fld, c_loc(sim%emf%bnd_con%type), sim%emf%smooth_type, sim%emf%smooth_niter, sim%emf%smooth_nmax, &
			c_loc(sim%emf%smooth%type),c_loc(sim%emf%smooth%order),sim%emf%smooth%max_order,c_loc(sim%emf%smooth%scoef),c_loc(sim%emf%smooth%digital_A),c_loc(sim%emf%smooth%digital_fc) )

			call init_cuda_current(sim%jay%coordinates, c_loc(sim%jay%gix_pos), c_loc(sim%jay%bc_type), sim%jay%interpolation, &
						c_loc(sim%jay%smooth%type),c_loc(sim%jay%smooth%order),sim%jay%smooth%max_order,c_loc(sim%jay%smooth%scoef),c_loc(sim%jay%smooth%digital_A),c_loc(sim%jay%smooth%digital_fc) )
		  if( ASSOCIATED(sim%jay%pf ) ) then
			  do i=1,size(sim%jay%pf)
			  	temp_vdf = sim%jay%pf(i)
			  	call init_cuda_current_vdf(	temp_vdf%x_dim, temp_vdf%f_dim, c_loc(temp_vdf%f1), c_loc(temp_vdf%f2), c_loc(temp_vdf%f3), c_loc(temp_vdf%nx), c_loc(temp_vdf%gc_num), c_loc(temp_vdf%dx))
				end do
			end if

		call check_c_state()

		! Now do the main GPU code setup/inialization
		!		And tell the CUDA code the inital timestep... this is to identify if we have
		!		restarted a simulation.
		call osiris_cuda__quickening( )

	end subroutine
end module m_walk_data_structures
