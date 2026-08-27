
#include "os-cuda-state.h"

// the container that holds all the CUDA C state.
t_simulation_cuda CUDA_STATE( createInputDeckFactory() );
my_gpu_environment gpu_environment;

t_simulation_cuda &getCurrentState() {
	return CUDA_STATE;
}

// Makes assumption that node_id = my_aid - 1 (my_aid is and Osiris Variable)
//
extern "C" void cuda_basic_configuration( int node_id, int timestep_at_initalization, char* input_deck_filename, size_t length) {
	int node_index = node_id;
	
	CUDA_STATE.CUDA_CONFIG->timestep_at_initalization = timestep_at_initalization;

	std::string nose;
	std::string filename = std::string( input_deck_filename, length);
	
	CUDA_STATE.CUDA_CONFIG->set_node_id( node_index - 1 );
	
	if( CUDA_STATE.CUDA_CONFIG->is_root() ) {
		std::cout << "   Obtained the input-deck filename '" << filename << "' from Osiris." << std::endl;
	}
	CUDA_STATE.CUDA_CONFIG->read_configuration_from_file(  filename );
	
	if( CUDA_STATE.CUDA_CONFIG->is_root() ) {
		std::cout << "   Cuda Configured to use a GPU 'block size' of " << CUDA_STATE.CUDA_CONFIG->default_cuda_block_size << " threads." << std::endl;
		std::cout << "   Cuda Configured to use a tile size of (" << CUDA_STATE.CUDA_CONFIG->tile_size__x << " cells," << CUDA_STATE.CUDA_CONFIG->tile_size__y << " cells) " << std::endl;
	}
}


void init_cuda_nodeconf_cpp(int no_num, int*nx,  int my_aid, int*my_ngp, int*neighbor, int*ifpr, int mpi_comm,  t_node_conf& v) {
	int i;

	v.x_dim = p_x_dim;
	v.mpi_comm = mpi_comm;																	//cout << "nodeconf.no_num (total number nodes): " << v.no_num << std::std::endl;
	v.no_num = no_num;																		//cout << "nodeconf.nx (nodes along each sptial direction):" << std::std::endl;

	for(i=0;i<p_max_dim;i++) {
		v.nx[i] = nx[i];																	//cout << "\t\t" << nx[i] << std::std::endl;
	}
																																							
	v.my_aid = my_aid;																		//cout << "v.my_aid (bacially Fortran index version of mpi_rank (so rank+1): " << v.my_aid << std::std::endl;
	v.my_MPI_rank = v.my_aid - 1;															//cout << "My MPI rank (TODO: get this excplitly not guess): " << v.my_MPI_rank << std::std::endl;

	for(i=0;i<p_max_dim;i++) {
		v.my_ngp[i] = my_ngp[i];															//cout << "my_ngp ( N_ode G_rid P_osition (within the grid of nodes given by [nx(1),nx(2),nx(3)] )):" << std::std::endl;
																							//cout << "\t\t" << my_ngp[i] << std::std::endl;
	}
	
	for(i=0; i < p_max_dim; i++) {															//cout << "neighbor: provides list of neighboring cells for each dim neighbor[p_max_dim[0=low,1=upper]" << std::std::endl;
																							//cout << "MPI " << v.my_MPI_rank << ".. for node " << v.my_aid << ":" << "neighbors x: " << v.neighbors[0][0] << " and " << v.neighbors[0][1] << ". Along Y: " << v.neighbors[1][0] << " and " << v.neighbors[1][1] << std::std::endl;
		v.neighbors[i][0] = neighbor[i*2 + 0];
		v.neighbors[i][1] = neighbor[i*2 + 1];
	}
																																	
	for(i=0;i<p_max_dim;i++) {																//cout << "ifpr (flags for periodic node grid (periodic boundary conditions)):" << std::std::endl;
		v.ifpr[i] = ifpr[i];																//cout << "\t\t"  << ifpr[i] << std::std::endl;
	}
	
	/*v.ngp = new int[p_x_dim*v.no_num];														//cout << "ngp (ngp of all processes):" << std::std::endl;
	if( v.no_num > 1) {																		
		for(i=0;i<(p_x_dim*v.no_num);i++) {
			v.ngp[i] = ngp[i];																//cout << "\t\t"  << ngp[i] << std::std::endl;
		}
	} */
}

extern "C" void init_cuda_nodeconf(int no_num, int*nx,  int my_aid, int*my_ngp, int*neighbor, int*ifpr, int mpi_comm) {
		init_cuda_nodeconf_cpp( no_num, nx, my_aid, my_ngp, neighbor, ifpr, mpi_comm, CUDA_STATE.no_co);
}


void init_cuda_space_cpp(	double *x_bnd, double *x_bnd_old, int  *move_num,
								double *xmoved, int * total_i_moved, int * if_move,
								double *mov_vel, double * x_bnd_initial, t_space& v) {
	
	//int i,j;
	
	#if 0
	v.x_dim = p_x_dim;
	
	for(i=0;i<p_max_dim;i++) {
		
		v.x_bnd[i][0] = x_bnd[i*2 + 0];
		v.x_bnd[i][1] = x_bnd[i*2 + 1];
		
		v.x_bnd_old[i][0] = x_bnd_old[i*2 + 0];
		v.x_bnd_old[i][1] = x_bnd_old[i*2 + 1];
		
		v.move_num[i][0] = move_num[i*2 + 0];
		v.move_num[i][1] = move_num[i*2 + 1];
		
		v.xmoved[i][0] = xmoved[i*2 + 0];
		v.xmoved[i][1] = xmoved[i*2 + 1];
		
		v.x_bnd_initial[i][0] = x_bnd_initial[i*2 + 0];
		v.x_bnd_initial[i][1] = x_bnd_initial[i*2 + 1];
	
		v.total_i_moved[i] = total_i_moved[i];
		v.if_move[i] = if_move[i];
		
		// LOG << "x_bnd[" << i << "]:\t\t\t" <<  v.x_bnd[i][0] << "," << v.x_bnd[i][1] << std::std::endl;
		// LOG << "x_bnd_old[" << i <<"]:\t\t\t" << v.x_bnd_old[i][0] << "," << v.x_bnd_old[i][1] << std::std::endl;
		// LOG << "x_x_bnd_initial[" << i <<"]:\t\t\t" << v.x_bnd_initial[i][0] << "," << v.x_bnd_initial[i][1] << std::std::endl;
		
		/*
		cout << "x_bnd[" << i <<"]: " << v.x_bnd[i][0] << "," << v.x_bnd[i][1] << std::std::endl;
		cout << "x_bnd_old[" << i <<"]: " << v.x_bnd_old[i][0] << "," << v.x_bnd_old[i][1] << std::std::endl;
		cout << "x_x_bnd_initial[" << i <<"]: " << v.x_bnd_initial[i][0] << "," << v.x_bnd_initial[i][1] << std::std::endl;
		*/
	}
	
	#endif
	
	//v.x_bnd = *reinterpret_cast<double (*)[p_max_dim][2]>(x_bnd); // ( double(*)[p_max_dim][2] ) x_bnd;
	v._x_bnd = ( double(*)[p_max_dim][2] ) x_bnd;
	v._x_bnd_old = ( double(*)[p_max_dim][2] ) x_bnd_old;
	v._move_num = ( int (*)[p_max_dim][2] ) move_num;
	v._xmoved = ( double(*)[p_max_dim][2] ) xmoved;
	v._total_i_moved = ( int(*)[p_max_dim][2] ) total_i_moved;
	v._if_move = ( int (*)[p_max_dim] ) if_move;
	v._mov_vel = ( double(*)[p_max_dim] ) mov_vel;
	v._x_bnd_initial = ( double(*)[p_max_dim][2] ) x_bnd_initial;

	/*
	// assigning points to c++ arrays...
	double fred = 1.0;
	
	
	double (&shot) [p_max_dim][2] = *reinterpret_cast<double (*)[p_max_dim][2]>(x_bnd);
	fred = shot[1][1];
	
	double (*shot1) [p_max_dim][2] = reinterpret_cast<double (*)[p_max_dim][2]>(x_bnd);
	fred = shot1[0][1][1];
	
	*/
	
	
	
	/*
	
		s->my_nx_p[i][0] = my_nx_p[3*i + 0];
		s->my_nx_p[i][1] = my_nx_p[3*i + 1];
		s->my_nx_p[i][2] = my_nx_p[3*i + 2];
	}
	
	for(i=0;i<p_x_dim;i++) {
		s->g_box[i][0] = g_box[i*2 + 0];
		s->g_box[i][1] = g_box[i*2 + 1];
	
	*/
}

extern "C" void init_cuda_space(	double *x_bnd, double *x_bnd_old, int  *move_num,
								double *xmoved, int * total_i_moved, int * if_move,
								double *mov_vel, double * x_bnd_initial) {

	t_space &g_space = getCurrentState().g_space;
	init_cuda_space_cpp(x_bnd, x_bnd_old,move_num,
											xmoved, total_i_moved, if_move,
											mov_vel, x_bnd_initial, g_space);
}

void init_cuda_grid_cpp(int * nnodes, int * g_nx, int * min_nx, int * max_nx, int* my_nx, int* nx_node, t_grid& v) {
	int i;

	v.x_dim = p_x_dim;
	v.lb_type = 0;
	v.dynamic = 0;
	v.max_unbalance = 0.0;
	v.cell_weight = 0.0;
	// was p_max_dim which is wrong
																							// LOG << "nnodes:" << std::std::endl;
	for(i=0;i<p_x_dim;i++) {
		v.nnodes[i] = nnodes[i];			// total number of 'nodes' in each direction
																							// LOG << "\t\t"  << nnodes[i] << "\t";
	}
			 																				// LOG << std::std::endl;
																							// LOG << "g_nx:" << std::std::endl;
	for(i=0;i<p_x_dim;i++) {
		v.g_nx[i] = g_nx[i];				// global number of cells in each direction.
																							// LOG << "\t\t"  << g_nx[i] <<  "\t";
	}
																							// LOG << std::std::endl;
																							// LOG << "min_nx:" << std::std::endl;
	for(i=0;i<p_x_dim;i++) {
		v.min_nx[i] = min_nx[i];			//  minimum number of cells in each partition (not used.. take out)
																							// LOG << "\t\t"  << max_nx[i] << "\t";
	}
																							// LOG << std::std::endl;
																							// LOG << "min_nx:" << std::std::endl;
	for(i=0;i<p_x_dim;i++) {
		v.max_nx[i] = max_nx[i];			//  maximum number of cells in each partition (not used.. take out)
																							// LOG << "\t\t"  << min_nx[i] << "\t";
	}
			
																							// LOG << std::std::endl;
	/*
	if(v.x_dim >= 1) {
		for(int i=0;i< v.nnodes[0]; i++) {
			vector<int> new_row;
			new_row.push_back(x1_nx[3*i + 0]);
			new_row.push_back(x1_nx[3*i + 1]);
			new_row.push_back(x1_nx[3*i + 2]);
			//cout << " X for node " << i << ": " << x1_nx[3*i + 0] << "," << x1_nx[3*i + 1] << "," << x1_nx[3*i + 2] << std::std::endl;
			v.x1_nx.push_back(new_row);
		}
	}
	
	if(v.x_dim >= 2) {
		for(int i=0;i< v.nnodes[1]; i++) {
			vector<int> new_row;
			new_row.push_back(x2_nx[3*i + 0]);
			new_row.push_back(x2_nx[3*i + 1]);
			new_row.push_back(x2_nx[3*i + 2]);
			//cout << " Y for node " << i << ": " << x2_nx[3*i + 0] << "," << x2_nx[3*i + 1] << "," << x2_nx[3*i + 2] << std::std::endl;
			v.x2_nx.push_back(new_row);
		}
	}	
	
	if(v.x_dim == 3) {
		for(int i=0;i< v.nnodes[2]; i++) {
			vector<int> new_row;
			new_row.push_back(x3_nx[3*i + 0]);
			new_row.push_back(x3_nx[3*i + 1]);
			new_row.push_back(x3_nx[3*i + 2]);
			v.x3_nx.push_back(new_row);
		}
	}*/
}

extern "C" void init_cuda_grid(int * nnodes, int * g_nx, int * min_nx, int * max_nx, int* my_nx, int * nx_node) {
	init_cuda_grid_cpp( nnodes, g_nx, min_nx, max_nx, my_nx, nx_node, CUDA_STATE.grid);
}

void init_cuda_particles_cpp(int * number_of_species, int * number_of_hybrid, t_simulation_cuda STATE) {
	STATE.parts.num_species = *number_of_species;
	STATE.parts.num_hybrid = *number_of_hybrid;
}

extern "C" void  init_cuda_particles(int * number_of_species, int * number_of_hybrid) {
	init_cuda_particles_cpp(number_of_species, number_of_hybrid, getCurrentState() );
}

void init_cuda_species_cpp(	int * species_index, int* particles_per_cell, double* g_box, double *dx, int *my_nx_p, int* coordinates,
									int* num_par_max, int * num_par, int *num_created, float__part *rqm, float__part *rq_real, float__part *q_real, float__part *m_real,
									int* interpolation, int* ix, float__part *x, float__part *p, float__part *q, double* total_kinetic_4energy, double *total_number_particles_in_system_initially,
									int* boundry_conditions__type, int* boundry_conditions__thermal_type, 
									float__part* boundry_conditions__pth_bnd, float__part* boundry_conditions__pfl_bnd, double* total_xmoved, t_simulation_cuda& STATE )
{
	
	int s_index;
	int i;

	t_species s_obj;
	
	// take care of index mapping from Fortran
	s_index = (*species_index) - 1;
	STATE.parts.species.push_back(s_obj);
	
	t_species& s = STATE.parts.species[s_index];
	

	// make a species object on the stack to fill in
	//--t_species s;

																								// LOG << "species " << s_index << ":" << std::std::endl;
	s.sp_id = s_index;
																								// LOG << "\tdx: (cell size)"<< std::std::endl;
	for(i=0;i<3;i++) {
		s.dx[i] = dx[i];
		//out << "<<<<<<<<<<<<< dx["<<i<<"]=" << dx[i] << std::std::endl;
																								// LOG << "\t\t"<< dx[i] << std::std::endl;
	}
	
																								// LOG << "\tnum_par_x (ppc):"<< std::std::endl;
	for(i=0;i<p_x_dim;i++) {
		s.num_par_x[i] = particles_per_cell[i];
																								// LOG << "\t\t"<< s.num_par_x[i] << std::std::endl;																	
	}

																								// LOG << "\tmy_nx_p: (local grid position in global grid)"<< std::std::endl;
	for(i=0;i<3;i++) {
		s.my_nx_p[i][0] = my_nx_p[3*i + 0];
		s.my_nx_p[i][1] = my_nx_p[3*i + 1];
		s.my_nx_p[i][2] = my_nx_p[3*i + 2];
																								// LOG << "\t\tfortan indices: "<< (i+1) << "space direction): << lower index:" <<  my_nx_p[3*i + 0] << ",upper index:" << my_nx_p[3*i + 1] << ",length:" << my_nx_p[3*i + 2] << std::std::endl;			
																								// LOG << "\t\t     c indices: "<< (i) << "space direction): << lower index:" <<  (my_nx_p[3*i + 0]-1) << ",upper index:" << (my_nx_p[3*i + 1]-1) << ",length:" << my_nx_p[3*i + 2] << std::std::endl;			
	}

																								// LOG << "\tg_box: (simulation box dimensions)"<< std::std::endl;
	//for(i=0;i<p_x_dim;i++) {
	//	s.g_box[i][0] = g_box[i*2 + 0];
	//	s.g_box[i][1] = g_box[i*2 + 1];
	//																							// LOG << "\t\t"<< g_box[2*i + 0] << "," << g_box[2*i + 1] << std::std::endl;					
	//}
	
	s._g_box = ( double(*)[p_max_dim][2] ) g_box;
	
	//s._total_moved =  ( double(*)[p_max_dim][2] ) total_moved;
	s._total_xmoved =  ( double(*)[p_x_dim] )  total_xmoved;
	
	s.coordinates = *coordinates; 
	
	s.num_par_max = *num_par_max;
	s.num_par_max__ptr = num_par_max;
	
	s.number_of_particles = *num_par;
																								// LOG << "\tnumber_of_particles:" << s.number_of_particles << std::std::endl;
	s.total_number_particles_in_system_initially = *total_number_particles_in_system_initially;
																								// LOG << "\tnumber_of_particles in sysmtes (sum of 'number_of_particles' across all nodes):" << s.total_number_particles_in_system_initially << std::std::endl;
																								
	s.num_created = *num_created;
																								// LOG << "s.num_created :" <<  s.num_created <<std::std::endl;
	s.rqm = *rqm;
																								// LOG << "s.rqm :"<< s.rqm <<std::std::endl;
	s.rq_real = *rq_real;
																								// LOG << "s.rq_real: " << s.rq_real <<std::endl;
	s.q_real = *q_real;
																								// LOG << "s.q_real: " << s.q_real <<std::endl;
	s.m_real = *m_real;
																								// LOG << "s.m_real: " << s.m_real  <<std::endl;
	s.interpolation = *interpolation;
																								// LOG << "\tinterpolation:" << s.interpolation << std::endl;
	s.ix = ix;
	s.x = x;
	s.p = p;
	s.q = q;
    s.total_kinetic_4energy = total_kinetic_4energy;
	s.num_par = num_par;
	
	
  
	
	//cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
	//cout << s.total_kinetic_4energy[0] << std::endl;
	//cout << s.total_kinetic_4energy[1] << std::endl;
	//cout << s.total_kinetic_4energy[2] << std::endl;
	//cout << s.total_kinetic_4energy[3] << std::endl;
	//cout << s.total_kinetic_4energy[3] << std::endl;
	//cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
	
	//s.is_energy = false;								// TODO: WTF is this?
	//s.pos_type = species_pos_type___p_low;				// TODO: WTF ids this?
	
	// create buffer for the transpose objects...
	//s.p_transpose = (double*) malloc(s.num_par_max *  sizeof( double ) * 3);
	//s.p_transpose = new particle_float_type[s.num_par_max *  sizeof( particle_float_type ) * 3];


	// now add this to the 'particles' object. (This is a global object)
	//-- parts.species.push_back(s);
	
	
	// now read in the particle's boundry conditions...
	//cout << "slupring in species boundires!!!!!! " << p_max_dim << std::endl << std::endl;
	for(int i=0;i< p_max_dim;i++) {
		STATE.parts.species.back().bnd_con.type[i][0] = boundry_conditions__type[2*i + 0]; 
		STATE.parts.species.back().bnd_con.type[i][1] = boundry_conditions__type[2*i + 1]; 
		
		//cout << " for component " << i << " we gots " << s.bnd_con.type[i][0] << "," << s.bnd_con.type[i][1] << std::endl;
		STATE.parts.species.back().bnd_con.thermal_type[i][0] = boundry_conditions__thermal_type[2*i + 0]; 
		STATE.parts.species.back().bnd_con.thermal_type[i][1] = boundry_conditions__thermal_type[2*i + 1]; 
		
		for(int j=0;j<p_p_dim;j++) {
			STATE.parts.species.back().bnd_con.pth_bnd[i][0][j] = boundry_conditions__pth_bnd[2*p_p_dim*i+0*p_p_dim + j];
			STATE.parts.species.back().bnd_con.pth_bnd[i][1][j] = boundry_conditions__pth_bnd[2*p_p_dim*i+1*p_p_dim + j];
		}
	}
	
	s.is_energy = false;								// TODO: WTF is this?
	s.pos_type = 1; //p_low				// TODO: WTF ids this?

	/*
	// print out stuff to test...
	for(i=0;i<p_x_dim;i++) {
		printf("g_box[dim1=%d][0] = %e", i, s.g_box[i][0]);
		printf(" g_box[dim1=%d][1] = %e\n", i, s.g_box[i][1]);
	}
	
	for(i=0;i<3;i++) {
		printf("my_nx_p[dim1=%d][0] = %d\n", i, s.my_nx_p[i][0]);
		printf("my_nx_p[dim1=%d][1] = %d\n", i, s.my_nx_p[i][1]);
		printf("my_nx_p[dim1=%d][2] = %d\n", i, s.my_nx_p[i][2]);
	} 
	*/
}

extern "C" void init_cuda_species(	int * species_index, int* particles_per_cell, double* g_box, double *dx, int *my_nx_p, int* coordinates,
									int* num_par_max, int * num_par, int *num_created, float__fld *rqm, float__fld *rq_real, float__fld *q_real, float__fld *m_real,
									int* interpolation, int* ix, float__fld *x, float__fld *p, float__fld *q, double* total_kinetic_4energy, double *total_number_particles_in_system_initially,
									int* boundry_conditions__type, int* boundry_conditions__thermal_type, 
									float__part* boundry_conditions__pth_bnd, float__part* boundry_conditions__pfl_bnd, double* total_xmoved ) {
	init_cuda_species_cpp(	species_index, particles_per_cell, g_box, dx, my_nx_p, coordinates,
													num_par_max, num_par, num_created, rqm, rq_real, q_real, m_real,
										interpolation, ix, x, p, q, total_kinetic_4energy, total_number_particles_in_system_initially,
										boundry_conditions__type, boundry_conditions__thermal_type, 
										boundry_conditions__pth_bnd, boundry_conditions__pfl_bnd, total_xmoved, getCurrentState() );

}


void init_species_diagnostics_cpp( int *ndump_fac_ene, int *ndump_fac_temp,
                                      int *ndump_fac_raw, int *ndump_fac_tracks,
                                      int *phase_space__ndump_fac, int *phase_space__ndump_fac_tavg,
                                      int *phase_space__ndump_n_tavg, int *species_index, t_species_diagnostics& v) {
	//cout << "********************** ps out: " << *phase_space__ndump_fac << std::endl;
    //t_species_diagnostics& v = parts.species[(*species_index)-1].species_diagnostics;
    v.ndump_fac_ene = *ndump_fac_ene;
    v.ndump_fac_temp = *ndump_fac_temp;
    //v.n_tavg = *n_tavg;
    v.ndump_fac_raw = *ndump_fac_raw;
    v.ndump_fac_tracks = *ndump_fac_tracks;
    v.phase_space__ndump_fac= *phase_space__ndump_fac;
    v.phase_space__ndump_fac_tavg= *phase_space__ndump_fac_tavg;
    v.phase_space__ndump_n_tavg = *phase_space__ndump_n_tavg;
	
	//cout << "ndump_fac_ene:" << *ndump_fac_ene << std::endl;
	
}

extern "C" void init_species_diagnostics( int *ndump_fac_ene, int *ndump_fac_temp,
                                      int *ndump_fac_raw, int *ndump_fac_tracks,
                                      int *phase_space__ndump_fac, int *phase_space__ndump_fac_tavg,
                                      int *phase_space__ndump_n_tavg, int *species_index) {
	t_species_diagnostics& v = getCurrentState().parts.species[(*species_index)-1].species_diagnostics;
	init_species_diagnostics_cpp( ndump_fac_ene, ndump_fac_temp,
                                ndump_fac_raw, ndump_fac_tracks,
                                phase_space__ndump_fac, phase_space__ndump_fac_tavg,
                                phase_space__ndump_n_tavg, species_index, v);
}


void init_cuda_smooth(int*type,int*order,int*max_order,int*scoef,double*digital_A,double*digital_fc, t_smooth& smooth ) {
	smooth.max_order = *max_order;
	for(int i=0;i<p_x_dim;i++) {
		smooth.type[i] = type[i];
		smooth.order[i] = order[i];
		smooth.digital_A[i] = digital_A[i];
		smooth.digital_fc[i] = digital_fc[i];
	}
	int max_val_order = 0;
	for(int i=0;i<p_x_dim;i++) {
		if( smooth.order[i] > max_val_order ) max_val_order = smooth.order[i];
	}

	for(int i=0;i<p_x_dim;i++) {
		for(int j=0;j<(2*max_val_order+1);j++) {
			smooth.scoef[i].push_back(j);
		}
	}
}


void  setup_a_cuda_t_smooth(	int* ifsmooth, int* smooth_level, int* nsm, 
								float__fld* omega0, int* dir, float__fld* swfj,
								float__fld* scoef, t_smooth& v) {
/*
	v.enabled = 0; // Fortran has Flase == 0.. TRUE!= 0 (all values)
	for(int i=0;i<p_x_dim;i++) {
		v.ifsmooth[i] 		= ifsmooth[i];
		//cout << "I SEE enabled smooth " << i << " as " << ifsmooth[i] << std::endl;
		if(ifsmooth[i] != 0 ) v.enabled = 1;
	}

	if( v.enabled < 1) {
		// TODO: chek in other ways? messages?
		// The JAY we wwere linked is diabled.
		// Add function to our VDF (cal save the passed in one)
		//		so be able to callany of it's API
		return;
	}
	
	for(int i=0;i<p_x_dim;i++) {
		//v.ifsmooth[i] 		= ifsmooth[i];
		v.smooth_level[i] 	= smooth_level[i];
		v.nsm[i] 			= nsm[i];
	}
	v.omega0 = *omega0;
	v.dir = *dir;
	
	// the other two vaeiables are.. well.. variable... so let's figure out their size
	//	and get em.
	//int nscoef[ p_x_dim ];
	
	int max_value_in_nsm = -100;
	int max_value_nscoef = -100;
	
	for(int i=0;i<p_x_dim;i++) {
		int temp = v.nsm[i];
		if(temp > max_value_in_nsm) max_value_in_nsm = temp;
		
		int temp1 = 2*v.nsm[i] + 1;
		if(temp1 > max_value_nscoef) max_value_nscoef = temp1;
		v.nscoef[i] = temp1;
	}
	v.swfj__max_middle_idx = max_value_in_nsm;
	v.scoef__max_last_idx = max_value_nscoef;
	
	if(max_value_in_nsm  <= 0) {
		v.swfj__max_middle_idx = 1;
	}
	if(max_value_nscoef <= 0) {
		v.scoef__max_last_idx = 1;
	}
	
	// Having the max value for 'nsm' will allow us to get the
	//	mysterious size of swfj and then fill in from the passed in data.
	// the dimensions of swfj: (3, max_value_in_nsm, p_x_dim)
	//nt swfj_size = 3*max_value_in_nsm*p_x_dim;
	// I already cheated in allocation... middle exten is  ALWAYS 20.
	for(int i=0;i<p_x_dim;i++) {
		for(int j=0;j<v.nsm[i];j++) {
			v.swfj[i][j][0] =  swfj[3*i +20*j +0]; //swfj[i][j][0];
			v.swfj[i][j][1] =  swfj[3*i +20*j +1]; //swfj[i][j][0];
			v.swfj[i][j][2] =  swfj[3*i +20*j +2]; //swfj[i][j][0];
		}
	}
	
	for(int i=0;i<p_x_dim;i++) {
		v.nscoef[i] = 2*v.nsm[i]+1;
	}
	
	// for now, my stegy is repeat no setup calculations (with the exception
	//	of getting array size like we did above). I will just copy Osiris's smooth
	//	coeffecients/paramters in the short term and we can move calcualations in her..
	// further setup...
	
	// All we really need are filled in 'swfj', 'scoef' arrays.These are the centroids
	//
	*/
	/*
	for(int i=0;i<p_x_dim;i++) {
		if( !v.ifsmooth[i] ) {
			
			v.nsm[i]= 0;
			v.smooth_level[i]= 0;
			
		} else if( v.nsm[i] < 0) {
			v.ifsmooth[i] = 0; 
		}
	}*/
}

/*
	Setup function (not called directly from Fortran)
*/
void create_vdf( int *xdim, int * fdim, p_k_fld* f1, p_k_fld* f2, p_k_fld* f3, int * nx, int * gc_num, double * dx, t_vdf& v) {

	int i;
	int spatial_dims =  *xdim;
	v.x_dim = spatial_dims;
	v.f_dim = *fdim;
	
	// copy over the array paramters.
	//for(i=0;i<p_max_dim;i++) {
																							// LOG << "\t\tVdf config ("<<v.x_dim << " dim, " << v.f_dim << "field dim):" << std::endl;
	v.stride_x=0;v.stride_y=0;v.stride_z=0;
	for(i=0;i<p_max_dim;i++) {
		v.nx[i] =  0;
		v.gc_num[i][0] = 0;
		v.gc_num[i][1] = 0;
		v.dx[i] = 0;
		v.extent[i] = 0;
	}
	
	for(i=0;i<spatial_dims;i++) {
		v.nx[i] = nx[i];
		
		v.gc_num[i][0] = gc_num[i*2 + 0];
		v.gc_num[i][1] = gc_num[i*2 + 1];
		v.dx[i] = dx[i];
																						// LOG << "\t\tspace dim" << i << ":\tnx="<< nx[i] << "\tgc_num="<<gc_num[i*2 + 0]<<","<< gc_num[i*2 + 1]<<"\tdx="<< dx[i] << std::endl;
		v.extent[i] = nx[i] +  v.gc_num[i][p_lower] + v.gc_num[i][p_upper];
	}
	
	v.f1 = f1;
	v.f2 = f2;
	v.f3 = f3;
	//v.stride_x = 0; v.stride_y = 0; v.stride_z = 0;
	//v.stride_cells_x = 0; v.stride_y = 0; v.stride_z = 0;
	
	v.stride_y = v.f_dim;
	v.stride_x = v.f_dim*v.extent[0];
	v.stride_z = v.f_dim*v.extent[0]*v.extent[1];
	
	//v.stride_y = v.stride_x*v.extent[1];
	//v.stride_z = v.stride_y*v.extent[2];
	
	v.we_own_buffer = false;
	
	if(spatial_dims == 1) {
		if(v.f1 == NULL) {
			// make a buffer using the data we gots...
			v.f1 = new float__fld[v.buffer_size()]; //v.extent[0]];
			v.we_own_buffer = true;
		}
		v.field_data = v.f1;
		v.central_data = &( f1[ (v.gc_num[0][p_lower])*v.f_dim ] );
		// this is stupid but i am lazy.. it dont feel like fucking with c++ constructor ording.
		//Array<float, 2> all_data(v.field_data, shape( v.extent[0], p_f_dim), neverDeleteData);
		//v.f1_array = all_data;
		//v.f1_data = all_data(Range( v.gc_num[0][0],v.gc_num[0][0]+v.nx[0]));
	}
	if(spatial_dims == 2) {
	
		if(v.f2 == NULL) {
			// make a buffer using the data we gots...
			v.f2 = new float__fld[v.buffer_size()]; //v.extent[0]*v.extent[1]];
			v.we_own_buffer = true;
		}
		v.field_data = v.f2;
		
		// LOG << "\t\tSize of data buffer[y=" << (1-v.gc_num[1][0]) << ":" << (v.nx[1]+v.gc_num[1][1]) << "][x="<< (1-v.gc_num[0][0]) << ":" << (v.nx[0]+v.gc_num[0][1]) << "][fields_dims=" << v.f_dim << "]" << std::endl;
		// LOG << "\t\tSize of data buffer c index [y=" << (1-v.gc_num[1][0] - 1) << ":" << (v.nx[1]+v.gc_num[1][1]-1) << "][x="<< (1-v.gc_num[0][0]-1) << ":" << (v.nx[0]+v.gc_num[0][1]-1) << "][fields_dims=" << v.f_dim << "]" << std::endl;
		//// LOG << 
		//v.central_data = &( f2[(v.gc_num[0][p_lower]*v.stride_x) + (v.gc_num[1][p_lower])*v.stride_y] );
	 
		v.central_data = &( f2[ (v.gc_num[1][p_lower]*v.extent[0]+v.gc_num[0][p_lower])*v.f_dim ] );
		
		// this is stupid but i am lazy.. it dont feel like fucking with c++ constructor ording.
		//Array<p_k_fld, 3> all_data(v.field_data, shape( v.extent[1],  v.extent[0], p_f_dim), neverDeleteData);
		//Array<p_k_fld, 3> meat = all_data(	Range( v.gc_num[1][0],  v.gc_num[1][0]+v.nx[1]), 
		//								Range( v.gc_num[0][0],v.gc_num[0][0]+v.nx[0]));
										
		//cout << "---- " << all_data(0,0) << std::endl;
		//cout << "---- " << meat(0,0) << std::endl;
		
		//v.f2_array.reference(all_data);
		//v.f2_data.reference(meat);
										
										
		//cout << "---- " << v.f2_array(0,0) << std::endl;
	}
	if(spatial_dims == 3) {
		if(v.f3 == NULL) {
			// make a buffer using the data we gots...
			v.f3 = new float__fld[v.buffer_size()]; //v.extent[0]*v.extent[1]*v.extent[2]];
			v.we_own_buffer = true;
		}
		v.field_data = v.f3;
		v.central_data = &( f3[ ( v.gc_num[2][p_lower]*v.extent[1]*v.extent[0] + v.gc_num[1][p_lower]*v.extent[0]+v.gc_num[0][p_lower])*v.f_dim ] );
		//Array<p_k_fld, 4> all_data(v.field_data, shape(  v.extent[2], v.extent[1],  v.extent[0], p_f_dim), neverDeleteData);
		//v.f3_array = all_data;
		//v.f3_data = all_data (	Range( v.gc_num[2][0],  v.gc_num[2][0]+v.nx[2]), 
		//								Range( v.gc_num[1][0],  v.gc_num[1][0]+v.nx[1]), 
		//								Range( v.gc_num[0][0],v.gc_num[0][0]+v.nx[0]));
	}
	
	/*
	cout << "\t" << std::endl;
	cout << "\t spatial_dims" << spatial_dims<< std::endl;
	cout << "\t x_dim" << v.x_dim<< std::endl;
	cout << "\t v.f_dim" << v.f_dim<< std::endl;
	cout << "\t v.f1" << v.f1<< std::endl;
	cout << "\t v.f2" << v.f2<< std::endl;
	cout << "\t v.f3" << v.f3<< std::endl;
	cout << "\t v.field_data" << v.field_data<< std::endl;
	*/
}


void init_cuda_emf_cpp(	int *e_xdim, int * e_fdim, p_k_fld* e_f1, p_k_fld* e_f2, p_k_fld* e_f3, int * e_nx, int * e_gc_num, double * e_dx,
										int *b_xdim, int * b_fdim, p_k_fld* b_f1, p_k_fld* b_f2, p_k_fld* b_f3, int * b_nx, int * b_gc_num, double * b_dx,
										int *epart_xdim, int * epart_fdim, p_k_fld* epart_f1, p_k_fld* epart_f2, p_k_fld* epart_f3, int * epart_nx, int * epart_gc_num, double * epart_dx,
										int *bpart_xdim, int * bpart_fdim, p_k_fld* bpart_f1, p_k_fld* bpart_f2, p_k_fld* bpart_f3, int * bpart_nx, int * bpart_gc_num, double * bpart_dx,
										int *exte_xdim, int * exte_fdim, p_k_fld* exte_f1, p_k_fld* exte_f2, p_k_fld* exte_f3, int * exte_nx, int * exte_gc_num, double * exte_dx,
										int *extb_xdim, int * extb_fdim, p_k_fld* extb_f1, p_k_fld* extb_f2, p_k_fld* extb_f3, int * extb_nx, int * extb_gc_num, double * extb_dx,
										int solver, int coordinates, //p_k_fld* e_scratch_f1, p_k_fld* e_scratch_f2, p_k_fld* e_scratch_f3, p_k_fld* b_scratch_f1, p_k_fld* b_scratch_f2, p_k_fld* b_scratch_f3,
										int *ifextfld, int *emf_bound_types,
										int *smooth_type, int *smooth_niter, int *smooth_nmax,
										//int* ifsmooth, int* smooth_level, int* nsm, 
										//float__fld* omega0, int* dir, float__fld* swfj,
										//float__fld* scoef,
										t_emf &emf) { 
                    //double stencil_k1, double stencil_k2, double kark_k1, double kark_k2, double kark_k3,
										//p_k_fld marder_corr) {

		// LOG << "\t efield" << std::endl;
		create_vdf(e_xdim, e_fdim, e_f1, e_f2, e_f3, e_nx, e_gc_num, e_dx, emf.e);
		// LOG << "\t bfield" << std::endl;coordinates
		create_vdf(b_xdim, b_fdim, b_f1, b_f2, b_f3, b_nx, b_gc_num, b_dx, emf.b);
		// LOG << "\t e_part_field" << std::endl;
		create_vdf(epart_xdim, epart_fdim, epart_f1, epart_f2, epart_f3, epart_nx, epart_gc_num, epart_dx, emf.e_part);
		// LOG << "\t b_part_field" << std::endl;
		create_vdf(bpart_xdim, bpart_fdim, bpart_f1, bpart_f2, bpart_f3, bpart_nx, bpart_gc_num, bpart_dx, emf.b_part);
		// LOG << "\t e_part_field" << std::endl;
		create_vdf(exte_xdim, exte_fdim, exte_f1, exte_f2, exte_f3, exte_nx, exte_gc_num, exte_dx, emf.ext_e0);
		// LOG << "\t b_part_field" << std::endl;
		create_vdf(extb_xdim, extb_fdim, extb_f1, extb_f2, extb_f3, extb_nx, extb_gc_num, extb_dx, emf.ext_b0);

    // TODO: take out.
		//create_vdf(e_xdim, e_fdim, e_scratch_f1, e_scratch_f2, e_scratch_f3, e_nx, e_gc_num, e_dx, emf.e_scratch);
		//create_vdf(b_xdim, b_fdim, b_scratch_f1, b_scratch_f2, b_scratch_f3, b_nx, b_gc_num, b_dx, emf.b_scratch);
		//emf.e_scratch.zero_buffer();
	//	emf.b_scratch.zero_buffer();
		
		emf.solver = solver;
		if(*ifextfld==0) {
			emf.ifextfld = false;
		} else {
			emf.ifextfld = true;
		}
		
		// these mean nothing for now.. what we relly need is the type of external field (when we acaully genratte it in cuda)
		emf.is_ext_e = false;
		emf.is_ext_b = false;
		
		// setup the boundry conditions for the field...
		// loop over valid direction... aet the boundry type (upper and lower)
		for(int i=0;i<p_max_dim;i++) {
			emf.bnd_con.type[i][0] = emf_bound_types[2*i+0];
			emf.bnd_con.type[i][1] = emf_bound_types[2*i+1];
		}
		// furthur setup for lindman/vmpl will happen in subsequent setp function call.
		emf.smooth_type = *smooth_type;
		emf.smooth_niter = *smooth_niter;
		emf.smooth_type = *smooth_nmax;
		//setup_a_cuda_t_smooth( 	ifsmooth, smooth_level, nsm, omega0, dir, 
		//						swfj, scoef, ps.emf.smooth);
}

extern "C" void init_cuda_emf(	int *e_xdim, int * e_fdim, p_k_fld* e_f1, p_k_fld* e_f2, p_k_fld* e_f3, int * e_nx, int * e_gc_num, double * e_dx,
										int *b_xdim, int * b_fdim, p_k_fld* b_f1, p_k_fld* b_f2, p_k_fld* b_f3, int * b_nx, int * b_gc_num, double * b_dx,
										int *epart_xdim, int * epart_fdim, p_k_fld* epart_f1, p_k_fld* epart_f2, p_k_fld* epart_f3, int * epart_nx, int * epart_gc_num, double * epart_dx,
										int *bpart_xdim, int * bpart_fdim, p_k_fld* bpart_f1, p_k_fld* bpart_f2, p_k_fld* bpart_f3, int * bpart_nx, int * bpart_gc_num, double * bpart_dx,
										int *exte_xdim, int * exte_fdim, p_k_fld* exte_f1, p_k_fld* exte_f2, p_k_fld* exte_f3, int * exte_nx, int * exte_gc_num, double * exte_dx,
										int *extb_xdim, int * extb_fdim, p_k_fld* extb_f1, p_k_fld* extb_f2, p_k_fld* extb_f3, int * extb_nx, int * extb_gc_num, double * extb_dx,
										int solver, int coordinates,//, //p_k_fld* e_scratch_f1, p_k_fld* e_scratch_f2, p_k_fld* e_scratch_f3, p_k_fld* b_scratch_f1, p_k_fld* b_scratch_f2, p_k_fld* b_scratch_f3,
										int *ifextfld, int *emf_bound_types,
										int *smooth_type, int *smooth_niter, int *smooth_nmax,
										int*type,int*order,int*max_order,int*scoef,double*digital_A,double*digital_fc
										//int* ifsmooth, int* smooth_level, int* nsm, 
										//float__fld* omega0, int* dir, float__fld* swfj,
										//float__fld* scoef
										) { 
	
	t_emf& emf = getCurrentState().emf;
	init_cuda_emf_cpp(	e_xdim, e_fdim,  e_f1,  e_f2,  e_f3, e_nx, e_gc_num, e_dx,
										b_xdim, b_fdim,  b_f1,  b_f2,  b_f3, b_nx, b_gc_num, b_dx,
										epart_xdim, epart_fdim,  epart_f1,  epart_f2,  epart_f3, epart_nx, epart_gc_num, epart_dx,
										bpart_xdim, bpart_fdim,  bpart_f1,  bpart_f2,  bpart_f3, bpart_nx, bpart_gc_num, bpart_dx,
										exte_xdim, exte_fdim,  exte_f1,  exte_f2,  exte_f3, exte_nx, exte_gc_num, exte_dx,
										extb_xdim, extb_fdim,  extb_f1,  extb_f2,  extb_f3, extb_nx, extb_gc_num, extb_dx,
										solver, coordinates,  //e_scratch_f1, e_scratch_f2, e_scratch_f3, b_scratch_f1,b_scratch_f2,  b_scratch_f3,
										ifextfld, emf_bound_types,
										smooth_type, smooth_niter, smooth_nmax, emf);
										//ifsmooth, smooth_level, nsm, 
										//omega0, dir, swfj,
										//scoef, emf); 
t_smooth& emf_smooth = getCurrentState().emf.smooth;
init_cuda_smooth( type, order, max_order, scoef, digital_A, digital_fc, emf_smooth);

}

// TODO: this is badly named.. it is really setting up jay.. and it should be done when setting up t_particle structure.
void init_cuda_current_cpp( int *coordinates, int*gix_pos, int*bc_type, int*interpolation, t_current &jay) {
	//int* ifsmooth, int* smooth_level, int* nsm, float__fld* omega0, int* dir, float__fld* swfj, float__fld* scoef,
	jay.coordinates = *coordinates;
	jay.interpolation = *interpolation;

	int i=0;
	for(i=0;i<p_x_dim;i++) {
		jay.gix_pos[i] = gix_pos[i];
		jay.bc_type[i][0] = bc_type[2*1+0];
		jay.bc_type[i][1] = bc_type[2*i+1];
	}

	// create configuration data for smoothing of the Current.
	//setup_a_cuda_t_smooth( 	ifsmooth, smooth_level, nsm, omega0, dir, 
	//						swfj, scoef, jay.smooth);

}

extern "C" void init_cuda_current( int *coordinates, int*gix_pos, int*bc_type, int*interpolation,
	int*type,int*order,int*max_order,int*scoef,double*digital_A,double*digital_fc) {
																//int* ifsmooth, int* smooth_level, int* nsm, float__fld* omega0, int* dir, float__fld* swfj, float__fld* scoef) {
	t_current& jay =  getCurrentState().jay;
	t_smooth& jay_smooth = getCurrentState().jay.smooth;
	init_cuda_current_cpp( coordinates, gix_pos, bc_type, interpolation, jay );
	init_cuda_smooth( type, order, max_order, scoef, digital_A, digital_fc, jay_smooth);
}

// TODO: this is badly named.. it is really setting up jay.. and it should be done when setting up t_particle structure.
void init_cuda_current_vdf_cpp(	int *j_xdim, int * j_fdim, p_k_fld* j_f1, p_k_fld* j_f2, p_k_fld* j_f3, int * j_nx, int * j_gc_num, double * j_dx,
																std::vector<t_vdf>& pf){

	t_vdf s_obj;
	pf.push_back(s_obj);
	t_vdf& s = pf.back();

	create_vdf(j_xdim, j_fdim, j_f1, j_f2, j_f3, j_nx, j_gc_num, j_dx, s);
}

extern "C" void init_cuda_current_vdf(	int *j_xdim, int * j_fdim, p_k_fld* j_f1, p_k_fld* j_f2, p_k_fld* j_f3, int * j_nx, int * j_gc_num, double * j_dx) {
	t_current& jay =  getCurrentState().jay;
	init_cuda_current_vdf_cpp(	j_xdim, j_fdim, j_f1, j_f2, j_f3, j_nx, j_gc_num, j_dx, jay.pf);
}

// look for simple, known errors...
extern "C" void check_c_state( ) {
	if( getCurrentState().jay.pf.size() == 0) {
		std::cout << "No valid current Grid (i.e. VDF) objects were defined and mirrored in the C/C++ objects." << std::endl;
		exit(1);
	}

	if( getCurrentState().jay.pf.size() > 1) {
		std::cout << "More then 1 valid current Grid (i.e. VDF) objects were defined and mirrored in the C/C++ objects. The C/C++ code can only handle one at this time." << std::endl;
		exit(1);
	}
}



bool checkcuda(cudaError_t error, bool exit_on_failure) {
	if( cudaSuccess != error) {
		std::stringstream mess;
		mess << "   CUDA/GPU runtime error (generated by the CUDA driver)." << std::endl;
		mess << "   CUDA reports the error as: " << cudaGetErrorString(error)  << std::endl;

		FATAL_ERROR( mess.str());
	}
	return false;
}

bool checkcuda(cudaError_t error, std::string message, bool exit_on_failure) {
	if( cudaSuccess != error) {
		std::stringstream mess;
		mess << "   " << message << std::endl << std::endl;

		mess << "   The underlying cause of the error is a runtime CUDA/GPU runtime/driver error." << std::endl;
		mess << "   CUDA reports this underlying error as: " << cudaGetErrorString(error)  << std::endl;

		FATAL_ERROR( mess.str());
	}
	return false;
}

bool checkcuda(cudaError_t error, const char* file, int line, bool exit_on_failure) {
	if( cudaSuccess != error) {
		std::stringstream mess;
		mess << "   CUDA/GPU runtime error (generated by the CUDA driver)." << std::endl;
		mess << "   CUDA reports the error as: " << cudaGetErrorString(error)  << std::endl;

		FATAL_ERROR( mess.str(), file, line, exit_on_failure);
	}
	return false;
}

bool checkcuda(cudaError_t error, std::string message, const char* file, int line, bool exit_on_failure) {
	if( cudaSuccess != error) {
		std::stringstream mess;
		mess << "   " << message << std::endl << std::endl;

		mess << "   The underlying cause of the error is a runtime CUDA/GPU runtime/driver error." << std::endl;
		mess << "   CUDA reports this underlying error as: " << cudaGetErrorString(error)  << std::endl;

		FATAL_ERROR( mess.str(), file, line, exit_on_failure);
	}
	return false;
}

void FATAL_ERROR( std::string message, bool exit_on_failure ) {
	std::cout << "-------------------------------------------------------" << std::endl;
	std::cout << "A Fatal error occured in the C/C++ cuda code. Aborting." << std::endl;
	std::cout << "   Error information: " << std::endl;
	std::cout << "       MPI_NODE=" << CUDA_STATE.CUDA_CONFIG->node_id << std::endl;
	std::cout << "       HOSTNAME=" << CUDA_STATE.CUDA_CONFIG->mpi_hostname << std::endl;
	std::cout << "-------------------------------------------------------" << std::endl;
	std::cout << "The error message generated:" << std::endl;
	std::cout << "" 								<< std::endl;
	std::cout << message << std::endl;

	if(exit_on_failure) exit(-1);
}

void FATAL_ERROR( std::string message, const char* filename, int line_number, bool exit_on_failure) {
	std::cout << "-------------------------------------------------------" 					<< std::endl;
	std::cout << "A Fatal error occured in the C/C++ cuda code. Aborting." 					<< std::endl;
	std::cout << "   Error information: " 																					<< std::endl;
	std::cout << "       MPI_NODE=" << CUDA_STATE.CUDA_CONFIG->node_id 							<< std::endl;
	std::cout << "       HOSTNAME=" << CUDA_STATE.CUDA_CONFIG->mpi_hostname 				<< std::endl;
	std::cout << "       File Name=" << filename <<  " (line " << line_number << ")"  << std::endl;
	std::cout << "-------------------------------------------------------" 					<< std::endl;
	std::cout << "The error message generated:" << std::endl;
	std::cout << "" 														<< std::endl;
	std::cout << message 												<< std::endl;

	if(exit_on_failure) exit(-1);
}
void FATAL_ERROR( std::string message) {
	FATAL_ERROR(message, true);
}
void FATAL_ERROR( std::string message, const char* filename, int line_number) {
	FATAL_ERROR(message, filename, line_number, true);
}
