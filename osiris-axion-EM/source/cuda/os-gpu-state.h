// contents of
// ViktorCodeContext
// ViktorGlobalState


/*
	This is for the first step of integrating Osiris with Viktor's code.
*/
	struct ViktorStuff {
		
		// holds the edges (in node number)
		float__part * edges;
		
		// bulk buffers...
		float__part * part;
		float__part * part2;
		
		float__part * part_adam_temp1;
		float__part * ppart_adam_temp1;
		float__part * part_adam_temp2;
		float__part * ppart_adam_temp2;
		int * mpic_adam_temp1;
		int * mpic_adam_temp2;
		
		int * kpic_adam_temp1;
		int * kpic_adam_temp2;
		int * jpic_adam_temp1;
		int * jpic_adam_temp2;

		float__part * osiris_q2;
		int * osiris_ix2;
		//float__part * qe;
		//float__part * fxye;
		float__part * qt;
		int * ihole;
		//int * npic;
		int * mpic;
		float__part * sbufl_x;
		float__part * sbufr_x;
		float__part * rbufl_x;
		float__part * rbufr_x;
		
		float__part * sbufl_y;
		float__part * sbufr_y;
		float__part * rbufl_y;
		float__part * rbufr_y;
	
		float__part * sbufl_z;
		float__part * sbufr_z;
		float__part * rbufl_z;
		float__part * rbufr_z;
		
		int* sort__z__scratch1;
		int* sort__z__scratch2;
		int* sort__z__scratch3;
		int* sort__z__scratch4;
		
		//float__part * scr;
		
		float__part * ppart;
		float__part * ppart2;
		float__part * ppart3;
		float__part * ppbuff;
		int * kpic;
		int * jpic;
		int * ncl;
		int * iholep;
		int * ncll_y;
		int * nclr_y;
		int * mcll_y;
		int * mclr_y;
		
		int * ncll_x;
		int * nclr_x;
		int * mcll_x;
		int * mclr_x;

		int * ncll_z;
		int * nclr_z;
		int * mcll_z;
		int * mclr_z;
		
		int * kpic_debugging;
		
		// debugging..
		float__part * q2;
		int * lpic;
		int * ix__debug;
		float__part * x__debug;
		float__part * p__debug;
		float__part * q__debug;
		
		
		
	};
	
	struct gpu_buffers {
		//float__part * g_qe;
		//float__part * g_fxye;
		//float__part * g_fxyb;
		float__part * g_wke;
		
		float__part * g_ppart;
		//float__part * g_ppart_output;
		float__part * g_ppbuff;
		int * g_kpic;
		int * g_kpic_debug;
		
		int * g_mpic;
		float__part * g_sbufl_x;
		float__part * g_sbufr_x;
		float__part * g_sbufl_y;
		float__part * g_sbufr_y;
		float__part * g_sbufl_z;
		float__part * g_sbufr_z;
		
		int * g_ncl;
		int * g_ihole;
		int * g_ihole_fake;
		int * g_ncll_x;
		int * g_nclr_x;
		int * g_ncll_y;
		int * g_nclr_y;
		int * g_ncll_z;
		int * g_nclr_z;
		int * g_irc;
		int * fake_g_irc;
		
		// Cancer
		double* g_cancer_KE_temp;
		double* g_cancer_Q1_temp;
		double* g_cancer_Q2_temp;
		double* g_cancer_Q3_temp;
	};

// variables shared across species...
// this is a temp class used in development.
//class CudaOsirisSimulationState
class ViktorGlobalState {
	public:
		ParticleBuffer g_qe;
		ParticleBuffer g_qe_fake;  // for debugging TODO: take out
		
		//ParticleBuffer g_qe_terrible;
		//ParticleBuffer g_qe_terrible_parts;
		
		ParticleBuffer g_fxye;
		ParticleBuffer g_fxyb;
		
		ParticleBuffer g_fxye_ext;
		ParticleBuffer g_fxyb_ext;
		
		ParticleBuffer g_fxye__for_particle_push;
		ParticleBuffer g_fxyb__for_particle_push;
		
		cuda_pusher_function __pusher;
		cuda_deposit_function __deposit;
		
    void init_viktors_global_state();
		void init_viktors_global_state__moving_window();
		void destory_and_cleanup_viktors_global_state( );
	
	
    //StopWatchTimer* loop_timer;
    
    //vector<ScopeTimer*> fortran_timers;
    //vector<string> fortran_timer_names;
    
	
	// variables for moving boundrys
	bool is_moving_window_used;
	int moving_window__part_tile_size;
	int moving_window_inject__num_tiles;
	int *moving_window_inject__mpic_host;
	float__part *moving_window_inject__buffer_host;
	int *moving_window_inject__mpic_gpu;
	float__part *moving_window_inject__buffer_gpu;
	
	
	int* g_moving_windows__circular_tile_buffer_map;
	int* host_moving_windows__circular_tile_buffer_map1;
	int* host_moving_windows__circular_tile_buffer_map2;
	int *host_moving_windows__circular_tile_buffer_map__current;
	int *host_moving_windows__circular_tile_buffer_map__scratch;
	
	
	
    //vector<StopWatchTimer*> loop_timers;
   /* 
    StopWatchTimer* loop_timer1;
    StopWatchTimer* loop_timer2;
    StopWatchTimer* loop_timer3;
    StopWatchTimer* loop_timer4;
    StopWatchTimer* loop_timer5;*/
    
};

class ViktorCodeContext {
	public:
		ViktorCodeContext() : __counter(0), config__use_this_nppmx( 0 ) { };
			//__counter = 0;
		//};
		
		//~ViktorCodeContext();
		
		// TODO: remove..
		
		int __counter;		// TODO: beyond logghetto.. tkaae this out.
		int species_index;
		int config__use_this_nppmx;		// if set to <= 0 use this number.. if not then calculate
		
		// TODO move this..
		void* thermal_boundrires__random_number_device_pointers[3][2];
		int is_thermal_bath_bndry_used;
		int particle_boundry_type[p_x_dim][2];
		
		// dervied state variables for moving/absorbing boundries
		int any_correcting_boundry_deposits_needed;
		int bndry_discard__lower_x;	
		int bndry_discard__upper_x;
		int bndry_discard__lower_y;	
		int bndry_discard__upper_y;
		int bndry_discard__lower_z;	
		int bndry_discard__upper_z;
		
		float__part qme;
		int m1;
		float__part epsmax, eps;
		int ntime;
		int idimp;
		int idps;
		int mx;
		int my;
		int mz;
		int nblock;
		int nscache;
		//double np;
		int nx;
		int ny;
		int nz;
		//int nxe;
		//int nye;
		//int ny1;
		int num_cells__this_node;
		int num_cells__system;
		int num_tiles__this_node;
		int num_tiles__system;
		int dim_tiles_x__system;
		int dim_tiles_y__system;
		int dim_tiles_z__system;
		int mx1;
		//int nvp;
		int this_node__node_idx_along_x;
		int this_node__node_idx_along_y;
		int this_node__node_idx_along_z;
		int number_nodes_x;
		int number_nodes_y;
		int number_nodes_z;
		int total_number_nodes;
		int idproc;
		int kstrt;
		int nxp;
		int nyp;
		int nzp;
		int noff;  // badly names.. it's this node's starting global gird location along x
		int moff;  // badly names.. it's this node's starting global gird location along x
		int zoff;
		int nypmx;
		int nypmn;
		int npmax;
		int nbmax;
		int ntmax;
		int mxp1;
		int myp1;
		int mzp1;
		//int mxyp1;
		int npp;
		int nps;
		int nppmx;
		int ircc;
		int irccc;
		int ircccc;
		int irc;
		int irc_buffer[6];
		int irccccc;
		int ircccccc;
		int nppmx0;
		int ntmaxp;
		int ntmaxp_fake;
		int npbmx;
		int ierr;
		float__part qbme;
		float__part affp;
		int ipbc;
		float__part wke;
		float__part we;
		float__part wt;
		float__part gepsmax;
		float__part gqepsmax;
		int info[7];
		int sortime;
		//int ndim;
		
		int nppmx_initial;
		
		ParticleBuffer ppart;
		ParticleBuffer ppart2;
		ParticleBuffer ppart3;
		ParticleBuffer ppbuff;
		ParticleBuffer part;
		ParticleBuffer part2;
		
		ParticleBuffer part_adam_temp1;
		ParticleBuffer ppart_adam_temp1;
		ParticleBuffer part_adam_temp2;
		ParticleBuffer ppart_adam_temp2;
		
		//ParticleBuffer g_qe;
		//ParticleBuffer g_fxye;
		//ParticleBuffer g_fxyb;
		ParticleBuffer g_wke;
		
		ParticleBuffer g_ppart;
		//ParticleBuffer g_ppart_output;
		ParticleBuffer g_ppbuff;
		ParticleBuffer g_sbufl_x;
		ParticleBuffer g_sbufr_x;
		ParticleBuffer g_sbufl_y;
		ParticleBuffer g_sbufr_y;
		ParticleBuffer g_sbufl_z;
		ParticleBuffer g_sbufr_z;
		
		
		//float__part dt;
		ViktorStuff v_stuff;
		gpu_buffers v_gpu_stuff;

#if 0		
		// but the things in VikterStuff and gpu_buffer in namespace as ParticleBuffer Objects...
		
		void particle_sort_2d__attempt_recovery_new(int n);
		void print_particle_buffer_stats();
		int particle_sort_2d__cuda__part1__attempt_recovery(int n, bool do_sorting, bool is_retry_attempt, int* force_mpi_resize_mpi_buffers );
		int particle_sort_2d__cuda__part2__attempt_recovery(int n, bool do_sorting, bool is_retry_attempt );
		
		void validate_and_log_sort_1d__1_way_decomp(int n);
		void validate_and_log_sort_2d__2_way_decomp(int n);
		void validate_and_log_sort_3d__3_way_decomp(int n);
		
		
		void destory_and_cleanup_viktors_code__species( int species_index, bool keep_mpi_cpu_xfer_array, bool opt_out_of_reset );
		void init_viktors_code__species(int species_index, bool do_not_calc_node_configuration, bool suppress_output, bool keep_mpi_cpu_xfer_arrays, bool opt_out_of_reset);
		void reallocate_support_buffers(  int species_index);
		void reset_gpu_sort_error_codes( );
		//void refresh_GPU_particle_data_from_host();
		void copy_osriris_particle_buffer_to_gpu_effecient (	ParticleBuffer& target_gpu_buffer,  int* mpic,
																//int * target_gpu_tile_size_buffer,
																int* ix_source, float__part* x_source,
																float__part* p_source, 
																float__part* q_source,
																vector<int>& particle_ids);
				
		 void copy_osriris_particle_buffer_to_gpu(	ParticleBuffer& target_gpu_buffer, 
													int * target_gpu_tile_size_buffer,
													int* ix_source, float__part* x_source,
													float__part* p_source, 
													float__part* q_source,
													vector<int>& particle_ids);
	
															
		void copy_particles_from_gpu_to_osiris_buffer(	ParticleBuffer& source_gpu_buffer,
														int * target_gpu_tile_size_buffer,
														int* ix_target, 
														float__part* x_target, 
														float__part* p_target, 
														float__part* q_target,
														vector<int>& particle_ids);
                            
		int copy_particles_from_gpu_to_osiris_buffer_effecient(ParticleBuffer& g_ppart,
															int * g_kpic,
															int* ix, 
															float__part* x, 
															float__part* p,
															float__part* q,
															int target_buffer_size,
															int species_index,
															vector<int> & particle_ids ) ;
		void copy_osriris_particle_buffer_to_gpu_effecient__no_sort(ParticleBuffer& ppart,int* gpu_kpic,
															int* ix, float__part* x, 
															float__part* p, 
															float__part* q,
															vector<int> & particle_ids_target);
		
void print_particle( int tile_idx, int par_idx) {
			print_a_particle( v_gpu_stuff.g_ppart,
							  v_gpu_stuff.g_kpic,
							  nppmx0,
							  idimp,
							tile_idx,
							par_idx);
}	
	
		void map_linear_particle_to_tile_adr(int linear_index, int * tile_index, int *offset_within_tile) {
			//makes too many assumptions...
			//the standard mpic buffer had the addresing info we need.
			int linear_tile_index = -1;
			for(int i=1;i<num_tiles__this_node;i++) {
				printf( "%d:%d   \n", v_stuff.mpic[i], v_stuff.mpic[i]);
				if(linear_index < v_stuff.mpic[i]) {
					linear_tile_index = i-1;
					break;
				}
			}
			
			if(linear_tile_index < 0) {
				printf("ERROR. COULD NOT MAP LINEAR INDEX %d \n\n", linear_index);
				exit(-1);
			}
			*tile_index = linear_tile_index+1;
			*offset_within_tile = (linear_index - v_stuff.mpic[linear_tile_index]);
			
		};
		
		
		//void debug__print_state_1();
		//void debug__print_state_2(); 
		void __ADVANCE_DEPOSIT_2D_CUDA__CPP__(float__part dt, float__part t, int n, int ndump, int calc_total_species_KE4);
		// void current_deposit_2d_cuda();	
		void particle_push_2d_cuda(float__part dt, int n, int calc_total_species_KE4);
		
		void particle_sort_2d__attempt_recovery(int n);
		int particle_sort_2d_cuda__attempt_recovery(int n);
		void particle_sort_2d_cuda(	int particle_boundry_type__x_lower, int particle_boundry_type__x_upper,
									int particle_boundry_type__y_lower, int particle_boundry_type__y_upper );
		//void particle_sort_2d_cpu();
		
		// TODO: take this out.. it is only for low level deveolpment debugging...
		//			should be replaced by particle tracking.
		void debug_track_particle_gpu(float__part* values);
		void test_particle_assertions(int line_number, const char* message_passed, const char* message_failed);
		void test_particle_assertion__filled_below_ni(int is_k_updated_yet);
		void test_particle_assertions__pusher_output(int line_number, const char* message_passed, const char* message_failed);
		void test_particle_assertions__pre_pusher(int line_number, const char* message_passed, const char* message_failed);
		void test_particle_assertion__filled_above_ni();
		void test_particle_assertions__after_species_bndry( );
		void debug__particle_tracking(int n);
		void debug__particle_tracking__fields(int n);
		
		// PER PARTICLE ATTRIBURE DEBUGGING BUFFERS...
		//PPA_ATTRIBUTE_BUFFERS ppa_buffs;
		
		/*
		float__part *ppa__e1;
		float__part *g_ppa__e1;
		float__part *ppa__e2;
		float__part *g_ppa__e2;
		float__part *ppa__e3;
		float__part *g_ppa__e3;
		float__part *ppa__b1;
		float__part *g_ppa__b1;
		float__part *ppa__b2;
		float__part *g_ppa__b2;
		float__part *ppa__b3;
		float__part *g_ppa__b3;
		*/

#endif

};


