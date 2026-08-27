
#ifndef OS_CUDA_STATE_H
#define OS_CUDA_STATE_H


#include "os-cuda-input-deck.h"

#define float__part double
#define float__fld double
#define p_x_dim 2
#define p_p_dim 2
#define p_f_dim 2

#define p_max_dim 3
#define p_k_fld double

#define p_lower 0
#define p_upper 1

// throw-away tyedefs for functions that simulate c multidimensioal arrays that directly point to the live OSIRIS data structure..
// TODO: use ISO-C binding in Fortran 2003 for tighter integration.
typedef double __fixed_sizesmal_aray__live_osiris__crap1[p_max_dim][2] ;
typedef int __fixed_sizesmal_aray__live_osiris__crap2[p_max_dim][2];
typedef bool __fixed_sizesmal_aray__live_osiris__crap3[p_max_dim];
typedef double __fixed_sizesmal_aray__live_osiris__crap4[p_max_dim];
typedef int __fixed_sizesmal_aray__live_osiris__crap5[p_max_dim];
typedef double __fixed_sizesmal_aray__live_osiris__crap6[p_x_dim];


//using namespace std;

/*
	Osiris Mirroring state.
*/
class  t_node_conf {
	public:
		//t_node_conf();
		//~t_node_conf();
		int x_dim;								// number of dimensions..
		//integer :: comm = mpi_comm_null
		//logical :: free_comm = .true.
		int nx[p_max_dim];						// number of nodes in each dimension
		int no_num;								// total number of nodes
		int my_ngp[p_max_dim];					// My Node Grid Position (within number given  by nx)
		int neighbors[p_max_dim][2];			// my neighbors.
		int my_aid;								// array identification number of this node
		int my_MPI_rank;
		int *ngp;								// The 'my_ngp' for all processes... A map
		int mpi_comm;							// Mpi cummuncation ordinal
		bool ifpr[p_max_dim];					// should peroid boubry conditions be used along each direction?
};

class t_options	{
	public:
		//t_options();
		//~t_options();
		
		// members.. all public for now.
		bool test;
		bool restart;
		char* input_file;
		char* work_dir;
		int random_seed;
		int algorithm;
		double n0;
		double gamma;
};
class t_time {
	public:
		//t_time();
		//~t_time();

		double tmin;
		double tmax;
		double t;
};
class t_time_step {
	public:
		//t_time_step();
		//~t_time_step();
};
class t_antenna {
	public:
		//t_antenna();
		//~t_antenna();

		p_k_fld 	a0;
		p_k_fld 	t_rise;
		p_k_fld 	t_fall;
		p_k_fld 	t_flat;
		p_k_fld 	omega0;
		p_k_fld		pol;
};
class t_antenna_array {
	public:
		//t_antenna_array();
		//~t_antenna_array();
		int n_antenna;
		t_antenna* ant_array;
};

	class t_space {
	public:
		//t_space();
		//~t_space();
		int x_dim;								// number of dimensions..
		
		// Note:
		// 	Becuause of moving window, we will not copy these values.. rather
		//		we will make the member variables references.. This is bad c++
		//		because Orisiris might kill the original object leaving dangling refernece
		//		Deal with it later when the final integration of Osiris and c++ is thought out.
		
		// These are the the orginal (non-reference) members
		//double x_bnd[p_max_dim][2];					// Tells use the sim coordinates of the starting and ending cell of this node.
		//double x_bnd_old[p_max_dim][2];
		//double move_num[p_max_dim][2];			// for moving window....
		//double xmoved[p_max_dim][2];				// for moving window....	
		//int total_i_moved[p_max_dim];
		//bool if_move[p_max_dim];
		//double move_vel;
		//double x_bnd_initial[p_max_dim][2];
		
		// MOVING WINDOW DIRECT OSIRIS LINKED VARIABLES
		double (*_x_bnd)					[p_max_dim][2];
		double (*_x_bnd_old)			[p_max_dim][2];
		int    (*_move_num)				[p_max_dim][2];
		double (*_xmoved)					[p_max_dim][2];
		int    (*_total_i_moved)	[p_max_dim][2];
		int    (*_if_move)				[p_max_dim];
		double (*_mov_vel)				[p_max_dim];
		double (*_x_bnd_initial)	[p_max_dim][2];
		
		// a bit of a hack... just to make pointers to osiris objetts semi-tranparent

		
		__fixed_sizesmal_aray__live_osiris__crap1& x_bnd_() { return *_x_bnd; }
		__fixed_sizesmal_aray__live_osiris__crap1& x_bnd_old_() { return *_x_bnd_old; }
		__fixed_sizesmal_aray__live_osiris__crap2& move_num_() { return *_move_num; }
		__fixed_sizesmal_aray__live_osiris__crap1& xmoved_() { return *_xmoved; }
		__fixed_sizesmal_aray__live_osiris__crap2& total_i_moved_() { return *_total_i_moved; }
		__fixed_sizesmal_aray__live_osiris__crap5& if_move_() { return *_if_move; }
		// TODO: velovity is is just one double.. not an array
		__fixed_sizesmal_aray__live_osiris__crap4& mov_vel_() { return *_mov_vel; }
		__fixed_sizesmal_aray__live_osiris__crap1& x_bnd_initial_() { return *_x_bnd_initial; }
				
};

	class  t_grid {
		public:
			//t_grid();
			//~t_grid();

			int x_dim;
			int lb_type;
			int dynamic;
			float__part max_unbalance;
			float__part cell_weight;
			int nnodes[p_max_dim];					// number of nodes in each direction
			
			int g_nx[p_max_dim];					// global number of cells  on a node in each direction
			
			int min_nx[p_max_dim];					// minimum number of cells  on a node in each direction 
													//	(used by load balancing)
			int max_nx[p_max_dim];					// maximum number of cells on a node in each direction 
													//  (used by load balancing)
			
			int coordinates;						// type of coordinate system to use.

			std::vector<std::vector <int> > x1_nx;
			std::vector<std::vector <int> > x2_nx;
			std::vector<std::vector <int> > x3_nx;
			//int *x1_nx; 							// number of cells for each partition (along x1)
			//int *x2_nx;							// number of cells for each partition (along x1)
			//int *x3_nx; 							// number of cells for each partition (along x1)
};

class t_spe_bound {
	public:
		//t_spe_bound();
		//~t_spe_bound();
		
		int type [p_max_dim][2];
		int thermal_type[p_max_dim][2];
		float__part pth_bnd[p_p_dim][2][3];
		float__part pfl_bnd[p_p_dim][2][3];
		
};

class t_species_diagnostics {
  public:
    int ndump_fac_ene;
    int ndump_fac_temp;
    int ndump_fac_raw;
    int ndump_fac_tracks;
    int phase_space__ndump_fac;
    int phase_space__ndump_fac_tavg;
    int phase_space__ndump_n_tavg;
};


class t_species {
	public:
		//t_species();
		//~t_species();
		 
		int sp_id;									// this is a species id
		int num_par_x[p_x_dim];						// this is the Particles per cell
		double dx[p_max_dim];						// The cell size (1,2,3 dims]
		int my_nx_p[p_max_dim][3];					// which cells a node is repsonsible for.. [x][0] = start of dim_x [x][1]=end of dim_x [x][2]=number cel dim_x on this node.
		int coordinates;							// type of coordinates
		int num_par_max;							// This is the total size of the particles buffer.
		int* num_par_max__ptr;
		int number_of_particles;					// total number of particles we care about.
		integer64 total_number_particles_in_system_initially; // the intial total number fo particales!
		int num_created;							// total number of particles created.
		float__part rqm;
		float__part rq_real;
		float__part q_real;
		float__part m_real;
		bool is_energy;
		int pos_type;				// TODO: fill this in.
		
		t_spe_bound bnd_con;
		
		// MOVING WINDOW DIRECT OSIRIS LINKED VARIABLES
		double (*_g_box) [p_max_dim][2];				// box size.  [dim][0] = leftmost_cord alon dim.. [dim][1] = right most coord along dim 
		double (*_total_xmoved)[p_x_dim];
		
		// energy always accumlated as a double.
		double *total_kinetic_4energy; //energy[4];			// TODO: init this correctly
		
		float__part *x;
		float__part *p;
		float__part *q;
		int *ix;
		
		
		float__part *x_scratch;
		float__part *p_scratch;
		float__part *q_scratch;
		int *ix_scratch;
		
		
		float__part *x_transpose;
		float__part *p_transpose;
		float__part *q_transpose;
		int *ix_transpose;
		
		int interpolation;
		
		float__part* cudaX;
		float__part* cudaP;
		float__part* cudaQ;
		int*	cudaiX;
		
		// todo: probably should remove vectors in final version and
		//			synch better with Osirirs's buffer sizes. but maybe not.
		std::vector<int> particle_ids;
		std::vector<int> particle_ids_scratch;
		
	  // Place for the pointers it turns out I actually need.
	  int * num_par;
	  
		//Array<float__part, 2> nose;
		
		//vector<int> particle_keys_scratch;
		//vector<int> particle_keys;
		
		std::vector<std::pair< int, int> > particle_keys;
		std::vector<std::pair< int, int> > particle_keys_scratch;
		t_species_diagnostics species_diagnostics;
	  
		
		__fixed_sizesmal_aray__live_osiris__crap1& g_box_() { return *_g_box;}
		__fixed_sizesmal_aray__live_osiris__crap6& total_xmoved_() { return *_total_xmoved; }
	};

	class t_particles {
		public:
			//t_particles();
			//~t_particles();
			int num_species;
			//t_species_c* species;
			int num_hybrid;
			//t_hybrid_c* hybrid;
			
			// try using STL
			std::vector<t_species> species;
			//vector<t_hybrid> hybrid;
	};

/*
	notes about VDF (i.e. grids):
		'central data' means non-guard cell data.

		in fortran array syntax a vdf.f2 is:  vdf.f1[1:3][1:nx]							(not including guard cells)
		in fortran array syntax a vdf.f2 is:  vdf.f2[1:3][1:nx][1:ny] 			(not including guard cells)
		in fortran array syntax a vdf.f3 is:  vdf.f3[1:3][1:nx][1:ny][1:nz] (not including guard cells)
		
		including guard cells:
				vdf.f1:  f2[ 1:3, 1-gc_num[p_lower,1]:nx+gc_num[p_upper,1] ]
				vdf.f2:  f2[ 1:3, 1-gc_num[p_lower,1]:nx+gc_num[p_upper,1] ][ 1-gc_num[p_lower,2]:ny+gc_num[p_upper,2]]
				vdf.f3:  f2[ 1:3, 1-gc_num[p_lower,1]:nx+gc_num[p_upper,1] ][ 1-gc_num[p_lower,2]:ny+gc_num[p_upper,2]][ 1-gc_num[p_lower,3]:nz+gc_num[p_upper,3]]

			In c-arrays:
					vdf.f1[0:nx-1][0:2], vdf.f2[0:ny-1][0:nx-1][0:2], vdf.f3[0:nz-1][0:ny-1][0:nx-1][0:2]   (not including guard cells)
								
			including guard cells (x same for other dims):
				lower_bound = 1-gc_num[0][p_lower]-1 = -gc_num[0][p_lower]
					So using only positive array indices: central data starts at index = gc_num[0][p_lower]
				upper bound = nx + gc_num[0][p_uppper]-1 + gc_num[0][p_lower] = (nx-1) + gc_num[0][p_lower] + gc_num[0][p_uppper]
					So using only positve array indices	: Central data ends at index = (nx-1) + gc_num[0][p_lower]
																							: Final index of array is index = (nx-1) + gc_num[0][p_lower]+gc_num[0][p_upper]
				Using shiftd array where index=0 is at first 'Central Data' element:
						lower bound:  -gc_num[0][p_lower]
						upper bound:  nx+gc_num[0][p_upper] - 1

				For a 2D array, the first Central Data element is (using onlypositive indicies):
					gc_num[0][p_lower]+gc_num[1][p_lower]*(gc_num[0][p_lower]+gc_num[0][p_upper]+nx)		

		  // for a FROTRAN loop that looks like:
		  do y=lbound(jay%f2,3) ,ubound(jay%f2,3)
			do x=lbound(jay%f2,2) ,ubound(jay%f2,2)
				jay%f2(:,x,y) = float(x-1) + float(y-1)/100.0
			end do
		  end do 
			  
		for(int i=0;i<4;i++) {
			int base = (i)* ps.jay.pf.extent[0]*p_f_dim;
			printf("\t\t %f,%f,%f ",   ps.jay.pf.central_data[base+0*p_f_dim+0], ps.jay.pf.central_data[base+0*p_f_dim+1], ps.jay.pf.central_data[+base+0*p_f_dim+2]); 
			printf("\t\t %f,%f,%f \n", ps.jay.pf.central_data[base+1*p_f_dim+0], ps.jay.pf.central_data[base+1*p_f_dim+1], ps.jay.pf.central_data[base+1*p_f_dim+2]); 
		}
		will print    jay[y][x][comp]:			j[0][0][0]	j[0][0][1] j[0][0][2] j[0][1][0] j[0][1][1] j[0][1][2] ...
												j[1][0][0]	j[1][0][1] j[1][0][2] j[1][1][0] j[0][1][1] j[1][1][2] ...
												j[2][0][0]	j[2][0][1] j[2][0][2] j[2][1][0] j[2][1][1] j[2][1][2] ...

												0.0		0.0		0.0		1.0		1.0		1.0 	...
												0.01	0.01 	0.01	1.01	1.01	1.01	...
												0.02	0.02 	0.02	1.02	1.02	1.02	...
												0.03	0.03 	0.03	1.03	1.03	1.03	...
					This is in clostest locality mem access.
					
		// using Blitz++
		Array<float, 3> all_data(ps.jay.pf.field_data, shape( ps.jay.pf.extent[1],  ps.jay.pf.extent[0], p_f_dim), neverDeleteData);
		Array<float, 3> central_data = nose_all(Range( ps.jay.pf.gc_num[1][0],  ps.jay.pf.gc_num[1][0]+ps.jay.pf.nx[1]), Range( ps.jay.pf.gc_num[0][0],ps.jay.pf.gc_num[0][0]+ps.jay.pf.nx[0]));
	*/
		// nx[0,1,2]: number of data cells in each dimension
		// gc_num: number guard cells
		// extent, nx_with_gc: nx+gc_num_lower+gc_num_upper
		// stride: usually same as extent but not always.. the number needed to goto next index for a dimension. Helps with readablilty.
		class t_vdf {
			
			public:
				t_vdf() {};
				
				t_vdf(int *xdim, int * fdim, p_k_fld* f1, p_k_fld* f2, p_k_fld* f3, int * nx, int * gc_num, double * dx) {
					
				}
				
				//~t_vdf();

				int x_dim;
				int f_dim;
				int nx[p_max_dim];
				int gc_num[p_max_dim][2];
				double dx[p_max_dim];
				int extent[p_max_dim];
				int nx_with_gc[p_max_dim];

				float__fld*		f1;
				float__fld*		f2;
				float__fld*		f3;
				float__fld*		f1_transpose;
				float__fld*		field_data;
				
				// points to where the central data lives.
				float__fld*		central_data;
				
				// strides
				int				stride_x;
				int				stride_y;
				int				stride_z;	
				
				int 			length_buffer_x;
				int 			length_buffer_y;
				int 			length_buffer_z;
				
				//Array<float__fld, 2> f1_array;
				//Array<float__fld, 2> f1_data;
				
				//Array<float__fld, 3> f2_array;
				//Array<float__fld, 3> //;
				
				//Array<float__fld, 4> f3_array;
				//Array<float__fld, 4> f3_data;
				
				// TDO: smart pointer.
				bool we_own_buffer;
				
				void zero_buffer() {
				
					int size = 0;
					float__fld * buffer = 0;
					
					if(p_x_dim == 1) {
						size = f_dim*extent[0];
						buffer = f1;
					} else if(p_x_dim == 2) {
						size = f_dim*extent[0]*extent[1];
						buffer = f2;
					} else {
						size = f_dim*extent[0]*extent[1]*extent[2];
						buffer = f3;
					}
					
					for(int i=0;i<size;i++) {
						buffer[i] = 0.0;
					}
				}
				
				void eye_check_pattern() {
					int size = 0;
					float__fld * buffer = 0;
					
					if(p_x_dim == 1) {
						size = f_dim*extent[0];
						buffer = f1;
					} else if(p_x_dim == 2) {
						for(int j=0;j<extent[1];j++) {
							for(int i=0;i<extent[0];i++) {
								for(int n=0;n<f_dim;n++) {
									f2[(j*extent[0] + i)*f_dim + n] = (float__part) i + (float__part) (j/100.0);
								}
							}
						}
					} else {
						for(int k=0;k<extent[2];k++) {
							for(int j=0;k<extent[1];j++) {
								for(int i=0;k<extent[0];i++) {
									f2[(j*extent[0] + i)*f_dim + i] = (float__part) i + (float__part) (j/100.0);
								}
							}
						}
						size = f_dim*extent[0]*extent[1]*extent[2];
						buffer = f3;
					}
					
					for(int i=0;i<size;i++) {
						buffer[i] = 0.0;
					}		
				}
				
				// TODO: makes these constants (templates) dso no lookup overhead.
				int get_central_data_flat_index() {
					// cheap.. make nicer..
					if(x_dim == 1) {
						return gc_num[0][p_lower]*f_dim;
					} else if(x_dim == 2) {
						return (gc_num[1][p_lower]*extent[0] + gc_num[0][p_lower])*f_dim;
					} else if(x_dim == 3) {
						return  ( gc_num[2][p_lower]*extent[1]*extent[0] + gc_num[1][p_lower]*extent[0] + gc_num[0][p_lower])*f_dim;
					}
					return 0;
				}
				
				int buffer_size() {
					return f_dim*grid_size();
				}
				
				int grid_size() {
					if(x_dim == 1) {
						return extent[0];
					} else if(x_dim == 2) {
						return extent[0]*extent[1];
					} else if(x_dim == 3) {
						return extent[0]*extent[1]*extent[2];
					}
					return 0;			
				}
				
				bool same_underlying_buffer( t_vdf other ) {
					if( field_data == other.field_data) return true;
					return false;
				}
		};

		class t_smooth {
			public:
				//t_smooth();
				//~t_smooth();
				
				// Type of smoothing in each direction
				int type[p_x_dim];

				// Order of smoothing (number of neighboring cells included)
				int order[p_x_dim];
				int max_order;

				// Single pass filter kernel for each direction
				// scoef( 2*maxval(order)+1, p_x_dim ) 
				//p_k_fld, dimension(:,:), pointer :: scoef => NULL()
				std::vector<p_k_fld> scoef[3];

				double digital_A[p_x_dim];			// size of Gibbs phenomenon wiggles in -db
				double digital_fc[p_x_dim];		// cutoff frequency

/*
				int enabled;
				
				// switches for smoothing in each direction
				// integers instead of logical (0 is false, 1 true);
				int ifsmooth[p_x_dim];
				
				// level of smoothing (number of neighboring cells included)
				int smooth_level[p_x_dim];		
				int nsm[p_x_dim];
				
				
				// frequency for no-atenuation
				float__fld omega0;
				
				// smoothing direction
				int dir;
				
				// number of smoothing halo points needs for each direction, given the smoothing stencil paramters.
				// Call the largest value nscoef = X. The actual, final sencil coefficient buffer is of a size (2*X + 1, p_x_dim)
				int nscoef[p_x_dim];
				
				// smoothing weighting factors of each successive smoothing
				// step for the field in x(:) - these are given in the input dec
				// swfj( 1:3 , maxval(nsm), p_x_dim )
				// real(p_k_fld), dimension(:,:,:), pointer :: swfj => NULL()
				//float__fld** swfj;
				// try cheating...
				int swfj__max_middle_idx;
				float__fld swfj[p_x_dim][10][3];
				
				// coefficients of the one step smoothing calculated from swfj
				// scoef( maxval(nscoef)+1, p_x_dim ) 
				//real(p_k_fld), dimension(:,:), pointer :: scoef => NULL() 
				//float__fld** scoef;
				int scoef__max_last_idx;
				float__fld scoef[p_x_dim][25]; */
				
		};

		class t_current {
			public:
				//t_current();
				//~t_current();
				std::vector<t_vdf> pf;

				//t_vdf  pf;
				//t_current_bound_c bnd_con,
				t_smooth smooth;
				//t_current_diag_c diag,
				int coordinates;
				// Position on global grid
				int gix_pos[p_x_dim];
				int bc_type[p_x_dim][2];
				int interpolation;
				
				int jay_sizes__CUDA;
				float__part* jay__CUDA;
				float__part* jay_empty_cheat;
				
				t_vdf  jay_scratch;
		    t_vdf  jay_scratch_new;
		};

		class t_wall {
			public:
				t_wall() {
					idir = -1; idir = -1;
					range[0] = 0; range[1] = 1;
					x_dim = 1;
					f_dim = 1;
					
					for(int i=0;i<p_max_dim;i++) {
						gc_num[i][0] = 0;
						gc_num[i][1] = 0;
						
						nx[i] = 0;
					}
				};
				
				// perp. direction of the wall surface
				int idir;

				// boundary of the wall (lower=1,upper=2)
				int ibnd;

				// range of wall at the boundary i.e. bounds perpendicular
				// to the wall
				int range[2];

			#ifdef __HAS_ALLOC_EXT__
				// allocatable arrays for 1D, 2D and 3D field data
				//real(p_k_fld), dimension(:,:), allocatable :: f1 
				//real(p_k_fld), dimension(:,:,:), allocatable :: f2 
				//real(p_k_fld), dimension(:,:,:,:), allocatable :: f3
				float__fld *f1;
				float__fld *f2;
				float__fld *f3;
				
			#else
				// pointers to 1D, 2D and 3D field data 
				//real(p_k_fld), dimension(:,:), pointer :: f1 => null()
				//real(p_k_fld), dimension(:,:,:), pointer :: f2 => null()
				//real(p_k_fld), dimension(:,:,:,:), pointer :: f3 => null()
				float__fld* *f1;
				float__fld* *f2;
				float__fld* *f3;
			#endif
				
				// spacial dimensions
				int x_dim;

				// field dimensions
				int f_dim;
				
				// grid size
				// integer, dimension(p_max_dim) :: nx = 0
				int nx[ p_max_dim ];
				
				// guard cell information
				int gc_num[p_max_dim][2];
		};


		class t_lindman {
			public:
				//t_lindman();
				//~t_lindman();
				
				int idir;				// direction of lindman bc (1->p_x_dim)
				int orientation;		// upper or lower wall
				int perp[ p_f_dim-1 ];	// perpendicular field directions
				int range[2];			// range of buffers
				
				t_wall e_buffer;		// e and b buffers for storing the
				t_wall b_buffer;		// fields at previous timesteps
				
				float__fld dtdxi[p_max_dim];		// (dt/dx/2)
				float__fld disp_corr;				// dispersion correction
													// (1-dx/(dt*c))/(1+dx/(dt*c))
				int imin[ p_max_dim ];
				int imax[ p_max_dim ];
				int bound[2];
		};	

		class t_vpml {
			public:
			
			//t_vpml();
			//~t_vpml();
				
			//Boundaries are sorted as follows
			// x1_lower, x1_upper, x2_lower, x2_upper, x3_lower, x3_upper

			// Number of vpml boundaries
			int n_vpml;
			
			// Perp. direction of the boundary
			int* dir_vpml;						// needs allocation
			// Boundary location (p_lower, p_upper)
			int* loc_vpml;						// needs allocation

			// Number of cells in vpml boundaries
			// Constant for all boundaries, for now...
			int* bnd_size;						// needs allocation
			
			// i am using int 0==false, 1 == true rather then a 'logical type'
			int if_diffuse[p_x_dim][2];
			
			// dt/dx
			float__fld dtdx[p_x_dim];
			
			// Contact walls for each corner
			//integer, dimension(:,:,:), allocatable :: pos_corner
			int*  pos_corner; 					// needs allocation   , dimension(:,:,:), allocatable ::
		  
		  // Coeficients of E and B in eq (27) of Vay's paper: JCP 165, 511 (2000)
		  // [
		  //  direction  => wall direction (p_x_dim elements)
		  //  coeff. idx => index of coefficients in eq. (27) of Vay's paper (3 elements)
		  //  array pos. => position in the vpml boundary (bnd_size elements) 
		  // ] 
		  /*real(p_k_fld), dimension(:,:,:), allocatable  :: coef_e_low                   
		  real(p_k_fld), dimension(:,:,:), allocatable  :: coef_b_low    
		  real(p_k_fld), dimension(:,:,:), allocatable  :: coef_e_up                   
		  real(p_k_fld), dimension(:,:,:), allocatable  :: coef_b_up*/
		  
		  // 2 walls per vpml boundary
		  /*type (t_wall), pointer, dimension(:) :: wall_array_e => NULL() 
		  type (t_wall), pointer, dimension(:) :: wall_array_b => NULL() */
		  float__fld* coef_e_low;
		  float__fld* coef_b_low;
		  float__fld* coef_e_up;
		  float__fld* coef_b_up;
		  
		  t_wall**  wall_array_e;
		  t_wall**  wall_array_b;
		};
		class  t_emf_bound {

		  // type of boundary conditions for each boundary
		  // p_bc_other_node : boundary to other node
		  // p_bc_periodic   : periodic b.c.
		  // p_bc_conducting : conducting boundary with particle absorption
		  // p_bc_axial      : axial b.c. for cylindrical coordinates
		  // p_bc_lindman    : Lindman open-space boundary - limited implementation
		  // p_bc_vpml       : Vay Perfectly Matched Layer

			public:
				//t_emf_bound();
				//~t_emf_bound();
				
				int type[p_max_dim][2];
				
				t_lindman lower_lindman;
				t_lindman upper_lindman;
				
			#ifdef __HAS_PML__
				t_vpml vpml_all;
			#endif
				
				int vpml_bnd_size;
				int vpml_diffuse[p_x_dim][2];
		};

		class t_emf  {
			public:
				//t_emf();
				//~t_emf();

				t_vdf	b;
				t_vdf	e;
				t_vdf b_part;
				t_vdf e_part;
				bool part_fld_alloc;
				int solver;
				t_vdf	ext_b0;
				t_vdf	ext_e0;
				bool ifextfld;

				// this is state I created and need to fill out better...
				bool is_ext_e;	
				bool is_ext_b;

				t_emf_bound bnd_con;			// emf boundry constions.
				
				// smoothing (before advancing particles_
				int smooth_type;
				int smooth_niter;
				int smooth_nmax;
				t_smooth smooth;
				
				//t_diag_emf
				int coordinates;
				
				// CUDA device pointers.....
				float__fld * e__CUDA;
				float__fld * e_part__CUDA;
				float__fld * b__CUDA;
				float__fld * b_part__CUDA;
				int field_sizes__CUDA;
		    
		    // for debuugin / early development.
		    t_vdf b_scratch;
		    t_vdf e_scratch;
			
		};
	class t_simulation {
		public:

	};

class t_simulation_cuda {
	public:
		CudaInputDeck* CUDA_CONFIG;
		t_node_conf no_co;
		t_grid grid;
		t_emf emf;
		t_particles parts;
		t_space g_space;
		t_current jay;

	t_simulation_cuda(CudaInputDeck* inputDeck) : CUDA_CONFIG(inputDeck) {

	};
};


#include "simple-cpp-profiler.h"

extern t_simulation_cuda CUDA_STATE;
CudaInputDeck* createInputDeckFactory( );

// TODO: HACK move these CUDA references into their own place.
#include <cuda_runtime_api.h>

void FATAL_ERROR( std::string message );
void FATAL_ERROR( std::string message, bool exit_on_failure);
void FATAL_ERROR( std::string message, const char* filename, int line_number  );
void FATAL_ERROR( std::string message, const char* filename, int line_number, bool exit_on_failure);

bool checkcuda(cudaError_t error, bool exit_on_failure);
bool checkcuda(cudaError_t error, std::string message, bool exit_on_failure);
bool checkcuda(cudaError_t error, const char* file, int line, bool exit_on_failure);
bool checkcuda(cudaError_t error, std::string message, const char* file, int line, bool exit_on_failure);

#define CUDA_CHECK(api_call) checkcuda(api_call, __FILE__, __LINE__,true);
#define CUDA_CHECK_DESC(api_call, message) checkcuda(api_call, message, __FILE__, __LINE__, true);

void setup_cuda_device_environment();
int setup_multi_gpu_mpi_mappings( int number_of_mpi_nodes );
void report_GPU_memory_usage( int device_index, size_t * total_mem, size_t* free_mem,size_t* used_mem );


// TODO: move and make nicer.
class my_gpu_environment {
	public:
		//my_gpu_environment();
		//~my_gpu_environment();

		int deviceCount;
		int numberDevicesToUse;
		int p_cache_size_cuda;
		
		int gpu_device_id;
		// this will be an array hold the properties read for the card.
		std::vector<cudaDeviceProp> deviceProps;

		int version;
		int maxgsx;
};

// ----------- HOST allocation functions -----------
#define host_alloc( A, B ) allocate_on_host( A, B, __LINE__, __FILE__)
#define host_alloc( A, B ) allocate_on_host( A, B, __LINE__, __FILE__)

#define host_dealloc( A )  allocate_on_host( A, __LINE__, __FILE__)
template<typename T> T* allocate_on_host( T** ptr, size_t sizeInBytes);
template<typename T> T* allocate_on_host( T** ptr, size_t sizeInBytes,int line_number, const char* filename);
template<typename T> T* deallocate_on_host( T* ptr, size_t sizeInBytes);
template<typename T> T* deallocate_on_host( T* ptr, size_t sizeInBytes,int line_number, const char* filename);

/*
#define host_alloc( A, B ) host_alloc1( A, B, __LINE__, __FILE__)
#define host_dealloc( A ) host_dealloc1( A, __LINE__, __FILE__)

int* host_alloc1( int** buffer, size_t size, int line_number, const char* );
float* host_alloc1( float** buffer, size_t size,  int line_number, const char* );
double* host_alloc1( double** buffer, size_t size,  int line_number, const char* );

void host_dealloc1( int** buffer,  int line_number, const char* );
void host_dealloc1( float** buffer, int line_number, const char* );
void host_dealloc1( double** buffer, int line_number, const char* );
*/

// ----------- DEVICE allocation functions -----------
void setup_heap_manager();
#define gpu_alloc( A, B, C ) 					allocate_on_gpu( A, B, C, #A, __LINE__, __FILE__)
#define gpu_alloc_desc( A, B, C, D ) 	allocate_on_gpu( A, B, C, D, __LINE__, __FILE__)
//#define gpu_alloc( A, B, C, D ) allocate_on_gpu( A, B, C, D, __LINE__, __FILE__)

#define gpu_free( A, B )  deallocate_on_gpu( A, B, __LINE__, __FILE__)
template<typename T> void allocate_on_gpu(	T** ptr, size_t size, int*irc,
																						const char* description,
																						int line_number, const char* filename);
//template<typename T> void allocate_on_gpu(T** ptr, size_t size, int*irc, int line_number, const char* filename);
//template<typename T> void allocate_on_gpu(T** ptr, size_t size, int*irc);
template<typename T> void deallocate_on_gpu(T* ptr, int*irc, int line_number, const char* filename);
//template<typename T> void deallocate_on_gpu(T* ptr, int*irc);

/*
#define gpu_alloc( A, B )   gpu_alloc_( A, B, __LINE__, __FILE__)
#define gpu_dealloc( A )  gpu_dealloc_( A, __LINE__, __FILE__)

int* gpu_alloc_( int** buffer, size_t size, int line_number, const char* );
float* gpu_alloc_( float** buffer, size_t size,  int line_number, const char* );
double* gpu_alloc_( double** buffer, size_t size,  int line_number, const char* );

void gpu_dealloc_( int** buffer,  int line_number, const char* );
void gpu_dealloc_( float** buffer, int line_number, const char* );
void gpu_dealloc_( double** buffer, int line_number, const char* );
*/


extern my_gpu_environment gpu_environment;

#endif