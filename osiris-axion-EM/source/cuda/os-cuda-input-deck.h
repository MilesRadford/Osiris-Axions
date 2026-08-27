#ifndef OS_CUDA_INOUT_DECK_H
#define OS_CUDA_INOUT_DECK_H
	
	#include <math.h>
	#include <vector>
	//#include <algorithm>  
	#include <string>
	// TODO: take followin 2 lines out
	
	#include <iostream>
	#include <iomanip>
	#include <sstream> 
	#include <fstream>
	
	#include "mpi.h"

	#include "os-cuda-basic-typedefs.h"

	class CudaInputDeck {
		public:
			std::string configuration_string;
			bool configured_from_file;
			std::string configuration_filename;
			int node_id;
			int timestep_at_initalization;
			
			int default_cuda_block_size;
			int tile_size;
			int tile_size__x;
			int tile_size__y;
			int tile_size__z;
			int force_cuda_device;

			std::string mpi_hostname;

			virtual int get_forced_cuda_device() = 0;
			virtual int get_opengl_display_mode() = 0;

			virtual void set_node_id( int id) = 0;
			virtual void read_configuration_from_file( std::string filename) = 0;
			virtual void read_configuration_from_string( std::string configuration) = 0;

			virtual void dump_config( ) = 0;
			virtual void set_value( std::string key, float__part val) = 0;
			virtual float get_value_float( std::string key) = 0;
			virtual double get_value_double( std::string key) = 0;
			virtual int get_value_int( std::string key) = 0;
			virtual std::string get_value_string( std::string key) = 0;

			virtual bool is_root() = 0;
			virtual bool is_restart_run() {
				// During initialization of a new simulation, the integer time (i.e stepstep) is n=-1
				//	timestep is advanced first thing in at the top of the main simulation loop, so n=0 during 
				//	first simulation timestep. If n>-1 then its because we are restarting the simulation from
				//	checkpoint files.
				if(timestep_at_initalization == -1) {
					return false;
				}
				return true;
			};
	};

#endif
