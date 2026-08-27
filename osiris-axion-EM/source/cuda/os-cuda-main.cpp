#include "os-cuda-state.h"

void dumpheap();
using namespace std;
 
void init_species() {
	int * gpu_ptr; int irc;
	float*b1;
	float*b2;
	float*b3;
	setup_heap_manager();
	gpu_alloc( &gpu_ptr, 1000, &irc);
	gpu_alloc( &b1, 1000, &irc);
	gpu_alloc( &b2, 1654, &irc);
	gpu_alloc( &b3, 6654, &irc);
	dumpheap();
	gpu_free(b2,&irc);
	dumpheap();
/*
	//hm.add_entry( species.pbuf, "pbuf", 34*100047+3, "particle" );
	//hm.add_entry( species.pbuf, "pbuf", 34*100047+3, "particle" );
	//hm.add_entry( species.mpi_x, "pbuf", 34*100047+3, "particle" );
	//hm.add_entry( species.pbuf, "pbuf", 34*100047+3, "particle" );

	hm.free_group( "particle");
	hm.free_groups( vector<string> );
	hm.free( vector<string> );
	hm.alloc( vector<string> );
	hm.size_requested();
	hm.size_requested_group();
*/

	//gpu_free( gpu_ptr, &irc);
}

void init_logs( int rank ) {
	if( CUDA_STATE.CUDA_CONFIG->is_root() ) {
		LOG_OBJECT().init(rank);
		TimerOracle::summon()->write_output = true;
	} else {
		TimerOracle::summon()->write_output = false;
	}

 	// enable timing...
 	TimerOracle::summon()->global_enable_timing();

 	// disable timing
 	//TimerOracle::summon()->global_disable_timing();

 	// test if enabled
 	//TimerOracle::summon()->is_enabled();
 
 	// get specific timer results
	//TimerOracle::summon()->get_result("cool timer name");

	/* to dump the timing logs do this:
			std::stringstream ss;
			ss << " {} ";
			TimerOracle::summon()->dump_nvpair( 0.0, initial_time_step,1000,1000.0, 1,ss.str(),  ss.str() );
	*/
}

void display_banner_configuration_options() {

}

/*
---------------------------------------------------------------
				Code that handles simulation RESTART events

							init_from_restart
							prepare_for_restart		-> calls write_restart_file
---------------------------------------------------------------
*/
void init_from_restart( ) {

	if(!CUDA_STATE.CUDA_CONFIG->is_restart_run()) return;

	if(CUDA_STATE.CUDA_CONFIG->is_root()) {
		std::cout << "\t This is a RESTART at time step " << CUDA_STATE.CUDA_CONFIG->timestep_at_initalization << " attempting to restart the CUDA state." << std::endl;
	}
	stringstream fname;
	int init_step = CUDA_STATE.CUDA_CONFIG->timestep_at_initalization;

	fname << "RE/CUDA_restart__" << CUDA_STATE.CUDA_CONFIG->node_id << "__" << init_step  << ".txt";
	//cout << " filename " << fname.str() << endl;
	//out.open(fname.str().c_str());
	
	ifstream restart_input_stream(fname.str().c_str()); 
	if( restart_input_stream ) {
		restart_input_stream.close();
		CudaInputDeck* restart_input = createInputDeckFactory();
		restart_input->read_configuration_from_file(  fname.str().c_str() );
		
		
		if(CUDA_STATE.CUDA_CONFIG->is_root()) {
			std::cout << "\t\t CUDA Osiris restart files found." << endl; 
		}

		for(int i=0;i<CUDA_STATE.parts.num_species;i++) {
			stringstream varname0; varname0 << "nppmx_" << i;
			int nppmx =  restart_input->get_value_int(varname0.str());
			
			stringstream varname1; varname1 << "Particle Array Size Factor " << i;
			float nppmx0 = ((float) restart_input->get_value_float(varname1.str()));

			stringstream varname2; varname2 << "Particle Sort Hole Array Size Factor " << i;
			float ntmaxp = ((float) restart_input->get_value_float(varname2.str()));

			stringstream varname3; varname3 << "Particle Sort Scratch Array Size Factor " << i;
			float npbmx  = ((float) restart_input->get_value_float(varname3.str()));

			stringstream varname4; varname4 << "Particle MPI Array Size Factor " << i;
			float nbmax = ((float) restart_input->get_value_float(varname4.str()));

			//viktor_contexts__species[i].config__use_this_nppmx = nppmx;
			//viktor_contexts__species[i].nppmx_initial = nppmx;
			
			CUDA_STATE.CUDA_CONFIG->set_value("Particle Array Size Factor", nppmx0);
			CUDA_STATE.CUDA_CONFIG->set_value("Particle MPI Array Size Factor", nbmax);
			CUDA_STATE.CUDA_CONFIG->set_value("Particle Sort Hole Array Size Factor",ntmaxp);
			CUDA_STATE.CUDA_CONFIG->set_value("Particle Sort Scratch Array Size Factor", npbmx);
		}

		delete restart_input;

		// also try to restore random number generator state....
		//ifstream myfile ("example.bin", ios::in | ios::binary);
	} else {
		if(CUDA_STATE.CUDA_CONFIG->is_root()) {
			cout << "\t\t CUDA Osiris restart files NOT found... CUDA configuration will be recomputed." << endl; 
		}
	}
}

void write_restart_file( int * n) {
	stringstream fname;

	fname << "RE/CUDA_restart__" << CUDA_STATE.CUDA_CONFIG->node_id << "__" << (*n) << ".txt";
	ofstream out;
	out.open(fname.str().c_str());
	out << "! BEGIN CUDA CONFIGURATION" << endl;
	out << "! { " << endl;
	//out << "! 'Particle Array Size Factor': " << ((float) CUDA_CONFIG.get_value("Particle Array Size Factor")) << "," << endl;
	//out << "! 'Particle Sort Hole Array Size Factor': " << ((float) CUDA_CONFIG.get_value("Particle Sort Hole Array Size Factor")) << "," << endl;
	//out << "! 'Particle Sort Scratch Array Size Factor': " << ((float) CUDA_CONFIG.get_value("Particle Sort Scratch Array Size Factor")) << "," << endl;
	//out << "! 'Particle MPI Array Size Factor': " << ((float) CUDA_CONFIG.get_value("Particle MPI Array Size Factor")) << "," << endl;
	out << "! 'Cuda Block Size': 12," << endl;
	out << "! 'Cuda Tile Size': 12," << endl;
	
	for(int i=0;i<CUDA_STATE.parts.num_species;i++) {
		out << "! 'Particle Array Size Factor "<< i <<"': " << ((float) CUDA_STATE.CUDA_CONFIG->get_value_float("Particle Array Size Factor")) << "," << endl;
		out << "! 'Particle Sort Hole Array Size Factor "<< i <<"': " << ((float) CUDA_STATE.CUDA_CONFIG->get_value_float("Particle Sort Hole Array Size Factor")) << "," << endl;
		out << "! 'Particle Sort Scratch Array Size Factor "<< i <<"': " << ((float) CUDA_STATE.CUDA_CONFIG->get_value_float("Particle Sort Scratch Array Size Factor")) << "," << endl;
		out << "! 'Particle MPI Array Size Factor "<< i <<"': " << ((float) CUDA_STATE.CUDA_CONFIG->get_value_float("Particle MPI Array Size Factor")) << "," << endl;
	}
	

	for(int i=0;i<CUDA_STATE.parts.num_species;i++) {
		stringstream varname0; varname0 << "nppmx_" << i;
		//int nppmx_initial = viktor_contexts__species[i].nppmx_initial;
		// TODO: HACK.. TAKE OUT
		int nppmx_initial = 100;
		if( i > 0) out << "," << endl;
		out << "! '" << varname0.str() << "': " << nppmx_initial;
	}
	out << endl;
	out << "! } " << endl;
	out << "! END CUDA CONFIGURATION" << endl;
	
	out.close();
	
}

void prepare_for_restart( int *n ) {
	cout << endl << "Preparing for Restart: " << endl;
	cout << "\t Copying fields to Osiris (if necessary): " << endl;
	//__COPY_FIELDS_FROM_GPU_TO_OSIRIS_FOR_DAIGNOSTIC__(n);
	cout << "\t Copying Current to Osiris (if necessary): " << endl;
	//__COPY_JAY_FROM_GPU_TO_OSIRIS_FOR_DAIGNOSTIC__(n);
	// TODO: keep track of last time step we sent particlesd and dont repeat if we already did it this timestep.
	cout << "\t Copying Particles to Osiris: " << endl;
	
	for(int i=0;i<CUDA_STATE.parts.num_species;i++) {
		int fortran_species_id = i + 1;
		//__SYNCH_PARTICLE_BUFFERS_CUDA__(&fortran_species_id);
		if( CUDA_STATE.CUDA_CONFIG->is_root() ) {
			cout << "\t\t Particle Buffer copied from GPU to Host for restart (species " << fortran_species_id <<  "). " << endl;
		}
	}
	write_restart_file( n );
}

extern "C"  void osiris_cuda__quickening( ) {

	//std::cout << "-------------- We Are in Osiris Cuda's Quickening ------------" << std::endl;
	//std::cout << "The inital time step is: " << (CUDA_STATE.CUDA_CONFIG->timestep_at_initalization) << std::endl;
	//std::cout << "On MPI hostname: " << CUDA_STATE.CUDA_CONFIG->mpi_hostname << std::endl;

	// initalize the profiling/logging for CUDA.
	//	Its very basic and we may want to move over to 3rd party libraries one day.
	init_logs(CUDA_STATE.no_co.my_MPI_rank);

	setup_cuda_device_environment();
	init_species();
}


 
/*
----------------------------------------------------------
TODO: Move this stuff out of this file to real resting place
---------------------------------------------------------- 
*/
template<typename T> T* allocate_on_host( T** ptr, size_t sizeInBytes) {
	(*ptr) = new T[ sizeInBytes ];
	return *ptr;
}
template<typename T> T* allocate_on_host( T** ptr, size_t sizeInBytes,int line_number, const char* filename) {
	(*ptr) = new T[ sizeInBytes ];
	return *ptr;
} 
template<typename T> void deallocate_on_host( T* ptr) {
	delete[] ptr;
}
template<typename T> void deallocate_on_host( T* ptr, int line_number, const char* filename) {
	delete[] ptr;
}