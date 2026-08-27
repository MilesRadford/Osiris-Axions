#include "os-cuda-state.h"
#include "os-cuda-device-runtime-interaction.h"
#include <cuda_runtime_api.h>

// part of c++ Standard Template Library that creates a 'dictionary'/hashtable
//	(i.e. something similiar to Python dictionaries.)
#include <unordered_map>
#include <list>
#include <map>
#include <vector>
#include <string>

class mem_entry {
	public:
		size_t	alloc_size;
		size_t	alloc_size_bytes;
		size_t	location;
		std::string 	name;
		int 		alloc_number;
};

class simple_heap_watcher {
	public:
		std::map<size_t, mem_entry> allocations;
		int number_allocations;

		simple_heap_watcher() : number_allocations(0) {};

		void record( mem_entry &entry ) {
			allocations[entry.location] = entry;
			allocations[entry.location].alloc_number = number_allocations++;
		};

		void record( size_t location, size_t size, size_t size_bytes, const char* desc = 0 ) {
			mem_entry entry;
			if(desc > 0)
				entry.name.assign(desc);

			entry.location = location;
			entry.alloc_size = size;
			entry.alloc_size_bytes = size_bytes;
			entry.alloc_number = number_allocations++;

			allocations[entry.location] = entry;
		};

		mem_entry* get_entry(size_t location) {
			if(allocations.count(location) > 0)
				return &allocations[location];

			return NULL;
		}

		bool contains(size_t location) {
				return (allocations.count(location) > 0);
		}

		void remove(size_t location) {
			size_t g = allocations.erase(location);
			//std::cout << "ERASED " << g << std::endl;
		}

		void dump_heap() {
			size_t previous_location = 0;
			size_t previous_size_bytes = 0;
			
			size_t delta = 0;
			size_t slack_space = 0;
			std::cout << std::endl;
			std::cout << "&&&&&&&&&&&&&&&&&&&&&&& " << "Dump of Heap Allocations" << " &&&&&&&&&&&&&&&&&&&&&&&" << std::endl;
			for (std::map<size_t, mem_entry>::iterator it=allocations.begin(); it!=allocations.end(); ++it) {
				mem_entry& entry = it->second;
				std::cout << "&&&&&&&&&&&&&&&&&&&&&&&" << std::endl;
				//std::cout << it->first << " => " << it->second.name << '\n';
				if( previous_location > 0) {
					delta = entry.location - previous_location;
					slack_space = delta - previous_size_bytes;
				}

				std::cout << entry.location << "\t" << entry.alloc_size 
									<< " (" << entry.alloc_size_bytes << " bytes) "
									<< "\t" << entry.alloc_number
									<< "\t" << entry.name;
				if( previous_location > 0 ) {
					std::cout << "\t" << delta << " (" << slack_space << ") \t";
				}
				std::cout << std::endl;

				previous_location = entry.location;
				previous_size_bytes = entry.alloc_size_bytes;

			}
		}

};

//std::map<size_t, mem_entry> allocations;
	simple_heap_watcher& get_gpu_heap_watcher() {
		static simple_heap_watcher temp;
		return temp;
	}


void setup_cuda_device_environment() {
	//cudaError_t error;

	int deviceCount;
	//int dev;
	//int p_cache_size_cuda;
	//struct cudaDeviceProp *devicePropArray;
 
	// count the number of valid CUDA devices on this processes's computer.
	CUDA_CHECK( cudaGetDeviceCount(&deviceCount) )

	if( CUDA_STATE.CUDA_CONFIG->is_root() ) {
		printf("       Found %d devices\n", deviceCount);
	}

	gpu_environment.deviceProps.resize( deviceCount );
	//cudaDeviceProp *deviceProps = new cudaDeviceProp[deviceCount]; 

	for (int dev = 0; dev < deviceCount; ++dev) {
		struct cudaDeviceProp *deviceProp = &( gpu_environment.deviceProps[dev] );
		cudaGetDeviceProperties(deviceProp, dev);

		if( CUDA_STATE.CUDA_CONFIG->is_root() ) {
			printf("          ------------------Device %d------------------------\n", dev);
			
			printf("          number mulitprocessors: %d\n", deviceProp->multiProcessorCount);
			printf("          total global memory: %lld\n", deviceProp->totalGlobalMem);
			printf("          shared memory per block: %lld\n", deviceProp->sharedMemPerBlock);
			printf("          registers per block: %d\n", deviceProp->regsPerBlock);
			printf("          warp size: %d\n", deviceProp->warpSize);
			printf("          threads per block: %d\n", deviceProp->maxThreadsPerBlock);
			printf("          clock rate: %d\n", deviceProp->clockRate);
			printf("          amount Constant memory : %lld\n", deviceProp->totalConstMem);
			printf("          number asynch actions allowed: %d\n", deviceProp->asyncEngineCount);
			printf("          L2 Cache: %d\n", deviceProp->l2CacheSize);
			printf("          Involved in unified addressing? (1=yes): %d\n", deviceProp->unifiedAddressing);
		}
	}
	if( CUDA_STATE.CUDA_CONFIG->is_root() ) {
		printf("          --------------------------------------------------\n\n"); 
	}

	// TODO: HACK HARDCODED. MAYBE JUST TAKE OUT SINCE UNUSED.
	//		We have set/hardcoded the number of GPUs to be used by each (MPI)process to one.
	//		at this point, usinng multiple GPUs per MPi Node is not supported.
	gpu_environment.numberDevicesToUse = 1;
	struct cudaDeviceProp& prop = gpu_environment.deviceProps[0];

	// thesee are not used as of now...
	// TODO: MAYBE TAKE THESE OUT
	gpu_environment.version = 10*prop.major + prop.minor;
	gpu_environment.maxgsx = prop.maxGridSize[0];

	// decide the p_cache_size_cuda
	// gpu_environment.deviceProps[0].maxThreadsPerBlock = 1024 on Dawson 2
	// for npw, chunk things into the size 1024*448 (number of multiprocessors.
	//p_cache_size_cuda = gpu_environment.deviceProps[0].maxThreadsPerBlock * gpu_environment.deviceProps[0].multiProcessorCount;
	//p_cache_size_cuda = 1024;

	// thesee are not used as of now...
	// TODO: MAYBE TAKE THESE OUT
	gpu_environment.p_cache_size_cuda = 1024;
	
	// see if the user specified that a certainn GPU be used.
	//		TODO: allow more flexability to select GPUs by name or CUDA version number
	int force_device = CUDA_STATE.CUDA_CONFIG->get_forced_cuda_device();
  if(force_device > -1) {
		std::cout <<  " CUDA Device Forced (via \"Use Cuda device\" input parameter) to " << force_device << std::endl;
  } else {
  	// the user did not force a certain GPU to be used, so attempt to assign the GPUs
  	//		in an orderly fashion. In particular, each MPI process should get 1 GPU to use.
		int gpu_idx_this_node_should_use = setup_multi_gpu_mpi_mappings( CUDA_STATE.no_co.no_num );
  	
  	// check that all nodes where assigned to GPUs successfully.
  	//		forst see if this node hd an error
  	int any_error = 0; int global_error_occured = 0;
  	if( gpu_idx_this_node_should_use < 0 ) any_error = 1;

  	// Now gossip and tell all nodes if any node had an error. Each node has an interger
  	//		of value 0 (no errors) or 1 (if error). Then we will use MPI get the sum of this integer across all nodes
  	//		and broadcast the results so that every node will get the sum (and thus all nodes get smae value)
  	//	 	If no error occured, then the sum will be zero. Otherwise the sum will be >0 and one or more nodes did have an error.
  	MPI_Allreduce(	&any_error,
  				&global_error_occured, 1, 
  				MPI_INT, MPI_SUM, MPI_COMM_WORLD );

  	if( global_error_occured > 0 && CUDA_STATE.CUDA_CONFIG->is_root() ) {
  		std::cout << " An error occured when assigning MPI nodes to GPUS. Aborting." << std::endl;
  		 MPI_Abort( MPI_COMM_WORLD, -1);
  	}
  	force_device = global_error_occured;
  }

  gpu_environment.gpu_device_id = force_device;
  CUDA_CHECK_DESC( cudaSetDevice(force_device), "cudaSetDevice: associating GPU with MPI node." )
}


//  ns.no_co.no_num
// TODO: do proper error checking
int setup_multi_gpu_mpi_mappings( int number_of_mpi_nodes ) {
	
	int this_node_gpu_assignment = -1;
	
	char * temp_buffer_root = 0;
	int * gpu_idx_assignments = 0;
	
	// get this node's hostname from MPI.
	char temp_buffer[ MPI_MAX_PROCESSOR_NAME+1];
	int temp_buffer_length = 0;

	// MPI_MAX_PROCESSOR_NAME:	constant defined by MPI to provide the maximum 
	//													possible character buffer size needed
	MPI_Get_processor_name( temp_buffer, & temp_buffer_length);

	if( CUDA_STATE.CUDA_CONFIG->is_root() ) {
		std::cout << std::endl << "\t ---------- GPU Assignments for all the MPI nodes -----------" << std::endl << std::endl;
		
		// create an array to hold all the host names from all processors. And zeros it for parinoia's sake.
		temp_buffer_root = new char[ number_of_mpi_nodes*MPI_MAX_PROCESSOR_NAME];

		for(int i=0;i<(number_of_mpi_nodes*MPI_MAX_PROCESSOR_NAME);i++) {
			temp_buffer_root[i] = 0;
		}
		gpu_idx_assignments = new int[ number_of_mpi_nodes ];
		for(int i=0;i<(number_of_mpi_nodes);i++) {
			gpu_idx_assignments[i] = -1;
		}		
	}
	
	std::string temp_string1;
	temp_string1.assign(temp_buffer );
	MPI_Gather( temp_buffer, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, temp_buffer_root, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, 0,  MPI_COMM_WORLD);
	
	if( CUDA_STATE.CUDA_CONFIG->is_root() ) {
		// terminate a C-strings just in case (worst case seniero)
		for(int i=1;i<=(number_of_mpi_nodes);i++) {
			temp_buffer_root[i*MPI_MAX_PROCESSOR_NAME - 1] = 0;
		}
		
		std::unordered_map< std::string, int >  mpi_host__gpu_idx__mapping;
		
		for(int i=0;i<(number_of_mpi_nodes);i++) {
			std::string temp_string;
			temp_string.assign( &(temp_buffer_root[i*MPI_MAX_PROCESSOR_NAME]) );
			
			std::cout << "\t\t MPI NODE " << i << " sent root the hostname '" <<  temp_string << "'" << std::endl;
			// if we have not seen this hostname yet, put in a default value of zero.
			if (!mpi_host__gpu_idx__mapping.count(temp_string)) {
				//cout << "\t\t hostname '" <<  temp_string << "' Not seen yet (not in hash table).. Setting hash value to 0" << endl;
				mpi_host__gpu_idx__mapping[temp_string] = 0;
			}
			
			// now pull out how many times we have seen the current hostname
			int times_host_name_seen = mpi_host__gpu_idx__mapping[ temp_string ];
			std::cout << "\t\t\t hostname '" <<  temp_string << "' has been seen " << times_host_name_seen << " times.. so mpi node " << i <<" get's the GPU assignment " << times_host_name_seen << std::endl;
			gpu_idx_assignments[i] = times_host_name_seen++;
			
			// save the incremented values back to the hash table.
			// cout << "\t\t hostname '" <<  temp_string << "' is putting the value " <<  times_host_name_seen << " into the hash for hostname " << temp_string << endl;			
			mpi_host__gpu_idx__mapping[temp_string] = times_host_name_seen;
		}
		std::cout << std::endl << "\t ------------------------------------------------" << std::endl << std::endl;	
	}
	
	MPI_Scatter( gpu_idx_assignments, 1, MPI_INT, &this_node_gpu_assignment, 1,  MPI_INT, 0, MPI_COMM_WORLD);		
	
	if( CUDA_STATE.CUDA_CONFIG->is_root() ) {
		delete [] temp_buffer_root;
		delete [] gpu_idx_assignments;
	}
	
	return this_node_gpu_assignment;
}

void report_GPU_memory_usage( int device_index, size_t * total_mem, size_t* free_mem,size_t* used_mem ) {
	*total_mem 	= 0;
	*free_mem 	= 0;
	*used_mem 	= 0;

	CUDA_CHECK_DESC( cudaMemGetInfo( free_mem, total_mem), "Calling 'cudaMemGetInfo' to get the amount of memory used on the GPU."  )
	*used_mem = (*total_mem - *free_mem);

}

/*
	Simple mem management/tracking of the GPU gloabl memory 'heap'
*/
	/*
template<typename T> void allocate_on_gpu(T** ptr, size_t size, int*irc, int line_number, const char* filename) {
	void *gptr;
	int crc = cudaMalloc( &gptr, size*sizeof(T) );
	if( crc ) *irc = 1;
	*ptr = (T*) gptr;
}
template<typename T> void allocate_on_gpu(T** ptr, size_t size, int*irc) {
	void *gptr;
	int crc = cudaMalloc( &gptr, size*sizeof(T) );
	if( crc ) *irc = 1;
	*ptr = (T*) gptr;
}
*/

template<typename T> void allocate_on_gpu(	T** ptr, size_t size, int*irc,
																						const char* description,
																						int line_number, const char* filename) {
	void *gptr;
	int crc = cudaMalloc( &gptr, size*sizeof(T) );
	if( crc ) *irc = 1;
	*ptr = (T*) gptr;

	get_gpu_heap_watcher().record( (size_t) gptr, size, size*sizeof(T), description);
}


template<typename T> void deallocate_on_gpu(T* ptr, int*irc, int line_number, const char* filename) {
	int crc = cudaFree( ptr );
	if( crc ) *irc = 1;
	get_gpu_heap_watcher().remove( (size_t) ptr );
}
template<typename T> void deallocate_on_gpu(T* ptr, int*irc) {
	int crc = cudaFree( ptr );
	if( crc ) *irc = 1;
	get_gpu_heap_watcher().remove(  (size_t) ptr );
}

template void allocate_on_gpu(int**,		size_t, int*, const char*, int, const char*);
template void allocate_on_gpu(float**,	size_t, int*, const char*, int, const char*);
template void allocate_on_gpu(double**,	size_t, int*, const char*, int, const char*);
/*
template void allocate_on_gpu(int**,size_t,int*,int, const char*);
template void allocate_on_gpu(float**,size_t,int*,int, const char*);
template void allocate_on_gpu(double**,size_t,int*,int, const char*);

template void allocate_on_gpu(int**,size_t,int*);
template void allocate_on_gpu(float**,size_t,int*);
template void allocate_on_gpu(double**,size_t,int*);
*/

template void deallocate_on_gpu(int*,int*,int, const char*);
template void deallocate_on_gpu(float*,int*,int, const char*);
template void deallocate_on_gpu(double*,int*,int, const char*);

template void deallocate_on_gpu(int*,int*);
template void deallocate_on_gpu(float*,int*);
template void deallocate_on_gpu(double*,int*);


void setup_heap_manager( ) {
}

void dumpheap() {
	get_gpu_heap_watcher().dump_heap();
}
