#include "os-cuda-input-deck.h"

	//#include "fortran.h"
	//#include "hdf5.h"

	#define OC_NEW_STYLE_INCLUDES

	#if defined(_MSC_VER) || defined(_WIN32) || defined(_WIN64)
		#if defined(_M_AMD64) || defined(_WIN64)
			#define OC_BYTES_IN_POINTER 8 
		#else
			#define OC_BYTES_IN_POINTER 4
			fff
		#endif
	#else
		#define OC_BYTES_IN_POINTER 8
	#endif

	//#else
	//	#define OC_BYTES_IN_POINTER 4
	//	#error Only Debug builds are supported
	//#endif

	////#include "chooseser.h"
	#include "ocval.h"
	#include "ocserialize.h"

	using namespace std;

	class CudaOsirisGlobalConfiguration : CudaInputDeck {
		
		private:
			static Val null_val;
			
			// bootstraping data...
			//std::string configuration_string;
			//bool configured_from_file;
			//std::string configuration_filename;
			//int node_id;
			Tab cuda_config_params;
			Tab* cuda_config_params_ptr;
			//string mpi_hostname;
			
			
			//Tab&  cuda_config_params_ref;
			//int default_cuda_block_size;
			//int tile_size;
			//int tile_size__x;
			//int tile_size__y;
			//int tile_size__z;
			//int force_cuda_device;
		public:		
			// NO COPYING
			CudaOsirisGlobalConfiguration(CudaOsirisGlobalConfiguration const&);
			// NO SETTING VIA '='
			CudaOsirisGlobalConfiguration& operator=(CudaOsirisGlobalConfiguration const&); 
			
			// CudaOsirisGlobalConfiguration() : node_id(0), default_cuda_block_size(12),configured_from_file( false ),force_cuda_device(-1)  {
			CudaOsirisGlobalConfiguration() {
				node_id = 0;
				default_cuda_block_size = 12;
				configured_from_file = false;
				force_cuda_device = -1;
			}
			
			int get_forced_cuda_device();
			int get_opengl_display_mode();
			string& get_hostname();
			
			void set_node_id( int id) { node_id = id; }
			void read_configuration_from_file( std::string filename);
			void read_configuration_from_string( std::string configuration);
			////Tab& operator () ();
			////Val& operator () (string key);
			void dump_config( );
			Tab& get_value();
			Val& get_value(string key);
			void set_value( string key, float__part val);
			
			float get_value_float( std::string key);
			double get_value_double( std::string key);
			int get_value_int( std::string key);
			std::string get_value_string( std::string key);

			bool is_root() {
				if( node_id == 0) return true;
				return false;
			}
	};

Val CudaOsirisGlobalConfiguration::null_val(0);

void CudaOsirisGlobalConfiguration::read_configuration_from_file( std::string filename) {
	
	configuration_filename = filename;
	
	// try to read in cuda paramters....
	std::stringstream input_config;
	string file_line;
	ifstream myReadFile;
	myReadFile.open(configuration_filename.c_str());
	Tab cuda_config_params;
	
	bool in_cuda_config_block = false;
	bool cuda_config_properly_closed = false;
	
	while(! (myReadFile.eof() ) ) {
		getline(myReadFile,file_line);
		if( in_cuda_config_block ) {
			size_t index = file_line.find("END CUDA CONFIGURATION");
			if(index != string::npos) {
				in_cuda_config_block = false;
				cuda_config_properly_closed = true;
				break;
			}
			
			// take out any leading '!' (fortran namelist comment chararacter)
			size_t first_bang = file_line.find("!");
			if(first_bang !=  string::npos) {
				file_line = file_line.substr( first_bang+1, string::npos);
			}
			
			
			input_config << file_line;
			
		} else {
			size_t index = file_line.find("BEGIN CUDA CONFIGURATION");
			if(index !=  string::npos) {
				in_cuda_config_block = true;
			}
		}
	}
	myReadFile.close();
	
	if( cuda_config_properly_closed ) {
		
		this->cuda_config_params = input_config.str();
		
		default_cuda_block_size = 	this->cuda_config_params("Cuda Block Size");
		tile_size = 				this->cuda_config_params("Cuda Tile Size");
		tile_size__x =				this->cuda_config_params("Cuda Tile Size");
		tile_size__y =				this->cuda_config_params("Cuda Tile Size");
		if( (this->cuda_config_params).contains( "Cuda Tile Y Size" )) tile_size__y = this->cuda_config_params(  "Cuda Tile Y Size" );
		else {
			if( p_x_dim > 1 ) {
				tile_size__y = this->cuda_config_params("Cuda Tile Size");
			} else {
				tile_size__y = 1;
			}
		}
		if( (this->cuda_config_params).contains( "Cuda Tile Z Size" )) tile_size__z = this->cuda_config_params(  "Cuda Tile Z Size" );
		else {
			if( p_x_dim > 2 ) {
				tile_size__z = this->cuda_config_params("Cuda Tile Size");
			} else {
				tile_size__z = 1;
			}
		}
		
		
		
		// save the configuration string...
		configuration_string = input_config.str();
		configured_from_file = true;
		configuration_filename = filename;
		
	} else {
		if(is_root()) { 
			cout << "      Warning: No Cuda Section was found in the Osiris input deck." << endl;
			cout << "      It is recommended that you place such a section in. A quick example is: " << endl;
			cout << "         !BEGIN CUDA CONFIGURATION" << endl;
			cout << "         !{      'Cuda Block Size': 64," << endl << "!       'Cuda Tile Size': 12 }" << endl;
			cout << "         !END CUDA CONFIGURATION" << endl << endl;
			cout << "      Since no configuration was found, baseline defaults will be used." << endl << endl;
		}
	}
	 
	// TODO: move this
	char * temp_buffer = new char[ MPI_MAX_PROCESSOR_NAME];
	int temp_buffer_len = 0;
	
	MPI_Get_processor_name( temp_buffer, & temp_buffer_len);
	mpi_hostname.assign( temp_buffer );
	delete [] temp_buffer;
}

void CudaOsirisGlobalConfiguration::read_configuration_from_string( string input_config) {
	configuration_string = input_config;
	this->cuda_config_params = input_config;
	
	default_cuda_block_size = 	this->cuda_config_params("Cuda Block Size");
	tile_size = 				this->cuda_config_params("Cuda Tile Size");
	tile_size__x = 				this->cuda_config_params("Cuda Tile Size");
	//tile_size__y = 				this->cuda_config_params("Cuda Tile Size");
	//tile_size__z =				this->cuda_config_params("Cuda Tile Z Size");
	if( (this->cuda_config_params).contains( "Cuda Tile Y Size" )) tile_size__y = this->cuda_config_params(  "Cuda Tile Y Size" );
	else {
		if( p_x_dim > 1 ) {
			tile_size__y = this->cuda_config_params("Cuda Tile Size");
		} else {
			tile_size__y = 1;
		}
	}
	if( (this->cuda_config_params).contains( "Cuda Tile Z Size" )) tile_size__z = this->cuda_config_params(  "Cuda Tile Z Size" );
	else {
		if( p_x_dim > 2 ) {
			tile_size__z = this->cuda_config_params("Cuda Tile Size");
		} else {
			tile_size__z = 1;
		}
	}
	
	// TODO: move this
	char * temp_buffer = new char[ MPI_MAX_PROCESSOR_NAME];
	int temp_buffer_len = 0;
	
	MPI_Get_processor_name( temp_buffer, & temp_buffer_len);
	mpi_hostname.assign( temp_buffer );
	delete[] temp_buffer;
}

void CudaOsirisGlobalConfiguration::dump_config( ) {
	cout << " ADDRESS of 'cuda_config_params' class is " << &(this->cuda_config_params) << " (as dump_config()))" << endl;
	cout << this->cuda_config_params << endl;
}

Tab& CudaOsirisGlobalConfiguration::get_value() {
	return this->cuda_config_params;
}

void CudaOsirisGlobalConfiguration::set_value( string key, float__part val) {
	this->cuda_config_params[key] = val;
}

Val& CudaOsirisGlobalConfiguration::get_value(string key) {
	//if( (this->cuda_config_params).contains( key )) {
	//	cout << "\t\t Value found for key '" << key << "' is " << this->cuda_config_params(key) << endl;  //Looking to return data for the key " << key << endl;
	//}
	
	if( (this->cuda_config_params).contains( key )) return this->cuda_config_params( key);
	
	//cout << "Error in CudaOsirisGlobalConfiguration. An unknown variable named '" << key << "' was asked for." << endl;
	return CudaOsirisGlobalConfiguration::null_val;
}

string& CudaOsirisGlobalConfiguration::get_hostname() {
	return mpi_hostname;
}

int CudaOsirisGlobalConfiguration::get_forced_cuda_device() {
	if( (this->cuda_config_params).contains( "Use Cuda device" )) return this->cuda_config_params( "Use Cuda device" );
	return -1;
}

int CudaOsirisGlobalConfiguration::get_opengl_display_mode() {
	if( (this->cuda_config_params).contains( "Opengl Display Mode" )) return this->cuda_config_params( "Opengl Display Mode" );
	return -1;
}

int CudaOsirisGlobalConfiguration::get_value_int( std::string key) {
	return this->get_value(key);
}
float CudaOsirisGlobalConfiguration::get_value_float( std::string key) {
	return this->get_value(key);
}
double CudaOsirisGlobalConfiguration::get_value_double( std::string key) {
	return this->get_value(key);
}
std::string CudaOsirisGlobalConfiguration::get_value_string( std::string key) {
	return this->get_value(key);
}


 
CudaInputDeck* createInputDeckFactory( ) { 
	return (CudaInputDeck*)  new CudaOsirisGlobalConfiguration();
}