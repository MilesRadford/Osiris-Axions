#ifndef OS_CUDA_BASE_INCLUDE_H
	#define OS_CUDA_BASE_INCLUDE_H
	
	#include <stddef.h>

	// This is the C99 header file with the fixed width data types.. it should be availiable on all systems so use this for now.
	#include <stdint.h>

	// This is the C++ header file with the fixed width data types.. it does all  the sathat stdint.h does and it also adds all
	//	the type defintions to the 'std'namespace. It's part of that standard as of C++11
	//#include <cstdint>
	
		
	// signed versus unsigned:
	//	unsigned is not same as positive... its error prone
	//	signed is slower for % modulus stuff
	
	// type defintions... for ease, some are actually compilier macros.. which amounts to basiclly same thing
	
	// Its a maze of standards-draft-speak and snarky comments on stack-overflow about the 'right', 'portable' way
	//	to gettypes of known, specific size..
	//	So I will abstract the C99 <stdint.h> types.. which for most Posix/windows systems will be just fine.. to be
	//	safe, I will abstract the type behind typedefs so they can easily be changed if a portability issue comes up.
	//	Also, to be paranoid, I will check to make sure that the type are of valid-size at compile time.
	//
	//	The issues is that the standard really only gauentees that a type will have at least this many bits.. and on some exocit system where
	//		there are 3 bits to a byte (a byte is a collection of N bits where N is defines by how many bits are in the 'char' data type), an exactly
	//		64-bit number will not exist in that implelmentation. But that is a corner case we can ignore for now.
	//
	// As written, these typedefs work on Linux (64bit), and windows (32 bit and 64 bit). If you need to update for a new/exocit arhitecture,
	//	find out if the system is ILP64,LP64, or LLP64 and adjust accordingly.
	//	Example: 	Linux is an LP64: 		int=32, long=64, pointer=64, long long=64
	//						Windows is an LLP64:	int=32, long=64, pointer=64, long long=64
	
	// uintptr_t					// unsigned integer that can hold a pointer (basically the int type for void*) (on 64-bit system, this is 64 bits)
	// uintptr_t					// int big enough to hold the size (in bytes) of memoery allocations. (on 64-bit system, this is 64 bits)
	// size_t: That which is returned by sizeof() (and malloc etc..).. datatype for holding size of objects. (usualy unsigned int or unsigned long etc..)
	// ptrdiff_t: sogned type that comes from subtracting 2 pointers. (usually int, long etc..)
	typedef long long integer64;
	typedef unsigned long long uinteger64;

	#ifndef P_X_DIM
	    #error You Must Define P_X_DIMM to be 1,3, or 3.. Right now, it is not even defined (at least c/c++ cant see it!)
	#endif
	
	// The number of spatial components.
	#if P_X_DIM == 1
		#define p_x_dim 1
	#elif P_X_DIM == 2	
		#define p_x_dim 2
	#elif P_X_DIM == 3
		#define p_x_dim 3
	#else
		#error You Must Define P_X_DIM to be 1,3, or 3.. It is set to an invalid value.
	#endif
	
	// The number of velocity(momentum) components.
	#define p_max_dim 3


	#ifdef PRECISION_SINGLE
		#define particle_float_type   float
		#define float__part           float
		#define float__fld float
		#define particle_int_type     int

		#define field_float_type      float
		#define field_int_type        int
		
		#define	p_k_fld               float
		#define p_k_part              float  
		#define p_k_jay               float
		
		// vector types used in CUDA code
		#define float__part2 float##2
		#define float__part3 float##3
		#define float__part4 float##4
		
		// define floating type used in MPI
		#define FLOAT_MPI_TYPE MPI_FLOAT

	#elif PRECISION_DOUBLE

	  #define particle_float_type   double
	  #define float__part           double
	  #define float__fld            double
	  #define particle_int_type     int
		
	  #define field_float_type      double
	  #define field_int_type        int
		
	  #define	p_k_fld               double
	  #define p_k_part              double
	  #define p_k_jay               double
	  
	  // vector types used in CUDA code
	  #define float__part2 double##2
	  #define float__part3 double##3
	  #define float__part4 double##4

		// define floating type used in MPI
		#define FLOAT_MPI_TYPE MPI_DOUBLE

	#else
		#error You must specify single percision (PRECISION_SINGLE) or double percision (PRECISION_DOUBLE) compiler flags.
	#endif

	// misc constants that, for continitity and convience, mirror those in Osiris
	//	While it looks like many of these constants could/should be declared as c/c++ 'enum' (enumerated types),
	//	they are implmented as #defines.. The reason for this choice is that it is functionaly the same thing given how
	//	we use these cosntants. Also, many of the constants are meant to be array indices or ealuate to ints which enums are not
	//	really alloweed to do.

	// Boundry types... These must match the types and values seen in the Fortran Osiris code.
	#define p_bc_other_node -1
	#define p_bc_periodic 0
	#define p_bc_absorbing 5
	#define p_bc_conducting 5
	#define p_bc_specular 40
	#define p_bc_thermal 50
	#define p_bc_move_c 1

	#define p_lower 0
	#define p_upper 1



// TODO: Split This off into a different Include

	#if p_x_dim == 1
		#define _IX_  0
		#define _X_   1						//p_x_dim
		#define _P1_  2						//p_x_dim+p_x_dim
		#define _P2_  3						//p_x_dim+p_x_dim+1
		#define _P3_  4						//p_x_dim+p_x_dim+2
		#define _Q_   5						//p_x_dim+p_x_dim+p_p_dim

		// TODO.. take this out... just here for now...
		#define _IY_  0
		#define _Y_  1
		#define _IZ_ 0
		#define _Z_ 1
		#define PHASESPACE_SIZE 6
		
	#elif p_x_dim == 2
		#define _IX_  0
		#define _IY_  1
		#define _X_   2						//p_x_dim
		#define _Y_   3						//p_x_dim+1
		#define _P1_  4						//p_x_dim+p_x_dim
		#define _P2_  5						//p_x_dim+p_x_dim+1
		#define _P3_  6						//p_x_dim+p_x_dim+2
		#define _Q_   7						//p_x_dim+p_x_dim+p_p_dim

		#define PHASESPACE_SIZE 8
		
		// TODO: take out
		#define _IZ_  1
		#define _Z_  3
		
	#elif p_x_dim == 3
		#define _IX_  0
		#define _IY_  1
		#define _IZ_  2
		#define _X_   3						//p_x_dim
		#define _Y_   4						//p_x_dim+1
		#define _Z_   5
		#define _P1_  6						//p_x_dim+p_x_dim
		#define _P2_  7						//p_x_dim+p_x_dim+1
		#define _P3_  8						//p_x_dim+p_x_dim+2
		#define _Q_   9						//p_x_dim+p_x_dim+p_p_dim
		
		#define PHASESPACE_SIZE 10
		
	#endif

#define _P_   p_x_dim+p_x_dim

#endif
