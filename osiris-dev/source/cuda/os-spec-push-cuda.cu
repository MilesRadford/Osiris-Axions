// if __TEMPLATE__ is defined then read the template definition at the end of the file
#ifndef __TEMPLATE__

#include <stdio.h>
#include <inttypes.h>
#include <cuda_runtime.h>

#include "os-spec-push-cuda.h"

// ---------------------------------------------------------------------------------------
// Atomic add for architectures which don't support atomic adds for floats
// ---------------------------------------------------------------------------------------
#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 600)
__device__ double atomicAdd(double* address, double val) {
    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                               __longlong_as_double(assumed)));

    // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);

    return __longlong_as_double(old);
}
#endif

// ---------------------------------------------------------------------------------------
// Perform exclusive parallel prefix sum on p_double array idata which is of size BLOCK_SIZE
// Note: Could implement avoiding bank conficts if this kernel is heavy
// ---------------------------------------------------------------------------------------
template <int BLOCK_SIZE>
__device__ void prescan(p_double *const idata) {

  int tx = threadIdx.x;
  int offset = 1;

  // build sum in place up the tree
  for (int d = BLOCK_SIZE>>1; d > 0; d >>= 1) {
    __syncthreads();
    if (tx < d) {
      int ai = offset*(2*tx+1)-1;
      int bi = offset*(2*tx+2)-1;
      idata[bi] += idata[ai];
    }
    offset *= 2;
  }

  // clear the last element
  if (tx == 0) {
    idata[BLOCK_SIZE - 1] = 0;
  }

  // traverse down tree & build scan
  for (int d = 1; d < BLOCK_SIZE; d *= 2) {
    offset >>= 1;
    __syncthreads();
    if (tx < d) {
      int ai = offset*(2*tx+1)-1;
      int bi = offset*(2*tx+2)-1;
      p_double t = idata[ai];
      idata[ai] = idata[bi];
      idata[bi] += t;
    }
  }

}

//----------------------------------------------------------------------------------------
// Splines
//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
// Linear
//----------------------------------------------------------------------------------------
__forceinline__ __device__ void spline_s1( const p_k_fld x, p_k_fld *const s ) {

  s[0] = _0_5_p_k_fld - x;
  s[1] = _0_5_p_k_fld + x;

}

__forceinline__ __device__ void splineh_s1( const p_k_fld x, const int h,
                                            p_k_fld *const s ) {

  s[0] = ( 1 - h ) - x;
  s[1] = (     h ) + x;

}

//----------------------------------------------------------------------------------------
// Quadratic
//----------------------------------------------------------------------------------------
__forceinline__ __device__ void spline_s2( const p_k_fld x, p_k_fld *const s ) {

  const p_k_fld t0 = _0_5_p_k_fld - x;
  const p_k_fld t1 = _0_5_p_k_fld + x;

  s[-1] = _0_5_p_k_fld * t0*t0;
  s[ 0] = _0_5_p_k_fld + t0*t1;
  s[ 1] = _0_5_p_k_fld * t1*t1;

}

__forceinline__ __device__ void splineh_s2( const p_k_fld x, const int h,
                                            p_k_fld *const s ) {

  const p_k_fld t0 = ( 1 - h ) - x;
  const p_k_fld t1 = (     h ) + x;

  s[-1] = _0_5_p_k_fld * t0*t0;
  s[ 0] = _0_5_p_k_fld + t0*t1;
  s[ 1] = _0_5_p_k_fld * t1*t1;

}

//----------------------------------------------------------------------------------------
// Cubic
//----------------------------------------------------------------------------------------
__forceinline__ __device__ void spline_s3( const p_k_fld x, p_k_fld *const s ) {

  const p_k_fld c1_6 = _1_0_p_k_fld / p_k_fld(6.0);
  const p_k_fld c2_3 = _2_0_p_k_fld / p_k_fld(3.0);

  p_k_fld t0 = _0_5_p_k_fld - x;
  p_k_fld t1 = _0_5_p_k_fld + x;

  const p_k_fld t2 = t0 * t0;
  const p_k_fld t3 = t1 * t1;

  t0 *= t2;
  t1 *= t3;

  s[-1] = c1_6 * t0;
  s[ 0] = ( c2_3 - t3 ) + _0_5_p_k_fld*t1;
  s[ 1] = ( c2_3 - t2 ) + _0_5_p_k_fld*t0;
  s[ 2] = c1_6 * t1;

}

__forceinline__ __device__ void splineh_s3( const p_k_fld x, const int h,
                                            p_k_fld *const s ) {

  const p_k_fld c1_6 = _1_0_p_k_fld / p_k_fld(6.0);
  const p_k_fld c2_3 = _2_0_p_k_fld / p_k_fld(3.0);

  p_k_fld t0 = ( 1 - h ) - x;
  p_k_fld t1 = (     h ) + x;

  const p_k_fld t2 = t0 * t0;
  const p_k_fld t3 = t1 * t1;

  t0 *= t2;
  t1 *= t3;

  s[-1] = c1_6 * t0;
  s[ 0] = ( c2_3 - t3 ) + _0_5_p_k_fld*t1;
  s[ 1] = ( c2_3 - t2 ) + _0_5_p_k_fld*t0;
  s[ 2] = c1_6 * t1;

}

//----------------------------------------------------------------------------------------
// Quartic
//----------------------------------------------------------------------------------------
__forceinline__ __device__ void spline_s4( const p_k_fld x, p_k_fld *const s ) {

  const p_k_fld c1_6   = _1_0_p_k_fld / p_k_fld(6.0);
  const p_k_fld c1_24  = _1_0_p_k_fld / p_k_fld(24.0);
  const p_k_fld c11_24 = p_k_fld(11.0) / p_k_fld(24.0);

  const p_k_fld t0 = _0_5_p_k_fld - x;
  const p_k_fld t1 = _0_5_p_k_fld + x;
  const p_k_fld t2 = t0 * t1;

  s[-2] = c1_24 * ((t0*t0)*(t0*t0));
  s[-1] = c1_6 * ( p_k_fld(0.25) + t0 * ( _1_0_p_k_fld + t0 * ( p_k_fld(1.5) + t2 ) ) );
  s[ 0] = c11_24 + t2 * ( _0_5_p_k_fld + p_k_fld(0.25) * t2 );
  s[ 1] = c1_6 * ( p_k_fld(0.25) + t1 * ( _1_0_p_k_fld + t1 * ( p_k_fld(1.5) + t2 ) ) );
  s[ 2] = c1_24 * ((t1*t1)*(t1*t1));

}

__forceinline__ __device__ void splineh_s4( const p_k_fld x, const int h,
                                            p_k_fld *const s ) {

  const p_k_fld c1_6   = _1_0_p_k_fld / p_k_fld(6.0);
  const p_k_fld c1_24  = _1_0_p_k_fld / p_k_fld(24.0);
  const p_k_fld c11_24 = p_k_fld(11.0) / p_k_fld(24.0);

  const p_k_fld t0 = ( 1 - h ) - x;
  const p_k_fld t1 = (     h ) + x;
  const p_k_fld t2 = t0 * t1;

  s[-2] = c1_24 * ((t0*t0)*(t0*t0));
  s[-1] = c1_6 * ( p_k_fld(0.25) + t0 * ( _1_0_p_k_fld + t0 * ( p_k_fld(1.5) + t2 ) ) );
  s[ 0] = c11_24 + t2 * ( _0_5_p_k_fld + p_k_fld(0.25) * t2 );
  s[ 1] = c1_6 * ( p_k_fld(0.25) + t1 * ( _1_0_p_k_fld + t1 * ( p_k_fld(1.5) + t2 ) ) );
  s[ 2] = c1_24 * ((t1*t1)*(t1*t1));

}

//----------------------------------------------------------------------------------------
// Signbit function
//----------------------------------------------------------------------------------------
__forceinline__ __device__ int signbit_custom( const p_k_fld x ) {

  if ( x < 0 ) {
    return 1;
  } else {
    return 0;
  }

}

//----------------------------------------------------------------------------------------
// Longitudinal current weights
//----------------------------------------------------------------------------------------
__forceinline__ __device__ void wl_s1( const p_k_fld qnx, const p_k_fld x0,
                                       const p_k_fld x1, p_k_fld *const wl ) {

  wl[0] = qnx * ( x1 - x0 );

}

__forceinline__ __device__ void wl_s2( const p_k_fld qnx, const p_k_fld x0,
                                       const p_k_fld x1, p_k_fld *const wl ) {

  const p_k_fld d    = x1 - x0;
  const p_k_fld s1_2 = _0_5_p_k_fld - x0;
  const p_k_fld p1_2 = _0_5_p_k_fld + x0;

  wl[-1] = s1_2 - _0_5_p_k_fld * d;
  wl[ 0] = p1_2 + _0_5_p_k_fld * d;

  const p_k_fld n = qnx * d;
  wl[-1] *= n;
  wl[ 0]  *= n;

}

__forceinline__ __device__ void wl_s3( const p_k_fld qnx, const p_k_fld x0,
                                       const p_k_fld x1, p_k_fld *const wl ) {

  const p_k_fld c1_3 = _1_0_p_k_fld / p_k_fld(3.0);
  const p_k_fld d    = x1 - x0;
  const p_k_fld s1_2 = _0_5_p_k_fld - x0;
  const p_k_fld p1_2 = _0_5_p_k_fld + x0;
  const p_k_fld d_3  = c1_3 * d;

  wl[-1] =  _0_5_p_k_fld * ( s1_2*s1_2 - d * ( s1_2 - d_3 ) );
  wl[ 0] =  ( p_k_fld(0.75) - x0*x0 ) - d * ( x0 + d_3 );
  wl[ 1] =  _0_5_p_k_fld * ( p1_2*p1_2 + d * ( p1_2 + d_3 ) );

  const p_k_fld n = qnx * d;
  wl[-1] *= n;
  wl[ 0] *= n;
  wl[ 1] *= n;

}

__forceinline__ __device__ void wl_s4( const p_k_fld qnx, const p_k_fld x0,
                                       const p_k_fld x1, p_k_fld *const wl ) {

  const p_k_fld c1_6 = _1_0_p_k_fld / p_k_fld(6.0);
  const p_k_fld c1_3 = _1_0_p_k_fld / p_k_fld(3.0);
  const p_k_fld c2_3 = _2_0_p_k_fld / p_k_fld(3.0);

  const p_k_fld d = x1 - x0;
  const p_k_fld s = _0_5_p_k_fld - x0;
  const p_k_fld p = _0_5_p_k_fld + x0;
  const p_k_fld t = x0 - c1_6;
  const p_k_fld u = x0 + c1_6;

  const p_k_fld s2 = s * s;
  const p_k_fld s3 = s2 * s;
  const p_k_fld p2 = p * p;
  const p_k_fld p3 = p2 * p;

  wl[-2] =  c1_6 * ( d * ( d * ( s - p_k_fld(0.25)*d ) - p_k_fld(1.5) * s2 ) + s3 );
  wl[-1] =    d * ( ( _0_5_p_k_fld * d ) * (t + p_k_fld(0.25)*d) + ( p_k_fld(0.75)*(t*t) - c1_3 ) ) +
            ( _0_5_p_k_fld * p3 + ( c2_3 - p2 ) );
  wl[ 0] =  - d * ( ( _0_5_p_k_fld * d ) * (u + p_k_fld(0.25)*d) + ( p_k_fld(0.75)*(u*u) - c1_3 ) ) +
            ( _0_5_p_k_fld * s3 + ( c2_3 - s2 ) );
  wl[ 1] =  c1_6 * ( d * ( d * ( p + p_k_fld(0.25)*d ) + p_k_fld(1.5) * p2 ) + p3 );

  const p_k_fld n = qnx * d;
  wl[-2] *= n;
  wl[-1] *= n;
  wl[ 0] *= n;
  wl[ 1] *= n;

}


//----------------------------------------------------------------------------------------
// Returns the integer shift (-1, 0 or +1) so that the coordinate remains in the [-0.5, 0.5]
// range. This is the fastest implementation (twice as fast as a sequence of ifs) because
// the two if structures compile as conditional moves and can be processed independently.
// This has no precision problem and is only 12% slower than the previous "int(x+1.5)-1"
// routine that would break for x = nearest( 0.5, -1.0 )
//----------------------------------------------------------------------------------------
__forceinline__ __device__ int ntrim( const p_k_part x ) {
  int a, b;

  if ( x < -.5 ) {
    a = -1;
  } else {
    a = 0;
  }

  if ( x >= .5 ) {
    b = +1;
  } else {
    b = 0;
  }

  return a + b;

}


//----------------------------------------------------------------------------------------
// Zero out the current array
//----------------------------------------------------------------------------------------
template <int BLOCK_SIZE>
__global__ void zero_jay_1d( const int *const prefix_size, const int avg_block_tile, const t_f1<p_k_fld> *const jay_f1){

  const int tx = threadIdx.x;
  const int bx = blockIdx.x;

  __shared__ int til, ind_start, size;
  __shared__ t_f1<p_k_fld> jay_f1_til;

  if (tx==0) {
    til = bx / avg_block_tile; // floor, best guess, could be wrong

    // Find the correct tile for this block number
    if ( bx < prefix_size[til] ) {
      while( bx < prefix_size[til] ) {
        til--;
      }
    } else {
      while( bx >= prefix_size[til+1] ) {
        til++;
      }
    }

    ind_start = (bx-prefix_size[til]) * BLOCK_SIZE;

    jay_f1_til = jay_f1[til];

    size = jay_f1_til.size();

  }
  __syncthreads();

  if (tx+ind_start < size) jay_f1_til.f[tx+ind_start] = 0;

}


//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
template <int BLOCK_SIZE>
__global__ void zero_jay_2d( const int *const prefix_size, const int avg_block_tile, const t_f2<p_k_fld> *const jay_f2){

  const int tx = threadIdx.x;
  const int bx = blockIdx.x;

  __shared__ int til, ind_start, size;
  __shared__ t_f2<p_k_fld> jay_f2_til;

  if (tx==0) {
    til = bx / avg_block_tile; // floor, best guess, could be wrong

    // Find the correct tile for this block number
    if ( bx < prefix_size[til] ) {
      while( bx < prefix_size[til] ) {
        til--;
      }
    } else {
      while( bx >= prefix_size[til+1] ) {
        til++;
      }
    }

    ind_start = (bx-prefix_size[til]) * BLOCK_SIZE;

    jay_f2_til = jay_f2[til];

    size = jay_f2_til.size();

  }
  __syncthreads();

  if (tx+ind_start < size) jay_f2_til.f[tx+ind_start] = 0;

}


//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
template <int BLOCK_SIZE>
__global__ void zero_jay_3d( const int *const prefix_size, const int avg_block_tile, const t_f3<p_k_fld> *const jay_f3){

  const int tx = threadIdx.x;
  const int bx = blockIdx.x;

  __shared__ int til, ind_start, size;
  __shared__ t_f3<p_k_fld> jay_f3_til;

  if (tx==0) {
    til = bx / avg_block_tile; // floor, best guess, could be wrong

    // Find the correct tile for this block number
    if ( bx < prefix_size[til] ) {
      while( bx < prefix_size[til] ) {
        til--;
      }
    } else {
      while( bx >= prefix_size[til+1] ) {
        til++;
      }
    }

    ind_start = (bx-prefix_size[til]) * BLOCK_SIZE;

    jay_f3_til = jay_f3[til];

    size = jay_f3_til.size();

  }
  __syncthreads();

  if (tx+ind_start < size) jay_f3_til.f[tx+ind_start] = 0;

}


//----------------------------------------------------------------------------------------
// Function stub to call kernel
//----------------------------------------------------------------------------------------
void zero_jay_wrapper( const int *const prefix_size, const int avg_block_tile, 
                       const t_current_c *const jay, const int nblocks,
                       const cudaStream_t stream ){

  switch (p_x_dim) {
    case 1:
      zero_jay_1d<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( prefix_size, avg_block_tile,
                                                              jay->d_jd1 );
      break;
    case 2:
      zero_jay_2d<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( prefix_size, avg_block_tile,
                                                              jay->d_jd2 );
      break;
    case 3:
      zero_jay_3d<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( prefix_size, avg_block_tile,
                                                              jay->d_jd3 );
      break;
  }

    chkLastErr("zero_jay");
#ifdef __DEBUG__
    chkErr(cudaDeviceSynchronize());
#endif

}


//----------------------------------------------------------------------------------------
// Shift particle positions for moving window
//----------------------------------------------------------------------------------------
template <int BLOCK_SIZE>
__global__ void shift_positions( t_chunk ***chunks, int *block_start_tiles,
                                 int *block_start_chunks, int *block_num_par,
                                 int nmove0, int nmove1, int nmove2, int chunk_size ){

  const int tx = threadIdx.x;
  const int bx = blockIdx.x;

  // Shared variables that are the same between all threads
  __shared__ int til;

  if (tx==0) {
    til = block_start_tiles[bx];
  }
  __syncthreads();

  int num_par_proc_tot = block_num_par[bx];
  t_chunk **chunk = chunks[til] + block_start_chunks[bx];

  int nmove[3] = { nmove0, nmove1, nmove2 };
  int num_par_proc_chunk, part_idx;

  t_f1<int> ix = t_f1<int>( chunk_size, p_x_dim );

  const size_t ix_offset = ( p_x_dim + p_p_dim + 1 ) * sizeof( p_k_part ) * chunk_size;

  while( num_par_proc_tot > 0 ) {

    ix.g = (int *) ( *chunk + ix_offset);

    // loop over chunks of particles until all particles in tile have been processed
    num_par_proc_chunk = min(num_par_proc_tot, chunk_size);
    part_idx = tx;

    while ( part_idx < num_par_proc_chunk ){

      for (int i=0; i<p_x_dim; i++) {
        ix.fd(part_idx,i) -= nmove[i];
      }

      // Note: this only works if the chunks size is divisible by the block size (which we enforce upon initialization)
      // See sort for how to do it more generally if we decide we want to do that here too
      part_idx += BLOCK_SIZE;

    } // end loop over chunk of particles ``while ( part_idx < num_par )``

    num_par_proc_tot -= num_par_proc_chunk;
    chunk++;

  } // end loop over chunks

}


//----------------------------------------------------------------------------------------
// Function stub to call kernel
//----------------------------------------------------------------------------------------
void shift_positions_wrapper( t_chunk ***chunks, int *block_start_tiles,
                              int *block_start_chunks, int *block_num_par,
                              int nmove0, int nmove1, int nmove2,
                              int nblocks, int chunk_size ){

  shift_positions<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__ >>>( chunks, block_start_tiles,
                                                   block_start_chunks, block_num_par,
                                                   nmove0, nmove1, nmove2, chunk_size );

  chkLastErr("shift_positions");
#ifdef __DEBUG__
  chkErr(cudaDeviceSynchronize());
#endif

}


#define __TEMPLATE__

//********************************* Linear interpolation *********************************

#define DUDT_BORIS_1D   dudt_boris_1d_s1
#define DUDT_BORIS_2D   dudt_boris_2d_s1
#define DUDT_BORIS_3D   dudt_boris_3d_s1

#define DUDT_BORIS_1D_CS   dudt_boris_1d_cs1
#define DUDT_BORIS_2D_CS   dudt_boris_2d_cs1
#define DUDT_BORIS_3D_CS   dudt_boris_3d_cs1

#define SPLINE spline_s1
#define SPLINEH splineh_s1

// Lower/upper point limits for dudt
#define LP 0
#define UP 1

#define WL wl_s1

// Perpendicular current cell  limits
#define P0 0
#define P1 1

#define ADVANCE_DEPOSIT_1D   advance_deposit_1d_s1
#define ADVANCE_DEPOSIT_2D   advance_deposit_2d_s1
#define ADVANCE_DEPOSIT_3D   advance_deposit_3d_s1

#include __FILE__

//******************************* Quadratic interpolation ********************************

#define DUDT_BORIS_1D   dudt_boris_1d_s2
#define DUDT_BORIS_2D   dudt_boris_2d_s2
#define DUDT_BORIS_3D   dudt_boris_3d_s2

#define DUDT_BORIS_1D_CS   dudt_boris_1d_cs2
#define DUDT_BORIS_2D_CS   dudt_boris_2d_cs2
#define DUDT_BORIS_3D_CS   dudt_boris_3d_cs2

#define SPLINE spline_s2
#define SPLINEH splineh_s2

// Lower/upper point limits for dudt
#define LP -1
#define UP 1

#define WL wl_s2

// Perpendicular current cell  limits
#define P0 -1
#define P1 1

#define ADVANCE_DEPOSIT_1D   advance_deposit_1d_s2
#define ADVANCE_DEPOSIT_2D   advance_deposit_2d_s2
#define ADVANCE_DEPOSIT_3D   advance_deposit_3d_s2

#include __FILE__

//********************************* Cubic interpolation **********************************

#define DUDT_BORIS_1D   dudt_boris_1d_s3
#define DUDT_BORIS_2D   dudt_boris_2d_s3
#define DUDT_BORIS_3D   dudt_boris_3d_s3

#define DUDT_BORIS_1D_CS   dudt_boris_1d_cs3
#define DUDT_BORIS_2D_CS   dudt_boris_2d_cs3
#define DUDT_BORIS_3D_CS   dudt_boris_3d_cs3

#define SPLINE spline_s3
#define SPLINEH splineh_s3

// Lower/upper point limits for dudt
#define LP -1
#define UP 2

#define WL wl_s3

// Perpendicular current cell  limits
#define P0 -1
#define P1 2

#define ADVANCE_DEPOSIT_1D   advance_deposit_1d_s3
#define ADVANCE_DEPOSIT_2D   advance_deposit_2d_s3
#define ADVANCE_DEPOSIT_3D   advance_deposit_3d_s3

#include __FILE__

//******************************** Quartic interpolation *********************************

#define DUDT_BORIS_1D   dudt_boris_1d_s4
#define DUDT_BORIS_2D   dudt_boris_2d_s4
#define DUDT_BORIS_3D   dudt_boris_3d_s4

#define DUDT_BORIS_1D_CS   dudt_boris_1d_cs4
#define DUDT_BORIS_2D_CS   dudt_boris_2d_cs4
#define DUDT_BORIS_3D_CS   dudt_boris_3d_cs4

#define SPLINE spline_s4
#define SPLINEH splineh_s4

// Lower/upper point limits for dudt
#define LP -2
#define UP 2

#define WL wl_s4

// Perpendicular current cell  limits
#define P0 -2
#define P1 2

#define ADVANCE_DEPOSIT_1D   advance_deposit_1d_s4
#define ADVANCE_DEPOSIT_2D   advance_deposit_2d_s4
#define ADVANCE_DEPOSIT_3D   advance_deposit_3d_s4

#include __FILE__

#undef __TEMPLATE__

#else


// ---------------------------------------------------------------------------------------
// Advance velocities and calculate time-centered total energy
// with fields on a staggered Yee mesh
// ---------------------------------------------------------------------------------------
template <int BLOCK_SIZE>
__global__ void DUDT_BORIS_1D( t_chunk ***chunks, int *block_start_tiles,
                               int *block_start_chunks, int *block_num_par,
                               p_k_part rqm, t_f1<p_k_fld> *b_f1,
                               t_f1<p_k_fld> *e_f1, p_double *energy,
                               bool report_energy, p_double dt, int chunk_size ){

  const int bx = blockIdx.x;
  const int tx = threadIdx.x;

  // Shared variables that are the same between all threads
  __shared__ int til;
  __shared__ p_k_part tem;
  __shared__ t_f1<p_k_fld> b_f1_til, e_f1_til;

  // Shared variables that will be prefix scanned
  __shared__ p_double accum_energy[BLOCK_SIZE];

  if (tx==0) {
    til = block_start_tiles[bx];

    b_f1_til = b_f1[til];
    e_f1_til = e_f1[til];

    // get factor that includes timestep and charge-to-mass ratio
    // note that the charge-to-mass ratio is used since the
    // momentum p is not the total momentum of a particles but
    // the momentum per unit restmass (which is the electron mass)
    tem = p_k_part( _0_5_p_double * dt / rqm );
  }
  __syncthreads();

  int num_par_proc_tot = block_num_par[bx];
  t_chunk **chunk = chunks[til] + block_start_chunks[bx];

  int i1, i1h, h1, num_par_proc_chunk, part_idx;

  p_k_part u2, gamma, gam_tem, otsq;
  p_k_part bp[p_p_dim], ep[p_p_dim], utemp[p_p_dim];

  p_k_fld dx1;
  p_k_fld f1, f2, f3;

  p_k_fld w1_buff[UP-LP+1], w1h_buff[UP-LP+1];

  p_k_fld *const w1  = w1_buff  - LP;
  p_k_fld *const w1h = w1h_buff - LP;

  t_f1<p_k_part> x = t_f1<p_k_part>( chunk_size, p_x_dim );
  t_f1<p_k_part> p = t_f1<p_k_part>( chunk_size, p_p_dim );
  t_f0<p_k_part> q = t_f0<p_k_part>( chunk_size );
  t_f1<int>     ix = t_f1<int>( chunk_size, p_x_dim );

  const size_t chunk_size_t = sizeof( p_k_part ) * chunk_size;
  const size_t x_offset = 0;
  const size_t p_offset = x_offset + p_x_dim * chunk_size_t;
  const size_t q_offset = p_offset +  p_p_dim * chunk_size_t;
  const size_t ix_offset = q_offset + chunk_size_t;

  p_double accum_energy_last, loc_ene = 0;

  while( num_par_proc_tot > 0 ) {

    x.g  = (p_k_part *) ( *chunk + x_offset );
    p.g  = (p_k_part *) ( *chunk + p_offset );
    q.g  = (p_k_part *) ( *chunk + q_offset );
    ix.g = (int *)      ( *chunk + ix_offset);

    // modify bp & ep to include timestep and charge-to-mass ratio
    // and perform half the electric field acceleration.
    // Result is stored in UTEMP.

    // loop over chunks of particles until all particles in tile have been processed
    num_par_proc_chunk = min(num_par_proc_tot, chunk_size);
    part_idx = tx;

    while ( part_idx < num_par_proc_chunk ){

      // --------------
      // dudt_boris()
      // --------------
      i1 = ix.fd(part_idx,0) - 1;
      dx1 = p_k_fld( x.fd(part_idx,0) );

      h1 = signbit_custom(dx1);
      i1h = i1 - h1;

      // get spline weitghts for x
      SPLINE( dx1, w1 );
      SPLINEH( dx1, h1, w1h );

      // Interpolate E - Field
      f1 = 0;
      f2 = 0;
      f3 = 0;

      for (int k1=LP; k1<=UP; ++k1){
        f1 += e_f1_til.fd(0,i1h+k1) * w1h[k1];
        f2 += e_f1_til.fd(1,i1 +k1) * w1[ k1];
        f3 += e_f1_til.fd(2,i1 +k1) * w1[ k1];
      }

      ep[0] = f1;
      ep[1] = f2;
      ep[2] = f3;

      // Interpolate B - Field
      f1 = 0;
      f2 = 0;
      f3 = 0;

      for (int k1=LP; k1<=UP; ++k1){
        f1 += b_f1_til.fd(0,i1 +k1) * w1[ k1];
        f2 += b_f1_til.fd(1,i1h+k1) * w1h[k1];
        f3 += b_f1_til.fd(2,i1h+k1) * w1h[k1];
      }

      bp[0] = f1;
      bp[1] = f2;
      bp[2] = f3;

      // --------------

      ep[0] *= tem;
      ep[1] *= tem;
      ep[2] *= tem;

      utemp[0] = p.fd(part_idx,0) + ep[0];
      utemp[1] = p.fd(part_idx,1) + ep[1];
      utemp[2] = p.fd(part_idx,2) + ep[2];

      // Get time centered gamma
      u2 = utemp[0]*utemp[0] + utemp[1]*utemp[1] + utemp[2]*utemp[2];

      gamma = sqrt(u2+1);

      // accumulate time centered energy
      // this is done in double precicion always
      if (report_energy) loc_ene += q.fd(part_idx) * u2 / (gamma + _1_0_p_double);

      gam_tem = tem / gamma;

      bp[0] *= gam_tem;
      bp[1] *= gam_tem;
      bp[2] *= gam_tem;

      p.fd(part_idx,0) = utemp[0] + utemp[1] * bp[2];
      p.fd(part_idx,1) = utemp[1] + utemp[2] * bp[0];
      p.fd(part_idx,2) = utemp[2] + utemp[0] * bp[1];

      p.fd(part_idx,0) -= utemp[2] * bp[1];
      p.fd(part_idx,1) -= utemp[0] * bp[2];
      p.fd(part_idx,2) -= utemp[1] * bp[0];

      otsq = _2_0_p_k_part / ( ( (_1_0_p_k_part + bp[0]*bp[0]) + bp[1]*bp[1] ) + bp[2]*bp[2] );

      bp[0] *= otsq;
      bp[1] *= otsq;
      bp[2] *= otsq;

      utemp[0] += p.fd(part_idx,1) * bp[2];
      utemp[1] += p.fd(part_idx,2) * bp[0];
      utemp[2] += p.fd(part_idx,0) * bp[1];

      utemp[0] -= p.fd(part_idx,2) * bp[1];
      utemp[1] -= p.fd(part_idx,0) * bp[2];
      utemp[2] -= p.fd(part_idx,1) * bp[0];

      // Perform second half of electric field acceleration.
      p.fd(part_idx,0) = utemp[0] + ep[0];
      p.fd(part_idx,1) = utemp[1] + ep[1];
      p.fd(part_idx,2) = utemp[2] + ep[2];

      // Note, this only works if the chunks size is divisible by the block size (which
      // we enforce upon initialization)
      // See sort for how to do it more generally if we decide we want to do that here too
      part_idx += BLOCK_SIZE;

    } // end loop over chunk of particles ``while ( part_idx < num_par )``

    num_par_proc_tot -= num_par_proc_chunk;
    chunk++;

  } // end loop over chunks

  if (report_energy) {
    accum_energy[tx] = loc_ene;
    if (tx==0) accum_energy_last = accum_energy[BLOCK_SIZE-1];
    prescan<BLOCK_SIZE>( accum_energy );
    __syncthreads();
    if (tx==0) {
      atomicAdd( energy+til, accum_energy[BLOCK_SIZE-1] + accum_energy_last );
    }
  }

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
template <int BLOCK_SIZE>
__global__ void DUDT_BORIS_2D( t_chunk ***chunks, int *block_start_tiles,
                               int *block_start_chunks, int *block_num_par,
                               p_k_part rqm, t_f2<p_k_fld> *b_f2,
                               t_f2<p_k_fld> *e_f2, p_double *energy,
                               bool report_energy, p_double dt, int chunk_size ){

  const int bx = blockIdx.x;
  const int tx = threadIdx.x;

  // Shared variables that are the same between all threads
  __shared__ int til;
  __shared__ p_k_part tem;
  __shared__ t_f2<p_k_fld> b_f2_til, e_f2_til;

  // Shared variables that will be prefix scanned
  __shared__ p_double accum_energy[BLOCK_SIZE];

  if (tx==0) {
    til = block_start_tiles[bx];

    b_f2_til = b_f2[til];
    e_f2_til = e_f2[til];

    // get factor that includes timestep and charge-to-mass ratio
    // note that the charge-to-mass ratio is used since the
    // momentum p is not the total momentum of a particles but
    // the momentum per unit restmass (which is the electron mass)
    tem = p_k_part( _0_5_p_double * dt / rqm );
  }
  __syncthreads();

  int num_par_proc_tot = block_num_par[bx];
  t_chunk **chunk = chunks[til] + block_start_chunks[bx];

  int i1, i2, i1h, i2h, h1, h2, num_par_proc_chunk, part_idx;

  p_k_part u2, gamma, gam_tem, otsq;
  p_k_part bp[p_p_dim], ep[p_p_dim], utemp[p_p_dim];

  p_k_fld dx1, dx2;
  p_k_fld f1, f2, f3, f1line, f2line, f3line;

  p_k_fld w1_buff[UP-LP+1], w2_buff[UP-LP+1], w1h_buff[UP-LP+1], w2h_buff[UP-LP+1];

  p_k_fld *const w1  = w1_buff  - LP;
  p_k_fld *const w2  = w2_buff  - LP;
  p_k_fld *const w1h = w1h_buff - LP;
  p_k_fld *const w2h = w2h_buff - LP;

  t_f1<p_k_part> x = t_f1<p_k_part>( chunk_size, p_x_dim );
  t_f1<p_k_part> p = t_f1<p_k_part>( chunk_size, p_p_dim );
  t_f0<p_k_part> q = t_f0<p_k_part>( chunk_size );
  t_f1<int>     ix = t_f1<int>( chunk_size, p_x_dim );

  const size_t chunk_size_t = sizeof( p_k_part ) * chunk_size;
  const size_t x_offset = 0;
  const size_t p_offset = x_offset + p_x_dim * chunk_size_t;
  const size_t q_offset = p_offset +  p_p_dim * chunk_size_t;
  const size_t ix_offset = q_offset + chunk_size_t;

  p_double accum_energy_last, loc_ene = 0;

  while( num_par_proc_tot > 0 ) {

    x.g  = (p_k_part *) ( *chunk + x_offset );
    p.g  = (p_k_part *) ( *chunk + p_offset );
    q.g  = (p_k_part *) ( *chunk + q_offset );
    ix.g = (int *)      ( *chunk + ix_offset);

    // modify bp & ep to include timestep and charge-to-mass ratio
    // and perform half the electric field acceleration.
    // Result is stored in UTEMP.

    // loop over chunks of particles until all particles in tile have been processed
    num_par_proc_chunk = min(num_par_proc_tot, chunk_size);
    part_idx = tx;

    while ( part_idx < num_par_proc_chunk ){

      // --------------
      // dudt_boris()
      // --------------
      i1 = ix.fd(part_idx,0) - 1;
      i2 = ix.fd(part_idx,1) - 1;

      dx1 = p_k_fld( x.fd(part_idx,0) );
      dx2 = p_k_fld( x.fd(part_idx,1) );

      h1 = signbit_custom(dx1);
      h2 = signbit_custom(dx2);

      i1h = i1 - h1;
      i2h = i2 - h2;

      // get spline weitghts for x and y
      SPLINE( dx1, w1 );
      SPLINEH( dx1, h1, w1h );

      SPLINE( dx2, w2 );
      SPLINEH( dx2, h2, w2h );

      // Interpolate E - Field
      f1 = 0;
      f2 = 0;
      f3 = 0;

      for (int k2=LP; k2<=UP; ++k2){
        f1line = 0;
        f2line = 0;
        f3line = 0;

        for (int k1=LP; k1<=UP; ++k1){
          f1line += e_f2_til.fd(0,i1h+k1,i2 +k2) * w1h[k1];
          f2line += e_f2_til.fd(1,i1 +k1,i2h+k2) * w1[ k1];
          f3line += e_f2_til.fd(2,i1 +k1,i2 +k2) * w1[ k1];
        }

        f1 += f1line * w2[ k2];
        f2 += f2line * w2h[k2];
        f3 += f3line * w2[ k2];
      }

      ep[0] = f1;
      ep[1] = f2;
      ep[2] = f3;

      // Interpolate B - Field
      f1 = 0;
      f2 = 0;
      f3 = 0;

      for (int k2=LP; k2<=UP; ++k2){
        f1line = 0;
        f2line = 0;
        f3line = 0;

        for (int k1=LP; k1<=UP; ++k1){
          f1line += b_f2_til.fd(0,i1 +k1,i2h+k2) * w1[ k1];
          f2line += b_f2_til.fd(1,i1h+k1,i2 +k2) * w1h[k1];
          f3line += b_f2_til.fd(2,i1h+k1,i2h+k2) * w1h[k1];
        }

        f1 += f1line * w2h[k2];
        f2 += f2line * w2[ k2];
        f3 += f3line * w2h[k2];
      }

      bp[0] = f1;
      bp[1] = f2;
      bp[2] = f3;

      // --------------

      ep[0] *= tem;
      ep[1] *= tem;
      ep[2] *= tem;

      utemp[0] = p.fd(part_idx,0) + ep[0];
      utemp[1] = p.fd(part_idx,1) + ep[1];
      utemp[2] = p.fd(part_idx,2) + ep[2];

      // Get time centered gamma
      u2 = utemp[0]*utemp[0] + utemp[1]*utemp[1] + utemp[2]*utemp[2];

      gamma = sqrt(u2+1);

      // accumulate time centered energy
      // this is done in double precicion always
      if (report_energy) loc_ene += q.fd(part_idx) * u2 / (gamma + _1_0_p_double);

      gam_tem = tem / gamma;

      bp[0] *= gam_tem;
      bp[1] *= gam_tem;
      bp[2] *= gam_tem;

      p.fd(part_idx,0) = utemp[0] + utemp[1] * bp[2];
      p.fd(part_idx,1) = utemp[1] + utemp[2] * bp[0];
      p.fd(part_idx,2) = utemp[2] + utemp[0] * bp[1];

      p.fd(part_idx,0) -= utemp[2] * bp[1];
      p.fd(part_idx,1) -= utemp[0] * bp[2];
      p.fd(part_idx,2) -= utemp[1] * bp[0];

      otsq = _2_0_p_k_part / ( ( (_1_0_p_k_part + bp[0]*bp[0]) + bp[1]*bp[1] ) + bp[2]*bp[2] );

      bp[0] *= otsq;
      bp[1] *= otsq;
      bp[2] *= otsq;

      utemp[0] += p.fd(part_idx,1) * bp[2];
      utemp[1] += p.fd(part_idx,2) * bp[0];
      utemp[2] += p.fd(part_idx,0) * bp[1];

      utemp[0] -= p.fd(part_idx,2) * bp[1];
      utemp[1] -= p.fd(part_idx,0) * bp[2];
      utemp[2] -= p.fd(part_idx,1) * bp[0];

      // Perform second half of electric field acceleration.
      p.fd(part_idx,0) = utemp[0] + ep[0];
      p.fd(part_idx,1) = utemp[1] + ep[1];
      p.fd(part_idx,2) = utemp[2] + ep[2];

      // Note: this only works if chunk size divisible by block size (which we enforce at initialization)
      // See sort for how to do this more generally, if we decide to do that eventually
      part_idx += BLOCK_SIZE;

    } // end loop over chunk of particles ``while ( part_idx < num_par )``

    num_par_proc_tot -= num_par_proc_chunk;
    chunk++;

  } // end loop over chunks

  if (report_energy) {
    accum_energy[tx] = loc_ene;
    if (tx==0) accum_energy_last = accum_energy[BLOCK_SIZE-1];
    prescan<BLOCK_SIZE>( accum_energy );
    __syncthreads();
    if (tx==0) {
      atomicAdd( energy+til, accum_energy[BLOCK_SIZE-1] + accum_energy_last );
    }
  }

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
template <int BLOCK_SIZE>
__global__ void DUDT_BORIS_3D( t_chunk ***chunks, int *block_start_tiles,
                               int *block_start_chunks, int *block_num_par,
                               p_k_part rqm, t_f3<p_k_fld> *b_f3,
                               t_f3<p_k_fld> *e_f3, p_double *energy,
                               bool report_energy, p_double dt, int chunk_size ){

  const int bx = blockIdx.x;
  const int tx = threadIdx.x;

  // Shared variables that are the same between all threads
  __shared__ int til;
  __shared__ p_k_part tem;
  __shared__ t_f3<p_k_fld> b_f3_til, e_f3_til;

  // Shared variables that will be prefix scanned
  __shared__ p_double accum_energy[BLOCK_SIZE];

  if (tx==0) {
    til = block_start_tiles[bx];

    b_f3_til = b_f3[til];
    e_f3_til = e_f3[til];

    // get factor that includes timestep and charge-to-mass ratio
    // note that the charge-to-mass ratio is used since the
    // momentum p is not the total momentum of a particles but
    // the momentum per unit restmass (which is the electron mass)
    tem = p_k_part( _0_5_p_double * dt / rqm );
  }
  __syncthreads();

  int num_par_proc_tot = block_num_par[bx];
  t_chunk **chunk = chunks[til] + block_start_chunks[bx];

  int i1, i2, i3, i1h, i2h, i3h, h1, h2, h3, num_par_proc_chunk, part_idx;

  p_k_part u2, gamma, gam_tem, otsq;
  p_k_part bp[p_p_dim], ep[p_p_dim], utemp[p_p_dim];

  p_k_fld dx1, dx2, dx3;
  p_k_fld f1, f2, f3, f1line, f2line, f3line, f1plane, f2plane, f3plane;

  p_k_fld  w1_buff[UP-LP+1],  w2_buff[UP-LP+1],  w3_buff[UP-LP+1];
  p_k_fld w1h_buff[UP-LP+1], w2h_buff[UP-LP+1], w3h_buff[UP-LP+1];

  p_k_fld *const w1  = w1_buff  - LP;
  p_k_fld *const w2  = w2_buff  - LP;
  p_k_fld *const w3  = w3_buff  - LP;
  p_k_fld *const w1h = w1h_buff - LP;
  p_k_fld *const w2h = w2h_buff - LP;
  p_k_fld *const w3h = w3h_buff - LP;

  t_f1<p_k_part> x = t_f1<p_k_part>( chunk_size, p_x_dim );
  t_f1<p_k_part> p = t_f1<p_k_part>( chunk_size, p_p_dim );
  t_f0<p_k_part> q = t_f0<p_k_part>( chunk_size );
  t_f1<int>     ix = t_f1<int>( chunk_size, p_x_dim );

  const size_t chunk_size_t = sizeof( p_k_part ) * chunk_size;
  const size_t x_offset = 0;
  const size_t p_offset = x_offset + p_x_dim * chunk_size_t;
  const size_t q_offset = p_offset +  p_p_dim * chunk_size_t;
  const size_t ix_offset = q_offset + chunk_size_t;

  p_double accum_energy_last, loc_ene = 0;

  while( num_par_proc_tot > 0 ) {

    x.g  = (p_k_part *) ( *chunk + x_offset );
    p.g  = (p_k_part *) ( *chunk + p_offset );
    q.g  = (p_k_part *) ( *chunk + q_offset );
    ix.g = (int *)      ( *chunk + ix_offset);

    // modify bp & ep to include timestep and charge-to-mass ratio
    // and perform half the electric field acceleration.
    // Result is stored in UTEMP.

    // loop over chunks of particles until all particles in tile have been processed
    num_par_proc_chunk = min(num_par_proc_tot, chunk_size);
    part_idx = tx;

    while ( part_idx < num_par_proc_chunk ){

      // --------------
      // dudt_boris()
      // --------------
      i1 = ix.fd(part_idx,0) - 1;
      i2 = ix.fd(part_idx,1) - 1;
      i3 = ix.fd(part_idx,2) - 1;

      dx1 = p_k_fld( x.fd(part_idx,0) );
      dx2 = p_k_fld( x.fd(part_idx,1) );
      dx3 = p_k_fld( x.fd(part_idx,2) );

      h1 = signbit_custom(dx1);
      h2 = signbit_custom(dx2);
      h3 = signbit_custom(dx3);

      i1h = i1 - h1;
      i2h = i2 - h2;
      i3h = i3 - h3;

      // get spline weitghts for x, y and z
      SPLINE( dx1, w1 );
      SPLINEH( dx1, h1, w1h );

      SPLINE( dx2, w2 );
      SPLINEH( dx2, h2, w2h );

      SPLINE( dx3, w3 );
      SPLINEH( dx3, h3, w3h );

      // Interpolate E - Field
      f1 = 0;
      f2 = 0;
      f3 = 0;

      for (int k3=LP; k3<=UP; ++k3){
        f1plane = 0;
        f2plane = 0;
        f3plane = 0;

        for (int k2=LP; k2<=UP; ++k2){
          f1line = 0;
          f2line = 0;
          f3line = 0;

          for (int k1=LP; k1<=UP; ++k1){
            f1line += e_f3_til.fd(0,i1h+k1,i2 +k2,i3 +k3) * w1h[k1];
            f2line += e_f3_til.fd(1,i1 +k1,i2h+k2,i3 +k3) * w1[ k1];
            f3line += e_f3_til.fd(2,i1 +k1,i2 +k2,i3h+k3) * w1[ k1];
          }

          f1plane += f1line * w2[ k2];
          f2plane += f2line * w2h[k2];
          f3plane += f3line * w2[ k2];
        }

        f1 += f1plane * w3[ k3];
        f2 += f2plane * w3[ k3];
        f3 += f3plane * w3h[k3];

      }

      ep[0] = f1;
      ep[1] = f2;
      ep[2] = f3;

      // Interpolate B - Field
      f1 = 0;
      f2 = 0;
      f3 = 0;

      for (int k3=LP; k3<=UP; ++k3){
        f1plane = 0;
        f2plane = 0;
        f3plane = 0;

        for (int k2=LP; k2<=UP; ++k2){
          f1line = 0;
          f2line = 0;
          f3line = 0;

          for (int k1=LP; k1<=UP; ++k1){
            f1line += b_f3_til.fd(0,i1 +k1,i2h+k2,i3h+k3) * w1[ k1];
            f2line += b_f3_til.fd(1,i1h+k1,i2 +k2,i3h+k3) * w1h[k1];
            f3line += b_f3_til.fd(2,i1h+k1,i2h+k2,i3 +k3) * w1h[k1];
          }

          f1plane += f1line * w2h[k2];
          f2plane += f2line * w2[ k2];
          f3plane += f3line * w2h[k2];
        }

        f1 += f1plane * w3h[k3];
        f2 += f2plane * w3h[k3];
        f3 += f3plane * w3[ k3];

      }

      bp[0] = f1;
      bp[1] = f2;
      bp[2] = f3;

      // --------------

      ep[0] *= tem;
      ep[1] *= tem;
      ep[2] *= tem;

      utemp[0] = p.fd(part_idx,0) + ep[0];
      utemp[1] = p.fd(part_idx,1) + ep[1];
      utemp[2] = p.fd(part_idx,2) + ep[2];

      // Get time centered gamma
      u2 = utemp[0]*utemp[0] + utemp[1]*utemp[1] + utemp[2]*utemp[2];

      gamma = sqrt(u2+1);

      // accumulate time centered energy
      // this is done in double precicion always
      if (report_energy) loc_ene += q.fd(part_idx) * u2 / (gamma + _1_0_p_double);

      gam_tem = tem / gamma;

      bp[0] *= gam_tem;
      bp[1] *= gam_tem;
      bp[2] *= gam_tem;

      p.fd(part_idx,0) = utemp[0] + utemp[1] * bp[2];
      p.fd(part_idx,1) = utemp[1] + utemp[2] * bp[0];
      p.fd(part_idx,2) = utemp[2] + utemp[0] * bp[1];

      p.fd(part_idx,0) -= utemp[2] * bp[1];
      p.fd(part_idx,1) -= utemp[0] * bp[2];
      p.fd(part_idx,2) -= utemp[1] * bp[0];

      otsq = _2_0_p_k_part / ( ( (_1_0_p_k_part + bp[0]*bp[0]) + bp[1]*bp[1] ) + bp[2]*bp[2] );

      bp[0] *= otsq;
      bp[1] *= otsq;
      bp[2] *= otsq;

      utemp[0] += p.fd(part_idx,1) * bp[2];
      utemp[1] += p.fd(part_idx,2) * bp[0];
      utemp[2] += p.fd(part_idx,0) * bp[1];

      utemp[0] -= p.fd(part_idx,2) * bp[1];
      utemp[1] -= p.fd(part_idx,0) * bp[2];
      utemp[2] -= p.fd(part_idx,1) * bp[0];

      // Perform second half of electric field acceleration.
      p.fd(part_idx,0) = utemp[0] + ep[0];
      p.fd(part_idx,1) = utemp[1] + ep[1];
      p.fd(part_idx,2) = utemp[2] + ep[2];

      // Note: this only works if chunk size divisible by block size (which we enforce at initialization)
      // See sort for how to do this more generally, if we decide to do that eventually
      part_idx += BLOCK_SIZE;

    } // end loop over chunk of particles ``while ( part_idx < num_par )``

    num_par_proc_tot -= num_par_proc_chunk;
    chunk++;

  } // end loop over chunks

  if (report_energy) {
    accum_energy[tx] = loc_ene;
    if (tx==0) accum_energy_last = accum_energy[BLOCK_SIZE-1];
    prescan<BLOCK_SIZE>( accum_energy );
    __syncthreads();
    if (tx==0) {
      atomicAdd( energy+til, accum_energy[BLOCK_SIZE-1] + accum_energy_last );
    }
  }

}


// ---------------------------------------------------------------------------------------
// Advance velocities and calculate time-centered total energy
// with fields centered at cell corner
// ---------------------------------------------------------------------------------------
template <int BLOCK_SIZE>
__global__ void DUDT_BORIS_1D_CS( t_chunk ***chunks, int *block_start_tiles,
                                  int *block_start_chunks, int *block_num_par,
                                  p_k_part rqm, t_f1<p_k_fld> *b_f1,
                                  t_f1<p_k_fld> *e_f1, p_double *energy,
                                  bool report_energy, p_double dt, int chunk_size ){

  const int bx = blockIdx.x;
  const int tx = threadIdx.x;

  // Shared variables that are the same between all threads
  __shared__ int til;
  __shared__ p_k_part tem;
  __shared__ t_f1<p_k_fld> b_f1_til, e_f1_til;

  // Shared variables that will be prefix scanned
  __shared__ p_double accum_energy[BLOCK_SIZE];

  if (tx==0) {
    til = block_start_tiles[bx];

    b_f1_til = b_f1[til];
    e_f1_til = e_f1[til];

    // get factor that includes timestep and charge-to-mass ratio
    // note that the charge-to-mass ratio is used since the
    // momentum p is not the total momentum of a particles but
    // the momentum per unit restmass (which is the electron mass)
    tem = p_k_part( _0_5_p_double * dt / rqm );
  }
  __syncthreads();

  // const t_f1<p_k_fld> b_f1_til = b_f1[til];
  // const t_f1<p_k_fld> e_f1_til = e_f1[til];

  int num_par_proc_tot = block_num_par[bx];
  t_chunk **chunk = chunks[til] + block_start_chunks[bx];

  int i1, num_par_proc_chunk, part_idx;

  p_k_part u2, gamma, gam_tem, otsq;
  p_k_part bp[p_p_dim], ep[p_p_dim], utemp[p_p_dim];

  p_k_fld dx1;
  p_k_fld f1, f2, f3;

  p_k_fld w1_buff[UP-LP+1];

  p_k_fld *const w1 = w1_buff - LP;

  t_f1<p_k_part> x = t_f1<p_k_part>( chunk_size, p_x_dim );
  t_f1<p_k_part> p = t_f1<p_k_part>( chunk_size, p_p_dim );
  t_f0<p_k_part> q = t_f0<p_k_part>( chunk_size );
  t_f1<int>     ix = t_f1<int>( chunk_size, p_x_dim );

  const size_t chunk_size_t = sizeof( p_k_part ) * chunk_size;
  const size_t x_offset = 0;
  const size_t p_offset = x_offset + p_x_dim * chunk_size_t;
  const size_t q_offset = p_offset +  p_p_dim * chunk_size_t;
  const size_t ix_offset = q_offset + chunk_size_t;

  p_double accum_energy_last, loc_ene = 0;

  while( num_par_proc_tot > 0 ) {

    x.g  = (p_k_part *) ( *chunk + x_offset );
    p.g  = (p_k_part *) ( *chunk + p_offset );
    q.g  = (p_k_part *) ( *chunk + q_offset );
    ix.g = (int *)      ( *chunk + ix_offset);

    // modify bp & ep to include timestep and charge-to-mass ratio
    // and perform half the electric field acceleration.
    // Result is stored in UTEMP.

    // loop over chunks of particles until all particles in tile have been processed
    num_par_proc_chunk = min(num_par_proc_tot, chunk_size);
    part_idx = tx;

    while ( part_idx < num_par_proc_chunk ){

      // --------------
      // dudt_boris()
      // --------------
      i1 = ix.fd(part_idx,0) - 1;
      dx1 = p_k_fld( x.fd(part_idx,0) );

      // get spline weitghts for x
      SPLINE( dx1, w1 );

      // Interpolate E - Field
      f1 = 0;
      f2 = 0;
      f3 = 0;

      for (int k1=LP; k1<=UP; ++k1){
        f1 += e_f1_til.fd(0,i1+k1) * w1[k1];
        f2 += e_f1_til.fd(1,i1+k1) * w1[k1];
        f3 += e_f1_til.fd(2,i1+k1) * w1[k1];
      }

      ep[0] = f1;
      ep[1] = f2;
      ep[2] = f3;

      // Interpolate B - Field
      f1 = 0;
      f2 = 0;
      f3 = 0;

      for (int k1=LP; k1<=UP; ++k1){
        f1 += b_f1_til.fd(0,i1+k1) * w1[k1];
        f2 += b_f1_til.fd(1,i1+k1) * w1[k1];
        f3 += b_f1_til.fd(2,i1+k1) * w1[k1];
      }

      bp[0] = f1;
      bp[1] = f2;
      bp[2] = f3;

      // --------------

      ep[0] *= tem;
      ep[1] *= tem;
      ep[2] *= tem;

      utemp[0] = p.fd(part_idx,0) + ep[0];
      utemp[1] = p.fd(part_idx,1) + ep[1];
      utemp[2] = p.fd(part_idx,2) + ep[2];

      // Get time centered gamma
      u2 = utemp[0]*utemp[0] + utemp[1]*utemp[1] + utemp[2]*utemp[2];

      gamma = sqrt(u2+1);

      // accumulate time centered energy
      // this is done in double precicion always
      if (report_energy) loc_ene += q.fd(part_idx) * u2 / (gamma + _1_0_p_double);

      gam_tem = tem / gamma;

      bp[0] *= gam_tem;
      bp[1] *= gam_tem;
      bp[2] *= gam_tem;

      p.fd(part_idx,0) = utemp[0] + utemp[1] * bp[2];
      p.fd(part_idx,1) = utemp[1] + utemp[2] * bp[0];
      p.fd(part_idx,2) = utemp[2] + utemp[0] * bp[1];

      p.fd(part_idx,0) -= utemp[2] * bp[1];
      p.fd(part_idx,1) -= utemp[0] * bp[2];
      p.fd(part_idx,2) -= utemp[1] * bp[0];

      otsq = _2_0_p_k_part / ( ( (_1_0_p_k_part + bp[0]*bp[0]) + bp[1]*bp[1] ) + bp[2]*bp[2] );

      bp[0] *= otsq;
      bp[1] *= otsq;
      bp[2] *= otsq;

      utemp[0] += p.fd(part_idx,1) * bp[2];
      utemp[1] += p.fd(part_idx,2) * bp[0];
      utemp[2] += p.fd(part_idx,0) * bp[1];

      utemp[0] -= p.fd(part_idx,2) * bp[1];
      utemp[1] -= p.fd(part_idx,0) * bp[2];
      utemp[2] -= p.fd(part_idx,1) * bp[0];

      // Perform second half of electric field acceleration.
      p.fd(part_idx,0) = utemp[0] + ep[0];
      p.fd(part_idx,1) = utemp[1] + ep[1];
      p.fd(part_idx,2) = utemp[2] + ep[2];

      // Note: this only works if chunk size divisible by block size (which we enforce at initialization)
      // See sort for how to do this more generally, if we decide to do that eventually
      part_idx += BLOCK_SIZE;

    } // end loop over chunk of particles ``while ( part_idx < num_par )``

    num_par_proc_tot -= num_par_proc_chunk;
    chunk++;

  } // end loop over chunks

  if (report_energy) {
    accum_energy[tx] = loc_ene;
    if (tx==0) accum_energy_last = accum_energy[BLOCK_SIZE-1];
    prescan<BLOCK_SIZE>( accum_energy );
    __syncthreads();
    if (tx==0) {
      atomicAdd( energy+til, accum_energy[BLOCK_SIZE-1] + accum_energy_last );
    }
  }

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
template <int BLOCK_SIZE>
__global__ void DUDT_BORIS_2D_CS( t_chunk ***chunks, int *block_start_tiles,
                                  int *block_start_chunks, int *block_num_par,
                                  p_k_part rqm, t_f2<p_k_fld> *b_f2,
                                  t_f2<p_k_fld> *e_f2, p_double *energy,
                                  bool report_energy, p_double dt, int chunk_size ){

  const int bx = blockIdx.x;
  const int tx = threadIdx.x;

  // Shared variables that are the same between all threads
  __shared__ int til;
  __shared__ p_k_part tem;
  __shared__ t_f2<p_k_fld> b_f2_til, e_f2_til;

  // Shared variables that will be prefix scanned
  __shared__ p_double accum_energy[BLOCK_SIZE];

  if (tx==0) {
    til = block_start_tiles[bx];

    b_f2_til = b_f2[til];
    e_f2_til = e_f2[til];

    // get factor that includes timestep and charge-to-mass ratio
    // note that the charge-to-mass ratio is used since the
    // momentum p is not the total momentum of a particles but
    // the momentum per unit restmass (which is the electron mass)
    tem = p_k_part( _0_5_p_double * dt / rqm );
  }
  __syncthreads();

  // const t_f2<p_k_fld> b_f2_til = b_f2[til];
  // const t_f2<p_k_fld> e_f2_til = e_f2[til];

  int num_par_proc_tot = block_num_par[bx];
  t_chunk **chunk = chunks[til] + block_start_chunks[bx];

  int i1, i2, num_par_proc_chunk, part_idx;

  p_k_part u2, gamma, gam_tem, otsq;
  p_k_part bp[p_p_dim], ep[p_p_dim], utemp[p_p_dim];

  p_k_fld dx1, dx2;
  p_k_fld f1, f2, f3, f1line, f2line, f3line;

  p_k_fld w1_buff[UP-LP+1], w2_buff[UP-LP+1];

  p_k_fld *const w1 = w1_buff - LP;
  p_k_fld *const w2 = w2_buff - LP;

  t_f1<p_k_part> x = t_f1<p_k_part>( chunk_size, p_x_dim );
  t_f1<p_k_part> p = t_f1<p_k_part>( chunk_size, p_p_dim );
  t_f0<p_k_part> q = t_f0<p_k_part>( chunk_size );
  t_f1<int>     ix = t_f1<int>( chunk_size, p_x_dim );

  const size_t chunk_size_t = sizeof( p_k_part ) * chunk_size;
  const size_t x_offset = 0;
  const size_t p_offset = x_offset + p_x_dim * chunk_size_t;
  const size_t q_offset = p_offset +  p_p_dim * chunk_size_t;
  const size_t ix_offset = q_offset + chunk_size_t;

  p_double accum_energy_last, loc_ene = 0;

  while( num_par_proc_tot > 0 ) {

    x.g  = (p_k_part *) ( *chunk + x_offset );
    p.g  = (p_k_part *) ( *chunk + p_offset );
    q.g  = (p_k_part *) ( *chunk + q_offset );
    ix.g = (int *)      ( *chunk + ix_offset);

    // modify bp & ep to include timestep and charge-to-mass ratio
    // and perform half the electric field acceleration.
    // Result is stored in UTEMP.

    // loop over chunks of particles until all particles in tile have been processed
    num_par_proc_chunk = min(num_par_proc_tot, chunk_size);
    part_idx = tx;

    while ( part_idx < num_par_proc_chunk ){

      // --------------
      // dudt_boris()
      // --------------
      i1 = ix.fd(part_idx,0) - 1;
      i2 = ix.fd(part_idx,1) - 1;

      dx1 = p_k_fld( x.fd(part_idx,0) );
      dx2 = p_k_fld( x.fd(part_idx,1) );

      // get spline weitghts for x and y
      SPLINE( dx1, w1 );
      SPLINE( dx2, w2 );

      // Interpolate E - Field
      f1 = 0;
      f2 = 0;
      f3 = 0;

      for (int k2=LP; k2<=UP; ++k2){
        f1line = 0;
        f2line = 0;
        f3line = 0;

        for (int k1=LP; k1<=UP; ++k1){
          f1line += e_f2_til.fd(0,i1+k1,i2+k2) * w1[k1];
          f2line += e_f2_til.fd(1,i1+k1,i2+k2) * w1[k1];
          f3line += e_f2_til.fd(2,i1+k1,i2+k2) * w1[k1];
        }

        f1 += f1line * w2[k2];
        f2 += f2line * w2[k2];
        f3 += f3line * w2[k2];
      }

      ep[0] = f1;
      ep[1] = f2;
      ep[2] = f3;

      // Interpolate B - Field
      f1 = 0;
      f2 = 0;
      f3 = 0;

      for (int k2=LP; k2<=UP; ++k2){
        f1line = 0;
        f2line = 0;
        f3line = 0;

        for (int k1=LP; k1<=UP; ++k1){
          f1line += b_f2_til.fd(0,i1+k1,i2+k2) * w1[k1];
          f2line += b_f2_til.fd(1,i1+k1,i2+k2) * w1[k1];
          f3line += b_f2_til.fd(2,i1+k1,i2+k2) * w1[k1];
        }

        f1 += f1line * w2[k2];
        f2 += f2line * w2[k2];
        f3 += f3line * w2[k2];
      }

      bp[0] = f1;
      bp[1] = f2;
      bp[2] = f3;

      // --------------

      ep[0] *= tem;
      ep[1] *= tem;
      ep[2] *= tem;

      utemp[0] = p.fd(part_idx,0) + ep[0];
      utemp[1] = p.fd(part_idx,1) + ep[1];
      utemp[2] = p.fd(part_idx,2) + ep[2];

      // Get time centered gamma
      u2 = utemp[0]*utemp[0] + utemp[1]*utemp[1] + utemp[2]*utemp[2];

      gamma = sqrt(u2+1);

      // accumulate time centered energy
      // this is done in double precicion always
      if (report_energy) loc_ene += q.fd(part_idx) * u2 / (gamma + _1_0_p_double);

      gam_tem = tem / gamma;

      bp[0] *= gam_tem;
      bp[1] *= gam_tem;
      bp[2] *= gam_tem;

      p.fd(part_idx,0) = utemp[0] + utemp[1] * bp[2];
      p.fd(part_idx,1) = utemp[1] + utemp[2] * bp[0];
      p.fd(part_idx,2) = utemp[2] + utemp[0] * bp[1];

      p.fd(part_idx,0) -= utemp[2] * bp[1];
      p.fd(part_idx,1) -= utemp[0] * bp[2];
      p.fd(part_idx,2) -= utemp[1] * bp[0];

      otsq = _2_0_p_k_part / ( ( (_1_0_p_k_part + bp[0]*bp[0]) + bp[1]*bp[1] ) + bp[2]*bp[2] );

      bp[0] *= otsq;
      bp[1] *= otsq;
      bp[2] *= otsq;

      utemp[0] += p.fd(part_idx,1) * bp[2];
      utemp[1] += p.fd(part_idx,2) * bp[0];
      utemp[2] += p.fd(part_idx,0) * bp[1];

      utemp[0] -= p.fd(part_idx,2) * bp[1];
      utemp[1] -= p.fd(part_idx,0) * bp[2];
      utemp[2] -= p.fd(part_idx,1) * bp[0];

      // Perform second half of electric field acceleration.
      p.fd(part_idx,0) = utemp[0] + ep[0];
      p.fd(part_idx,1) = utemp[1] + ep[1];
      p.fd(part_idx,2) = utemp[2] + ep[2];

      // Note: this only works if chunk size divisible by block size (which we enforce at initialization)
      // See sort for how to do this more generally, if we decide to do that eventually
      part_idx += BLOCK_SIZE;

    } // end loop over chunk of particles ``while ( part_idx < num_par )``

    num_par_proc_tot -= num_par_proc_chunk;
    chunk++;

  } // end loop over chunks

  if (report_energy) {
    accum_energy[tx] = loc_ene;
    if (tx==0) accum_energy_last = accum_energy[BLOCK_SIZE-1];
    prescan<BLOCK_SIZE>( accum_energy );
    __syncthreads();
    if (tx==0) {
      atomicAdd( energy+til, accum_energy[BLOCK_SIZE-1] + accum_energy_last );
    }
  }

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
template <int BLOCK_SIZE>
__global__ void DUDT_BORIS_3D_CS( t_chunk ***chunks, int *block_start_tiles,
                                  int *block_start_chunks, int *block_num_par,
                                  p_k_part rqm, t_f3<p_k_fld> *b_f3,
                                  t_f3<p_k_fld> *e_f3, p_double *energy,
                                  bool report_energy, p_double dt, int chunk_size ){

  const int bx = blockIdx.x;
  const int tx = threadIdx.x;

  // Shared variables that are the same between all threads
  __shared__ int til;
  __shared__ p_k_part tem;
  __shared__ t_f3<p_k_fld> b_f3_til, e_f3_til;

  // Shared variables that will be prefix scanned
  __shared__ p_double accum_energy[BLOCK_SIZE];

  if (tx==0) {
    til = block_start_tiles[bx];

    b_f3_til = b_f3[til];
    e_f3_til = e_f3[til];

    // get factor that includes timestep and charge-to-mass ratio
    // note that the charge-to-mass ratio is used since the
    // momentum p is not the total momentum of a particles but
    // the momentum per unit restmass (which is the electron mass)
    tem = p_k_part( _0_5_p_double * dt / rqm );
  }
  __syncthreads();

  // const t_f3<p_k_fld> b_f3_til = b_f3[til];
  // const t_f3<p_k_fld> e_f3_til = e_f3[til];

  int num_par_proc_tot = block_num_par[bx];
  t_chunk **chunk = chunks[til] + block_start_chunks[bx];

  int i1, i2, i3, num_par_proc_chunk, part_idx;

  p_k_part u2, gamma, gam_tem, otsq;
  p_k_part bp[p_p_dim], ep[p_p_dim], utemp[p_p_dim];

  p_k_fld dx1, dx2, dx3;
  p_k_fld f1, f2, f3, f1line, f2line, f3line, f1plane, f2plane, f3plane;

  p_k_fld w1_buff[UP-LP+1], w2_buff[UP-LP+1], w3_buff[UP-LP+1];

  p_k_fld *const w1 = w1_buff - LP;
  p_k_fld *const w2 = w2_buff - LP;
  p_k_fld *const w3 = w3_buff - LP;

  t_f1<p_k_part> x = t_f1<p_k_part>( chunk_size, p_x_dim );
  t_f1<p_k_part> p = t_f1<p_k_part>( chunk_size, p_p_dim );
  t_f0<p_k_part> q = t_f0<p_k_part>( chunk_size );
  t_f1<int>     ix = t_f1<int>( chunk_size, p_x_dim );

  const size_t chunk_size_t = sizeof( p_k_part ) * chunk_size;
  const size_t x_offset = 0;
  const size_t p_offset = x_offset + p_x_dim * chunk_size_t;
  const size_t q_offset = p_offset +  p_p_dim * chunk_size_t;
  const size_t ix_offset = q_offset + chunk_size_t;

  p_double accum_energy_last, loc_ene = 0;

  while( num_par_proc_tot > 0 ) {

    x.g  = (p_k_part *) ( *chunk + x_offset );
    p.g  = (p_k_part *) ( *chunk + p_offset );
    q.g  = (p_k_part *) ( *chunk + q_offset );
    ix.g = (int *)      ( *chunk + ix_offset);

    // modify bp & ep to include timestep and charge-to-mass ratio
    // and perform half the electric field acceleration.
    // Result is stored in UTEMP.

    // loop over chunks of particles until all particles in tile have been processed
    num_par_proc_chunk = min(num_par_proc_tot, chunk_size);
    part_idx = tx;

    while ( part_idx < num_par_proc_chunk ){

      // --------------
      // dudt_boris()
      // --------------
      i1 = ix.fd(part_idx,0) - 1;
      i2 = ix.fd(part_idx,1) - 1;
      i3 = ix.fd(part_idx,2) - 1;

      dx1 = p_k_fld( x.fd(part_idx,0) );
      dx2 = p_k_fld( x.fd(part_idx,1) );
      dx3 = p_k_fld( x.fd(part_idx,2) );

      // get spline weitghts for x, y and z
      SPLINE( dx1, w1 );
      SPLINE( dx2, w2 );
      SPLINE( dx3, w3 );

      // Interpolate E - Field
      f1 = 0;
      f2 = 0;
      f3 = 0;

      for (int k3=LP; k3<=UP; ++k3){
        f1plane = 0;
        f2plane = 0;
        f3plane = 0;

        for (int k2=LP; k2<=UP; ++k2){
          f1line = 0;
          f2line = 0;
          f3line = 0;

          for (int k1=LP; k1<=UP; ++k1){
            f1line += e_f3_til.fd(0,i1+k1,i2+k2,i3+k3) * w1[k1];
            f2line += e_f3_til.fd(1,i1+k1,i2+k2,i3+k3) * w1[k1];
            f3line += e_f3_til.fd(2,i1+k1,i2+k2,i3+k3) * w1[k1];
          }

          f1plane += f1line * w2[k2];
          f2plane += f2line * w2[k2];
          f3plane += f3line * w2[k2];
        }

        f1 += f1plane * w3[k3];
        f2 += f2plane * w3[k3];
        f3 += f3plane * w3[k3];

      }

      ep[0] = f1;
      ep[1] = f2;
      ep[2] = f3;

      // Interpolate B - Field
      f1 = 0;
      f2 = 0;
      f3 = 0;

      for (int k3=LP; k3<=UP; ++k3){
        f1plane = 0;
        f2plane = 0;
        f3plane = 0;

        for (int k2=LP; k2<=UP; ++k2){
          f1line = 0;
          f2line = 0;
          f3line = 0;

          for (int k1=LP; k1<=UP; ++k1){
            f1line += b_f3_til.fd(0,i1+k1,i2+k2,i3+k3) * w1[k1];
            f2line += b_f3_til.fd(1,i1+k1,i2+k2,i3+k3) * w1[k1];
            f3line += b_f3_til.fd(2,i1+k1,i2+k2,i3+k3) * w1[k1];
          }

          f1plane += f1line * w2[k2];
          f2plane += f2line * w2[k2];
          f3plane += f3line * w2[k2];
        }

        f1 += f1plane * w3[k3];
        f2 += f2plane * w3[k3];
        f3 += f3plane * w3[k3];

      }

      bp[0] = f1;
      bp[1] = f2;
      bp[2] = f3;

      // --------------

      ep[0] *= tem;
      ep[1] *= tem;
      ep[2] *= tem;

      utemp[0] = p.fd(part_idx,0) + ep[0];
      utemp[1] = p.fd(part_idx,1) + ep[1];
      utemp[2] = p.fd(part_idx,2) + ep[2];

      // Get time centered gamma
      u2 = utemp[0]*utemp[0] + utemp[1]*utemp[1] + utemp[2]*utemp[2];

      gamma = sqrt(u2+1);

      // accumulate time centered energy
      // this is done in double precicion always
      if (report_energy) loc_ene += q.fd(part_idx) * u2 / (gamma + _1_0_p_double);

      gam_tem = tem / gamma;

      bp[0] *= gam_tem;
      bp[1] *= gam_tem;
      bp[2] *= gam_tem;

      p.fd(part_idx,0) = utemp[0] + utemp[1] * bp[2];
      p.fd(part_idx,1) = utemp[1] + utemp[2] * bp[0];
      p.fd(part_idx,2) = utemp[2] + utemp[0] * bp[1];

      p.fd(part_idx,0) -= utemp[2] * bp[1];
      p.fd(part_idx,1) -= utemp[0] * bp[2];
      p.fd(part_idx,2) -= utemp[1] * bp[0];

      otsq = _2_0_p_k_part / ( ( (_1_0_p_k_part + bp[0]*bp[0]) + bp[1]*bp[1] ) + bp[2]*bp[2] );

      bp[0] *= otsq;
      bp[1] *= otsq;
      bp[2] *= otsq;

      utemp[0] += p.fd(part_idx,1) * bp[2];
      utemp[1] += p.fd(part_idx,2) * bp[0];
      utemp[2] += p.fd(part_idx,0) * bp[1];

      utemp[0] -= p.fd(part_idx,2) * bp[1];
      utemp[1] -= p.fd(part_idx,0) * bp[2];
      utemp[2] -= p.fd(part_idx,1) * bp[0];

      // Perform second half of electric field acceleration.
      p.fd(part_idx,0) = utemp[0] + ep[0];
      p.fd(part_idx,1) = utemp[1] + ep[1];
      p.fd(part_idx,2) = utemp[2] + ep[2];

      // Note: this only works if chunk size divisible by block size (which we enforce at initialization)
      // See sort for how to do this more generally, if we decide to do that eventually
      part_idx += BLOCK_SIZE;

    } // end loop over chunk of particles ``while ( part_idx < num_par )``

    num_par_proc_tot -= num_par_proc_chunk;
    chunk++;

  } // end loop over chunks

  if (report_energy) {
    accum_energy[tx] = loc_ene;
    if (tx==0) accum_energy_last = accum_energy[BLOCK_SIZE-1];
    prescan<BLOCK_SIZE>( accum_energy );
    __syncthreads();
    if (tx==0) {
      atomicAdd( energy+til, accum_energy[BLOCK_SIZE-1] + accum_energy_last );
    }
  }

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
template <int BLOCK_SIZE>
__global__ void ADVANCE_DEPOSIT_1D( t_chunk ***chunks, int *block_start_tiles,
                                    int *block_start_chunks, int *block_num_par,
                                    p_double dx1,
                                    t_f1<p_k_fld> *jay_f1, p_double dt,
                                    int chunk_size, int *n_hole, int *n_comm_ext,
                                    int *n_recv_loc ){

  // Block, thread index
  const int bx = blockIdx.x;
  const int tx = threadIdx.x;

  // Shared variables that are the same between all threads
  __shared__ int til;
  __shared__ t_f1<p_k_fld> jay_f1_til;
  __shared__ p_k_part dt_dx1;

  if (tx==0) {
    til = block_start_tiles[bx];

    jay_f1_til = jay_f1[til];

    dt_dx1 = p_k_part( dt/dx1 );
  }
  __syncthreads();

  int num_par_proc_tot = block_num_par[bx];
  t_chunk **chunk = chunks[til] + block_start_chunks[bx];

  p_k_part rgamma;
  p_k_part xbuf;

  int dxi, part_idx, num_par_proc_chunk;

  // data for cpu type :: t_vp1D, we just store in local per thread vars
  p_k_part x0[2], x1[2], qq[2], vy[2], vz[2];
  int ii[2];

  p_k_part xint, delta, vz_, vy_;

  p_k_fld qnx, qvy, qvz;
  p_k_fld pw1;

  p_k_fld S0x_buff[P1-P0+1], S1x_buff[P1-P0+1];
  p_k_fld wl1_buff[P1-P0+1];

  // Obtain offset pointers to the above buffers to enable non-zero indexing
  p_k_fld *const S0x = S0x_buff - P0;
  p_k_fld *const S1x = S1x_buff - P0;
  p_k_fld *const wl1 = wl1_buff - P0;

  t_f1<p_k_part> x = t_f1<p_k_part>( chunk_size, p_x_dim );
  t_f1<p_k_part> p = t_f1<p_k_part>( chunk_size, p_p_dim );
  t_f0<p_k_part> q = t_f0<p_k_part>( chunk_size );
  t_f1<int>     ix = t_f1<int>( chunk_size, p_x_dim );

  const size_t chunk_size_t = sizeof( p_k_part ) * chunk_size;
  const size_t x_offset = 0;
  const size_t p_offset = x_offset + p_x_dim * chunk_size_t;
  const size_t q_offset = p_offset +  p_p_dim * chunk_size_t;
  const size_t ix_offset = q_offset + chunk_size_t;

  const p_k_fld dx1_dt = p_k_fld( dx1/dt );

  const p_k_fld c1_3 = _1_0_p_k_fld / p_k_fld(3.0);

  // loop over chunks of particles until all particles in tile have been processed
  while( num_par_proc_tot > 0 ) {

    x.g  = (p_k_part *) ( *chunk + x_offset );
    p.g  = (p_k_part *) ( *chunk + p_offset );
    q.g  = (p_k_part *) ( *chunk + q_offset );
    ix.g = (int *)      ( *chunk + ix_offset);

    num_par_proc_chunk = min(num_par_proc_tot, chunk_size);
    part_idx = tx;

    // loop over particles in a chunk, BLOCK_SIZE at a time
    while ( part_idx < num_par_proc_chunk ) {

      rgamma = _1_0_p_k_part /
          sqrt( ( (_1_0_p_k_part + p.fd(part_idx,0)*p.fd(part_idx,0))
                                 + p.fd(part_idx,1)*p.fd(part_idx,1) )
                                 + p.fd(part_idx,2)*p.fd(part_idx,2) );

      xbuf = x.fd(part_idx,0) + ( p.fd(part_idx,0) * rgamma ) * dt_dx1;

      dxi = ntrim(xbuf);

      // --------------
      // dep_current_2d_s*()
      // --------------

      // --------------
      // split_2d()
      // --------------

      // Splits particle trajectories so that all virtual particles have a motion starting and ending
      // in the same cell. This routines also calculates vz for each virtual particle.

      // create virtual particles that correspond to a motion that starts and ends
      // in the same grid cell

      int k = 0;

      if ( dxi==0 ) {

        x0[k] = x.fd(part_idx,0);
        x1[k] = xbuf;
        qq[k]  = q.fd(part_idx);

        vy[k] = p.fd(part_idx,1)*rgamma;
        vz[k] = p.fd(part_idx,2)*rgamma;

        ii[k] = ix.fd(part_idx,0) - 1;

        k=k+1;

      } else {

        xint  = _0_5_p_k_fld * dxi;
        delta = ( xint - x.fd(part_idx,0) ) / ( xbuf - x.fd(part_idx,0) );

        vy_ = p.fd(part_idx,1)*rgamma;
        vz_ = p.fd(part_idx,2)*rgamma;

        x0[k] = x.fd(part_idx,0);
        x1[k] = xint;

        qq[k] = q.fd(part_idx);

        vy[k] = vy_ * delta;
        vz[k] = vz_ * delta;

        ii[k] = ix.fd(part_idx,0) - 1;

        k=k+1;

        x0[k] = -xint;

        x1[k] = xbuf - dxi;

        qq[k] = q.fd(part_idx);

        vy[k] = vy_ * (1-delta);
        vz[k] = vz_ * (1-delta);

        ii[k] = ix.fd(part_idx,0) + dxi - 1;

        k=k+1;

      }
      __syncwarp();

      int nsplit = k;

      // --------------

      // now accumulate jay iooping through all virtual particles
      // and shifting grid indexes

      for (int i=0; i<nsplit; ++i){

        // Normalize charge
        qnx = qq[i] * dx1_dt;
        qvy = qq[i] * vy[i];
        qvz = qq[i] * vz[i];

        // get spline weights for x and y
        SPLINE( x0[i], S0x );
        SPLINE( x1[i], S1x );

        // get longitudinal motion weights
        // the last vaiue is set to 0 so we can accumulate all current components in a single
        // pass
        WL( qnx, x0[i], x1[i], wl1 );
        wl1[P1] = 0;

        // accumulate j1, j2, j3 in a single pass
        for (int k1=P0; k1<=P1; ++k1){

          pw1 = _0_5_p_k_fld * (S0x[k1]+S1x[k1]);

          atomicAdd( jay_f1_til.fdp(0,ii[i]+k1), wl1[k1] );
          atomicAdd( jay_f1_til.fdp(1,ii[i]+k1), qvy * pw1 );
          atomicAdd( jay_f1_til.fdp(2,ii[i]+k1), qvz * pw1 );

        }

      }

      // --------------

      x.fd(part_idx,0) = xbuf - dxi;
      ix.fd(part_idx,0) += dxi;

      // Note: this only works if chunk size divisible by block size (which we enforce at initialization)
      // See sort for how to do this more generally, if we decide to do that eventually
      part_idx += BLOCK_SIZE;

    } // end loop over chunk of particles ``while ( part_idx < num_par_proc_chunk )``

    num_par_proc_tot -= num_par_proc_chunk;
    chunk++;

  } // end loop over chunks  ``while ( num_par_remainin > 0 )``

  // one thread on each tile zero out counters to be used in sort
  if ( tx==0 and block_start_chunks[bx]==0 ) {
    n_hole[til] = 0;
    n_comm_ext[til] = 0;
    n_recv_loc[til] = 0;
  }

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
template <int BLOCK_SIZE>
__global__ void ADVANCE_DEPOSIT_2D( t_chunk ***chunks, int *block_start_tiles,
                                    int *block_start_chunks, int *block_num_par,
                                    p_double dx1, p_double dx2,
                                    t_f2<p_k_fld> *jay_f2, p_double dt,
                                    int chunk_size, int *n_hole, int *n_comm_ext,
                                    int *n_recv_loc ){

  // Block, thread index
  const int bx = blockIdx.x;
  const int tx = threadIdx.x;

  // Shared variables that are the same between all threads
  __shared__ int til;
  __shared__ t_f2<p_k_fld> jay_f2_til;
  __shared__ p_k_part dt_dx1, dt_dx2;

  if (tx==0) {
    til = block_start_tiles[bx];

    jay_f2_til = jay_f2[til];

    dt_dx1 = p_k_part( dt/dx1 );
    dt_dx2 = p_k_part( dt/dx2 );
  }
  __syncthreads();

  int num_par_proc_tot = block_num_par[bx];
  t_chunk **chunk = chunks[til] + block_start_chunks[bx];

  p_k_part rgamma;
  p_k_part xbuf[2];

  int dxi[2], part_idx, num_par_proc_chunk;

  // data for cpu type :: t_vp2D, we just store in local per thread vars
  p_k_part x0[3], y0[3], x1[3], y1[3], qq[3], vz[3];
  int ii[3], jj[3];

  p_k_part xint, yint, xint2, yint2, delta, vz_, vzint;

  p_k_fld qnx, qny, qvz;
  p_k_fld tmp1, tmp2;

  p_k_fld S0x_buff[P1-P0+1], S1x_buff[P1-P0+1], S0y_buff[P1-P0+1], S1y_buff[P1-P0+1];
  p_k_fld wp1_buff[P1-P0+1], wp2_buff[P1-P0+1];
  p_k_fld wl1_buff[P1-P0+1], wl2_buff[P1-P0+1];

  // Obtain offset pointers to the above buffers to enable non-zero indexing
  p_k_fld *const S0x = S0x_buff - P0;
  p_k_fld *const S1x = S1x_buff - P0;
  p_k_fld *const S0y = S0y_buff - P0;
  p_k_fld *const S1y = S1y_buff - P0;
  p_k_fld *const wp1 = wp1_buff - P0;
  p_k_fld *const wp2 = wp2_buff - P0;
  p_k_fld *const wl1 = wl1_buff - P0;
  p_k_fld *const wl2 = wl2_buff - P0;

  t_f1<p_k_part> x = t_f1<p_k_part>( chunk_size, p_x_dim );
  t_f1<p_k_part> p = t_f1<p_k_part>( chunk_size, p_p_dim );
  t_f0<p_k_part> q = t_f0<p_k_part>( chunk_size );
  t_f1<int>     ix = t_f1<int>( chunk_size, p_x_dim );

  const size_t chunk_size_t = sizeof( p_k_part ) * chunk_size;
  const size_t x_offset = 0;
  const size_t p_offset = x_offset + p_x_dim * chunk_size_t;
  const size_t q_offset = p_offset +  p_p_dim * chunk_size_t;
  const size_t ix_offset = q_offset + chunk_size_t;

  const p_k_fld c1_3 = _1_0_p_k_fld / p_k_fld(3.0);
  const p_k_fld jnorm1 = p_k_fld( dx1/dt/2 );
  const p_k_fld jnorm2 = p_k_fld( dx2/dt/2 );

  // loop over chunks of particles until all particles in tile have been processed
  while( num_par_proc_tot > 0 ) {

    x.g  = (p_k_part *) ( *chunk + x_offset );
    p.g  = (p_k_part *) ( *chunk + p_offset );
    q.g  = (p_k_part *) ( *chunk + q_offset );
    ix.g = (int *)      ( *chunk + ix_offset);

    num_par_proc_chunk = min(num_par_proc_tot, chunk_size);
    part_idx = tx;

    // loop over particles in a chunk, BLOCK_SIZE at a time
    while ( part_idx < num_par_proc_chunk ) {

      rgamma = _1_0_p_k_part /
          sqrt( ( (_1_0_p_k_part + p.fd(part_idx,0)*p.fd(part_idx,0))
                                 + p.fd(part_idx,1)*p.fd(part_idx,1) )
                                 + p.fd(part_idx,2)*p.fd(part_idx,2) );

      xbuf[0] = x.fd(part_idx,0) + ( p.fd(part_idx,0) * rgamma ) * dt_dx1;
      xbuf[1] = x.fd(part_idx,1) + ( p.fd(part_idx,1) * rgamma ) * dt_dx2;

      dxi[0] = ntrim(xbuf[0]);
      dxi[1] = ntrim(xbuf[1]);

      // --------------
      // dep_current_2d_s*()
      // --------------
      // Linear interpolation
      // Accumulate all components in a single pass

      // --------------
      // split_2d()
      // --------------

      // Splits particle trajectories so that all virtual particles have a motion starting and ending
      // in the same cell. This routines also calculates vz for each virtual particle.

      // create virtual particles that correspond to a motion that starts and ends
      // in the same grid cell

      int k = 0;

      int cross = abs(dxi[0]) + 2 * abs(dxi[1]);

      switch (cross) {
        // no cross
        case 0:
          x0[k] = x.fd(part_idx,0);
          y0[k] = x.fd(part_idx,1);
          x1[k] = xbuf[0];
          y1[k] = xbuf[1];
          qq[k]  = q.fd(part_idx);

          vz[k] = p.fd(part_idx,2)*rgamma;

          ii[k] = ix.fd(part_idx,0) - 1;
          jj[k] = ix.fd(part_idx,1) - 1;

          // printf("--CUDA DEBUG: til %d, part %d k %d x0 %f, y0 %f, x1 %f, y1 %f, qq %f, vz %f, ii %d, jj %d \n",
          //     bx+1, part_idx+1, k, x0[k], y0[k], x1[k], y1[k], qq[k], vz[k], ii[k], jj[k]);

          k=k+1;

          break;

        // x cross only
        case 1:

          xint  = _0_5_p_k_fld * dxi[0];
          delta = ( xint - x.fd(part_idx,0) ) / ( xbuf[0] - x.fd(part_idx,0) );
          yint  =  x.fd(part_idx,1) + (xbuf[1] - x.fd(part_idx,1)) * delta;

          vz_ = p.fd(part_idx,2)*rgamma;

          x0[k] = x.fd(part_idx,0);
          y0[k] = x.fd(part_idx,1);

          x1[k] = xint;
          y1[k] = yint;

          qq[k] = q.fd(part_idx);
          vz[k] = vz_ * delta;

          ii[k] = ix.fd(part_idx,0) - 1;
          jj[k] = ix.fd(part_idx,1) - 1;

          k=k+1;

          x0[k] = -xint;
          y0[k] = yint;

          x1[k] = xbuf[0] - dxi[0];
          y1[k] = xbuf[1];

          qq[k] = q.fd(part_idx);
          vz[k] = vz_ * (1-delta);

          ii[k] = ix.fd(part_idx,0) + dxi[0] - 1;
          jj[k] = ix.fd(part_idx,1) - 1;

          k=k+1;

          break;

        // y cross only
        case 2:

          yint = _0_5_p_k_fld * dxi[1];
          delta = ( yint - x.fd(part_idx,1) ) / ( xbuf[1] - x.fd(part_idx,1));
          xint =  x.fd(part_idx,0) + (xbuf[0] - x.fd(part_idx,0)) * delta;

          vz_ = p.fd(part_idx,2)*rgamma;

          x0[k] = x.fd(part_idx,0);
          y0[k] = x.fd(part_idx,1);

          x1[k] = xint;
          y1[k] = yint;
          qq[k] = q.fd(part_idx);
          vz[k] = vz_ * delta;

          ii[k] = ix.fd(part_idx,0) - 1;
          jj[k] = ix.fd(part_idx,1) - 1;

          k=k+1;

          x0[k] = xint;
          y0[k] = -yint;

          x1[k] = xbuf[0];
          y1[k] = xbuf[1] - dxi[1];
          qq[k] = q.fd(part_idx);
          vz[k] = vz_ * (1-delta);

          ii[k] = ix.fd(part_idx,0) - 1;
          jj[k] = ix.fd(part_idx,1) + dxi[1] - 1;

          k=k+1;

          break;

        // x, y cross
        case 3:

          // split in x direction first
          xint = _0_5_p_k_fld * dxi[0];
          delta = ( xint - x.fd(part_idx,0) ) / ( xbuf[0] - x.fd(part_idx,0));
          yint =  x.fd(part_idx,1) + ( xbuf[1] - x.fd(part_idx,1)) * delta;

          vz_ = p.fd(part_idx,2)*rgamma;

          // check if y intersection occured for 1st or 2nd split
          if ( yint >= -_0_5_p_k_fld && yint < _0_5_p_k_fld ){

            // no y cross on 1st vp
            x0[k] = x.fd(part_idx,0);
            y0[k] = x.fd(part_idx,1);

            x1[k] = xint;
            y1[k] = yint;

            qq[k] = q.fd(part_idx);
            vz[k] = vz_ * delta;

            vzint = vz_*(1-delta);

            ii[k] = ix.fd(part_idx,0) - 1;
            jj[k] = ix.fd(part_idx,1) - 1;

            k=k+1;

            // y split 2nd vp
            yint2 = _0_5_p_k_fld * dxi[1];
            delta = ( yint2 - yint ) / ( xbuf[1] - yint );
            xint2 = -xint + ( xbuf[0] - xint ) * delta;

            x0[k] = -xint;
            y0[k] =  yint;

            x1[k] = xint2;
            y1[k] = yint2;

            qq[k] = q.fd(part_idx);
            vz[k] = vzint * delta;

            ii[k] = ix.fd(part_idx,0) + dxi[0] - 1;
            jj[k] = ix.fd(part_idx,1) - 1;

            k=k+1;

            x0[k] = xint2;
            y0[k] = -yint2;

            x1[k] = xbuf[0] - dxi[0];
            y1[k] = xbuf[1] - dxi[1];
            qq[k] = q.fd(part_idx);
            vz[k] = vzint * (1-delta);

            ii[k] = ix.fd(part_idx,0) + dxi[0] - 1;
            jj[k] = ix.fd(part_idx,1) + dxi[1] - 1;

            k=k+1;

          } else {

            vzint = vz_ * delta;

            // y split 1st vp
            yint2 = _0_5_p_k_fld * dxi[1];
            delta = ( yint2 - x.fd(part_idx,1) ) / ( yint - x.fd(part_idx,1));
            xint2 = x.fd(part_idx,0) + (xint - x.fd(part_idx,0)) * delta;

            x0[k] = x.fd(part_idx,0);
            y0[k] = x.fd(part_idx,1);

            x1[k] = xint2;
            y1[k] = yint2;

            qq[k] = q.fd(part_idx);
            vz[k] = vzint*delta;

            ii[k] = ix.fd(part_idx,0) - 1;
            jj[k] = ix.fd(part_idx,1) - 1;

            k=k+1;

            x0[k] = xint2;
            y0[k] = -yint2;

            x1[k] = xint;
            y1[k] = yint - dxi[1];

            qq[k] = q.fd(part_idx);
            vz[k] = vzint*(1-delta);

            ii[k] = ix.fd(part_idx,0) - 1;
            jj[k] = ix.fd(part_idx,1) + dxi[1] - 1;

            k=k+1;

            // no y cross on second vp
            x0[k] = - xint;
            y0[k] = yint - dxi[1];
            x1[k] = xbuf[0] - dxi[0];
            y1[k] = xbuf[1] - dxi[1];
            qq[k] = q.fd(part_idx);
            vz[k] = vz_ - vzint;

            ii[k] = ix.fd(part_idx,0) + dxi[0] - 1;
            jj[k] = ix.fd(part_idx,1) + dxi[1] - 1;

            k=k+1;

          }

          break;

      }
      __syncwarp();

      int nsplit = k;

      // --------------

      // now accumulate jay iooping through all virtual particles
      // and shifting grid indexes

      for (int i=0; i<nsplit; ++i){

        // Normalize charge
        qnx = qq[i] * jnorm1;
        qny = qq[i] * jnorm2;
        qvz = c1_3 * qq[i] * vz[i];

        // get spline weights for x and y
        SPLINE( x0[i], S0x );
        SPLINE( x1[i], S1x );

        SPLINE( y0[i], S0y );
        SPLINE( y1[i], S1y );

        // get longitudinal motion weights
        // the last vaiue is set to 0 so we can accumulate all current components in a single
        // pass
        WL( qnx, x0[i], x1[i], wl1 );
        wl1[P1] = 0;

        WL( qny, y0[i], y1[i], wl2 );
        wl2[P1] = 0;

        // get perpendicular motion weights
        // (a division by 2 was moved to the jnorm vars. above)
        for (int k1=P0; k1<=P1; ++k1){
          wp1[k1] = S0y[k1] + S1y[k1];
          wp2[k1] = S0x[k1] + S1x[k1];
        }

        // accumulate j1, j2, j3 in a single pass
        for (int k2=P0; k2<=P1; ++k2){
          for (int k1=P0; k1<=P1; ++k1){

            tmp1 = S0x[k1]*S0y[k2] + S1x[k1]*S1y[k2];
            tmp2 = S0x[k1]*S1y[k2] + S1x[k1]*S0y[k2];

            atomicAdd( jay_f2_til.fdp(0,ii[i]+k1,jj[i]+k2), wl1[k1] * wp1[k2] );
            atomicAdd( jay_f2_til.fdp(1,ii[i]+k1,jj[i]+k2), wp2[k1] * wl2[k2] );
            atomicAdd( jay_f2_til.fdp(2,ii[i]+k1,jj[i]+k2), qvz * ( tmp1 + _0_5_p_k_fld * tmp2 ) );

          }
        }

      }

      // --------------

      x.fd(part_idx,0) = xbuf[0] - dxi[0];
      x.fd(part_idx,1) = xbuf[1] - dxi[1];
      ix.fd(part_idx,0) += dxi[0];
      ix.fd(part_idx,1) += dxi[1];

      // Note: this only works if chunk size divisible by block size (which we enforce at initialization)
      // See sort for how to do this more generally, if we decide to do that eventually
      part_idx += BLOCK_SIZE;

    } // end loop over chunk of particles ``while ( part_idx < num_par_proc_chunk )``

    num_par_proc_tot -= num_par_proc_chunk;
    chunk++;

  } // end loop over chunks  ``while ( num_par_remainin > 0 )``

  // one thread on each tile zero out counters to be used in sort
  if ( tx==0 and block_start_chunks[bx]==0 ) {
    n_hole[til] = 0;
    n_comm_ext[til] = 0;
    n_recv_loc[til] = 0;
  }

}


// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
template <int BLOCK_SIZE>
__global__ void ADVANCE_DEPOSIT_3D( t_chunk ***chunks, int *block_start_tiles,
                                    int *block_start_chunks, int *block_num_par,
                                    p_double dx1, p_double dx2, p_double dx3,
                                    t_f3<p_k_fld> *jay_f3, p_double dt,
                                    int chunk_size, int *n_hole, int *n_comm_ext,
                                    int *n_recv_loc ){

#define NB (P1-P0+1)

  // Block, thread index
  const int bx = blockIdx.x;
  const int tx = threadIdx.x;

  // Shared variables that are the same between all threads
  __shared__ int til;
  __shared__ t_f3<p_k_fld> jay_f3_til;
  __shared__ p_k_part dt_dx1, dt_dx2, dt_dx3;

  if (tx==0) {
    til = block_start_tiles[bx];

    jay_f3_til = jay_f3[til];

    dt_dx1 = p_k_part( dt/dx1 );
    dt_dx2 = p_k_part( dt/dx2 );
    dt_dx3 = p_k_part( dt/dx3 );
  }
  __syncthreads();

  int num_par_proc_tot = block_num_par[bx];
  t_chunk **chunk = chunks[til] + block_start_chunks[bx];

  p_k_part rgamma;
  p_k_part xbuf[3];

  int dxi[3], part_idx, num_par_proc_chunk;

  // data for cpu type :: t_vp3D, we just store in local per thread vars
  p_k_part x0[4], y0[4], z0[4], x1[4], y1[4], z1[4], qq[4];
  int ii[4], jj[4], kk[4];

  p_k_part xint, yint, zint, xint2, yint2, zint2, delta;

  p_k_fld qnx, qny, qnz;

  p_k_fld S0x_buff[NB], S1x_buff[NB], S0y_buff[NB], S1y_buff[NB], S0z_buff[NB], S1z_buff[NB];
  p_k_fld wp1_buff[NB*NB], wp2_buff[NB*NB], wp3_buff[NB*NB];
  p_k_fld wl1_buff[NB], wl2_buff[NB], wl3_buff[NB];

  // Obtain offset pointers to the above buffers to enable non-zero indexing
  p_k_fld *const S0x = S0x_buff - P0;
  p_k_fld *const S1x = S1x_buff - P0;
  p_k_fld *const S0y = S0y_buff - P0;
  p_k_fld *const S1y = S1y_buff - P0;
  p_k_fld *const S0z = S0z_buff - P0;
  p_k_fld *const S1z = S1z_buff - P0;
  p_k_fld *const wp1 = wp1_buff - P0 - P0*NB;
  p_k_fld *const wp2 = wp2_buff - P0 - P0*NB;
  p_k_fld *const wp3 = wp3_buff - P0 - P0*NB;
  p_k_fld *const wl1 = wl1_buff - P0;
  p_k_fld *const wl2 = wl2_buff - P0;
  p_k_fld *const wl3 = wl3_buff - P0;

  t_f1<p_k_part> x = t_f1<p_k_part>( chunk_size, p_x_dim );
  t_f1<p_k_part> p = t_f1<p_k_part>( chunk_size, p_p_dim );
  t_f0<p_k_part> q = t_f0<p_k_part>( chunk_size );
  t_f1<int>     ix = t_f1<int>( chunk_size, p_x_dim );

  const size_t chunk_size_t = sizeof( p_k_part ) * chunk_size;
  const size_t x_offset = 0;
  const size_t p_offset = x_offset + p_x_dim * chunk_size_t;
  const size_t q_offset = p_offset +  p_p_dim * chunk_size_t;
  const size_t ix_offset = q_offset + chunk_size_t;

  const p_k_fld jnorm1 = p_k_fld( dx1/dt/3 );
  const p_k_fld jnorm2 = p_k_fld( dx2/dt/3 );
  const p_k_fld jnorm3 = p_k_fld( dx3/dt/3 );

  // loop over chunks of particles until all particles in tile have been processed
  while( num_par_proc_tot > 0 ) {

    x.g  = (p_k_part *) ( *chunk + x_offset );
    p.g  = (p_k_part *) ( *chunk + p_offset );
    q.g  = (p_k_part *) ( *chunk + q_offset );
    ix.g = (int *)      ( *chunk + ix_offset);

    num_par_proc_chunk = min(num_par_proc_tot, chunk_size);
    part_idx = tx;

    // loop over particles in a chunk, BLOCK_SIZE at a time
    while ( part_idx < num_par_proc_chunk ) {

      rgamma = _1_0_p_k_part /
          sqrt( ( (_1_0_p_k_part + p.fd(part_idx,0)*p.fd(part_idx,0) )
                                 + p.fd(part_idx,1)*p.fd(part_idx,1) )
                                 + p.fd(part_idx,2)*p.fd(part_idx,2) );

      xbuf[0] = x.fd(part_idx,0) + ( p.fd(part_idx,0) * rgamma ) * dt_dx1;
      xbuf[1] = x.fd(part_idx,1) + ( p.fd(part_idx,1) * rgamma ) * dt_dx2;
      xbuf[2] = x.fd(part_idx,2) + ( p.fd(part_idx,2) * rgamma ) * dt_dx3;

      dxi[0] = ntrim(xbuf[0]);
      dxi[1] = ntrim(xbuf[1]);
      dxi[2] = ntrim(xbuf[2]);

      // --------------
      // dep_current_2d_s*()
      // --------------
      // Linear interpolation
      // Accumulate all components in a single pass

      // --------------
      // split_2d()
      // --------------

      // Splits particle trajectories so that all virtual particles have a motion starting and ending
      // in the same cell. This routines also calculates vz for each virtual particle.

      // create virtual particles that correspond to a motion that starts and ends
      // in the same grid cell

      int k = 0;

      int cross = abs(dxi[0]) + 2 * abs(dxi[1]) + 4 * abs(dxi[2]);

      switch (cross) {
        // no cross
        case 0:
          x0[k] = x.fd(part_idx,0);
          y0[k] = x.fd(part_idx,1);
          z0[k] = x.fd(part_idx,2);

          x1[k] = xbuf[0];
          y1[k] = xbuf[1];
          z1[k] = xbuf[2];

          qq[k]  = q.fd(part_idx);

          ii[k] = ix.fd(part_idx,0) - 1;
          jj[k] = ix.fd(part_idx,1) - 1;
          kk[k] = ix.fd(part_idx,2) - 1;

          k=k+1;

          break;

        // x cross only
        case 1:

          xint  = _0_5_p_k_fld * dxi[0];
          delta = ( xint - x.fd(part_idx,0) ) / ( xbuf[0] - x.fd(part_idx,0) );
          yint  =  x.fd(part_idx,1) + (xbuf[1] - x.fd(part_idx,1)) * delta;
          zint  =  x.fd(part_idx,2) + (xbuf[2] - x.fd(part_idx,2)) * delta;

          x0[k] = x.fd(part_idx,0);
          y0[k] = x.fd(part_idx,1);
          z0[k] = x.fd(part_idx,2);

          x1[k] = xint;
          y1[k] = yint;
          z1[k] = zint;

          qq[k] = q.fd(part_idx);

          ii[k] = ix.fd(part_idx,0) - 1;
          jj[k] = ix.fd(part_idx,1) - 1;
          kk[k] = ix.fd(part_idx,2) - 1;

          k=k+1;

          x0[k] = -xint;
          y0[k] = yint;
          z0[k] = zint;
          x1[k] = xbuf[0] - dxi[0];
          y1[k] = xbuf[1];
          z1[k] = xbuf[2];

          qq[k] = q.fd(part_idx);

          ii[k] = ix.fd(part_idx,0) + dxi[0] - 1;
          jj[k] = ix.fd(part_idx,1) - 1;
          kk[k] = ix.fd(part_idx,2) - 1;

          k=k+1;

          break;

        // y cross only
        case 2:

          yint = _0_5_p_k_fld * dxi[1];
          delta = ( yint - x.fd(part_idx,1) ) / ( xbuf[1] - x.fd(part_idx,1));
          xint =  x.fd(part_idx,0) + (xbuf[0] - x.fd(part_idx,0)) * delta;
          zint =  x.fd(part_idx,2) + (xbuf[2] - x.fd(part_idx,2)) * delta;

          x0[k] = x.fd(part_idx,0);
          y0[k] = x.fd(part_idx,1);
          z0[k] = x.fd(part_idx,2);

          x1[k] = xint;
          y1[k] = yint;
          z1[k] = zint;

          qq[k] = q.fd(part_idx);

          ii[k] = ix.fd(part_idx,0) - 1;
          jj[k] = ix.fd(part_idx,1) - 1;
          kk[k] = ix.fd(part_idx,2) - 1;

          k=k+1;

          x0[k] = xint;
          y0[k] = -yint;
          z0[k] = zint;
          x1[k] = xbuf[0];
          y1[k] = xbuf[1] - dxi[1];
          z1[k] = xbuf[2];

          qq[k] = q.fd(part_idx);

          ii[k] = ix.fd(part_idx,0) - 1;
          jj[k] = ix.fd(part_idx,1) + dxi[1] - 1;
          kk[k] = ix.fd(part_idx,2) - 1;

          k=k+1;

          break;

        // x, y cross
        case 3:

          xint = _0_5_p_k_fld * dxi[0];
          delta = ( xint - x.fd(part_idx,0) ) / ( xbuf[0] - x.fd(part_idx,0));
          yint =  x.fd(part_idx,1) + ( xbuf[1] - x.fd(part_idx,1)) * delta;
          zint =  x.fd(part_idx,2) + ( xbuf[2] - x.fd(part_idx,2)) * delta;

          if ( yint >= -_0_5_p_k_fld && yint < _0_5_p_k_fld ){

            // yspi = 1

            // no y cross on 1st vp
            x0[k] = x.fd(part_idx,0);
            y0[k] = x.fd(part_idx,1);
            z0[k] = x.fd(part_idx,2);

            x1[k] = xint;
            y1[k] = yint;
            z1[k] = zint;

            qq[k] = q.fd(part_idx);

            ii[k] = ix.fd(part_idx,0) - 1;
            jj[k] = ix.fd(part_idx,1) - 1;
            kk[k] = ix.fd(part_idx,2) - 1;

            // y split 2nd vp
            yint2 = _0_5_p_k_fld * dxi[1];

            delta = ( yint2 - yint ) / (xbuf[1] - yint);
            xint2 =  -xint + (xbuf[0] - xint) * delta;
            zint2 =   zint + (xbuf[2] - zint) * delta;

            k=k+1;

            x0[k] = -xint;
            y0[k] =  yint;
            z0[k] =  zint;

            x1[k] = xint2;
            y1[k] = yint2;
            z1[k] = zint2;

            qq[k] = q.fd(part_idx);

            ii[k] = ix.fd(part_idx,0) + dxi[0] - 1;
            jj[k] = ix.fd(part_idx,1) - 1;
            kk[k] = ix.fd(part_idx,2) - 1;

            k=k+1;

            x0[k] = xint2;
            y0[k] = -yint2;
            z0[k] = zint2;

            x1[k] = xbuf[0] - dxi[0];
            y1[k] = xbuf[1] - dxi[1];
            z1[k] = xbuf[2];

            qq[k] = q.fd(part_idx);

            ii[k] = ix.fd(part_idx,0) + dxi[0] - 1;
            jj[k] = ix.fd(part_idx,1) + dxi[1] - 1;
            kk[k] = ix.fd(part_idx,2) - 1;

            k=k+1;

          } else {

            // yspi = 2

            // y split 1st vp
            yint2 = _0_5_p_k_fld * dxi[1];
            delta = ( yint2 - x.fd(part_idx,1) ) / ( yint - x.fd(part_idx,1));
            xint2 = x.fd(part_idx,0) + (xint - x.fd(part_idx,0)) * delta;
            zint2 = x.fd(part_idx,2) + (zint - x.fd(part_idx,2)) * delta;

            x0[k] = x.fd(part_idx,0);
            y0[k] = x.fd(part_idx,1);
            z0[k] = x.fd(part_idx,2);

            x1[k] = xint2;
            y1[k] = yint2;
            z1[k] = zint2;

            qq[k] = q.fd(part_idx);

            ii[k] = ix.fd(part_idx,0) - 1;
            jj[k] = ix.fd(part_idx,1) - 1;
            kk[k] = ix.fd(part_idx,2) - 1;

            k=k+1;

            x0[k] = xint2;
            y0[k] = -yint2;
            z0[k] = zint2;

            x1[k] = xint;
            y1[k] = yint - dxi[1];
            z1[k] = zint;

            qq[k] = q.fd(part_idx);

            ii[k] = ix.fd(part_idx,0) - 1;
            jj[k] = ix.fd(part_idx,1) + dxi[1] - 1;
            kk[k] = ix.fd(part_idx,2) - 1;

            k=k+1;

            x0[k] = -xint;
            y0[k] = yint - dxi[1];
            z0[k] = zint;

            x1[k] = xbuf[0] - dxi[0];
            y1[k] = xbuf[1] - dxi[1];
            z1[k] = xbuf[2];

            qq[k] = q.fd(part_idx);

            ii[k] = ix.fd(part_idx,0) + dxi[0] - 1;
            jj[k] = ix.fd(part_idx,1) + dxi[1] - 1;
            kk[k] = ix.fd(part_idx,2) - 1;

            k=k+1;

          }

          break;

        // z cross only
        case 4:

          zint = _0_5_p_k_fld * dxi[2];
          delta = ( zint - x.fd(part_idx,2) ) / ( xbuf[2] - x.fd(part_idx,2));
          xint =  x.fd(part_idx,0) + (xbuf[0] - x.fd(part_idx,0)) * delta;
          yint =  x.fd(part_idx,1) + (xbuf[1] - x.fd(part_idx,1)) * delta;

          x0[k] = x.fd(part_idx,0);
          y0[k] = x.fd(part_idx,1);
          z0[k] = x.fd(part_idx,2);

          x1[k] = xint;
          y1[k] = yint;
          z1[k] = zint;

          qq[k] = q.fd(part_idx);

          ii[k] = ix.fd(part_idx,0) - 1;
          jj[k] = ix.fd(part_idx,1) - 1;
          kk[k] = ix.fd(part_idx,2) - 1;

          k=k+1;

          x0[k] = xint;
          y0[k] = yint;
          z0[k] = -zint;

          x1[k] = xbuf[0];
          y1[k] = xbuf[1];
          z1[k] = xbuf[2] - dxi[2];

          qq[k] = q.fd(part_idx);

          ii[k] = ix.fd(part_idx,0) - 1;
          jj[k] = ix.fd(part_idx,1) - 1;
          kk[k] = ix.fd(part_idx,2) + dxi[2] - 1;

          k=k+1;

          break;

        // x, z split - 3 vp
        case 5:

          xint = _0_5_p_k_fld * dxi[0];
          delta = ( xint - x.fd(part_idx,0) ) / ( xbuf[0] - x.fd(part_idx,0));
          yint =  x.fd(part_idx,1) + ( xbuf[1] - x.fd(part_idx,1)) * delta;
          zint =  x.fd(part_idx,2) + ( xbuf[2] - x.fd(part_idx,2)) * delta;

          zint2 = _0_5_p_k_fld * dxi[2];

          if ( zint >= -_0_5_p_k_fld && zint < _0_5_p_k_fld ){

            // zspi = 1

            // no z cross on 1st vp
            x0[k] = x.fd(part_idx,0);
            y0[k] = x.fd(part_idx,1);
            z0[k] = x.fd(part_idx,2);

            x1[k] = xint;
            y1[k] = yint;
            z1[k] = zint;

            qq[k] = q.fd(part_idx);

            ii[k] = ix.fd(part_idx,0) - 1;
            jj[k] = ix.fd(part_idx,1) - 1;
            kk[k] = ix.fd(part_idx,2) - 1;

            // z split 2nd vp
            delta = ( zint2 - zint ) / (xbuf[2] - zint);
            xint2 =  -xint + (xbuf[0] - xint) * delta;
            yint2 =   yint + (xbuf[1] - yint) * delta;

            k=k+1;

            x0[k] = -xint;
            y0[k] =  yint;
            z0[k] =  zint;

            x1[k] = xint2;
            y1[k] = yint2;
            z1[k] = zint2;

            qq[k] = q.fd(part_idx);

            ii[k] = ix.fd(part_idx,0) + dxi[0] - 1;
            jj[k] = ix.fd(part_idx,1) - 1;
            kk[k] = ix.fd(part_idx,2) - 1;

            k=k+1;

            x0[k] = xint2;
            y0[k] = yint2;
            z0[k] = -zint2;

            x1[k] = xbuf[0] - dxi[0];
            y1[k] = xbuf[1];
            z1[k] = xbuf[2] - dxi[2];

            qq[k] = q.fd(part_idx);

            ii[k] = ix.fd(part_idx,0) + dxi[0] - 1;
            jj[k] = ix.fd(part_idx,1) - 1;
            kk[k] = ix.fd(part_idx,2) + dxi[2] - 1;

            k=k+1;

          } else {

            // zspi = 2

            // z split 1st vp
            delta = ( zint2 - x.fd(part_idx,2) ) / ( zint - x.fd(part_idx,2));
            xint2 = x.fd(part_idx,0) + (xint - x.fd(part_idx,0)) * delta;
            yint2 = x.fd(part_idx,1) + (yint - x.fd(part_idx,1)) * delta;

            x0[k] = x.fd(part_idx,0);
            y0[k] = x.fd(part_idx,1);
            z0[k] = x.fd(part_idx,2);

            x1[k] = xint2;
            y1[k] = yint2;
            z1[k] = zint2;

            qq[k] = q.fd(part_idx);

            ii[k] = ix.fd(part_idx,0) - 1;
            jj[k] = ix.fd(part_idx,1) - 1;
            kk[k] = ix.fd(part_idx,2) - 1;

            k=k+1;

            x0[k] = xint2;
            y0[k] = yint2;
            z0[k] = -zint2;

            x1[k] = xint;
            y1[k] = yint;
            z1[k] = zint - dxi[2];

            qq[k] = q.fd(part_idx);

            ii[k] = ix.fd(part_idx,0) - 1;
            jj[k] = ix.fd(part_idx,1) - 1;
            kk[k] = ix.fd(part_idx,2) + dxi[2] - 1;

            k=k+1;

            x0[k] = -xint;
            y0[k] = yint;
            z0[k] = zint - dxi[2];

            x1[k] = xbuf[0] - dxi[0];
            y1[k] = xbuf[1];
            z1[k] = xbuf[2] - dxi[2];

            qq[k] = q.fd(part_idx);

            ii[k] = ix.fd(part_idx,0) + dxi[0] - 1;
            jj[k] = ix.fd(part_idx,1) - 1;
            kk[k] = ix.fd(part_idx,2) + dxi[2] - 1;

            k=k+1;

          }

          break;

        // y, z split - 3 vp
        case 6:

          yint = _0_5_p_k_fld * dxi[1];
          delta = ( yint - x.fd(part_idx,1) ) / ( xbuf[1] - x.fd(part_idx,1));
          xint =  x.fd(part_idx,0) + ( xbuf[0] - x.fd(part_idx,0)) * delta; // no x cross
          zint =  x.fd(part_idx,2) + ( xbuf[2] - x.fd(part_idx,2)) * delta;

          if ( zint >= -_0_5_p_k_fld && zint < _0_5_p_k_fld ){

            // zspi = 1

            // no z cross on 1st vp
            x0[k] = x.fd(part_idx,0);
            y0[k] = x.fd(part_idx,1);
            z0[k] = x.fd(part_idx,2);

            x1[k] = xint;
            y1[k] = yint;
            z1[k] = zint;

            qq[k] = q.fd(part_idx);

            ii[k] = ix.fd(part_idx,0) - 1;
            jj[k] = ix.fd(part_idx,1) - 1;
            kk[k] = ix.fd(part_idx,2) - 1;

            // z split 2nd vp
            zint2 = _0_5_p_k_fld * dxi[2];
            delta = ( zint2 - zint ) / (xbuf[2] - zint);
            xint2 =  xint + (xbuf[0] - xint) * delta;
            yint2 =  -yint + (xbuf[1] - yint) * delta;

            k=k+1;

            x0[k] = xint;
            y0[k] = -yint;
            z0[k] = zint;

            x1[k] = xint2;
            y1[k] = yint2;
            z1[k] = zint2;

            qq[k] = q.fd(part_idx);

            ii[k] = ix.fd(part_idx,0) - 1;
            jj[k] = ix.fd(part_idx,1) + dxi[1] - 1;
            kk[k] = ix.fd(part_idx,2) - 1;

            k=k+1;

            x0[k] = xint2;
            y0[k] = yint2;
            z0[k] = -zint2;

            x1[k] = xbuf[0];
            y1[k] = xbuf[1] - dxi[1];
            z1[k] = xbuf[2] - dxi[2];

            qq[k] = q.fd(part_idx);

            ii[k] = ix.fd(part_idx,0) - 1;
            jj[k] = ix.fd(part_idx,1) + dxi[1] - 1;
            kk[k] = ix.fd(part_idx,2) + dxi[2] - 1;

            k=k+1;

          } else {

            // zspi = 2

            // z split 1st vp
            zint2 = _0_5_p_k_fld * dxi[2];
            delta = ( zint2 - x.fd(part_idx,2) ) / ( zint - x.fd(part_idx,2));

            xint2 = x.fd(part_idx,0) + (xint - x.fd(part_idx,0)) * delta;
            yint2 = x.fd(part_idx,1) + (yint - x.fd(part_idx,1)) * delta;

            x0[k] = x.fd(part_idx,0);
            y0[k] = x.fd(part_idx,1);
            z0[k] = x.fd(part_idx,2);

            x1[k] = xint2;
            y1[k] = yint2;
            z1[k] = zint2;

            qq[k] = q.fd(part_idx);

            ii[k] = ix.fd(part_idx,0) - 1;
            jj[k] = ix.fd(part_idx,1) - 1;
            kk[k] = ix.fd(part_idx,2) - 1;

            k=k+1;

            x0[k] = xint2;
            y0[k] = yint2;
            z0[k] = -zint2;

            x1[k] = xint;
            y1[k] = yint;
            z1[k] = zint - dxi[2];

            qq[k] = q.fd(part_idx);

            ii[k] = ix.fd(part_idx,0) - 1;
            jj[k] = ix.fd(part_idx,1) - 1;
            kk[k] = ix.fd(part_idx,2) + dxi[2] - 1;

            k=k+1;

            x0[k] = xint;
            y0[k] = -yint;
            z0[k] = zint - dxi[2];

            x1[k] = xbuf[0];
            y1[k] = xbuf[1] - dxi[1];
            z1[k] = xbuf[2] - dxi[2];

            qq[k] = q.fd(part_idx);

            ii[k] = ix.fd(part_idx,0) - 1;
            jj[k] = ix.fd(part_idx,1) + dxi[1] - 1;
            kk[k] = ix.fd(part_idx,2) + dxi[2] - 1;

            k=k+1;

          }

          break;

        // x, y, z cross
        case 7:

          xint = _0_5_p_k_fld * dxi[0];
          delta = ( xint - x.fd(part_idx,0) ) / ( xbuf[0] - x.fd(part_idx,0));
          yint =  x.fd(part_idx,1) + ( xbuf[1] - x.fd(part_idx,1)) * delta;
          zint =  x.fd(part_idx,2) + ( xbuf[2] - x.fd(part_idx,2)) * delta;

          if ( yint >= -_0_5_p_k_fld && yint < _0_5_p_k_fld ){

            // no y cross on 1st vp
            x0[k] = x.fd(part_idx,0);
            y0[k] = x.fd(part_idx,1);
            z0[k] = x.fd(part_idx,2);

            x1[k] = xint;
            y1[k] = yint;
            z1[k] = zint;

            qq[k] = q.fd(part_idx);

            ii[k] = ix.fd(part_idx,0) - 1;
            jj[k] = ix.fd(part_idx,1) - 1;
            kk[k] = ix.fd(part_idx,2) - 1;

            // y split 2nd vp
            yint2 = _0_5_p_k_fld * dxi[1];

            delta = ( yint2 - yint ) / (xbuf[1] - yint);
            xint2 =  -xint + (xbuf[0] - xint) * delta;
            zint2 =   zint + (xbuf[2] - zint) * delta;

            k=k+1;

            x0[k] = -xint;
            y0[k] =  yint;
            z0[k] =  zint;

            x1[k] = xint2;
            y1[k] = yint2;
            z1[k] = zint2;

            qq[k] = q.fd(part_idx);

            ii[k] = ix.fd(part_idx,0) + dxi[0] - 1;
            jj[k] = ix.fd(part_idx,1) - 1;
            kk[k] = ix.fd(part_idx,2) - 1;

            k=k+1;

            x0[k] = xint2;
            y0[k] = -yint2;
            z0[k] = zint2;

            x1[k] = xbuf[0] - dxi[0];
            y1[k] = xbuf[1] - dxi[1];
            z1[k] = xbuf[2];

            qq[k] = q.fd(part_idx);

            ii[k] = ix.fd(part_idx,0) + dxi[0] - 1;
            jj[k] = ix.fd(part_idx,1) + dxi[1] - 1;
            kk[k] = ix.fd(part_idx,2) - 1;

            k=k+1;

          } else {

            // y split 1st vp
            yint2 = _0_5_p_k_fld * dxi[1];
            delta = ( yint2 - x.fd(part_idx,1) ) / ( yint - x.fd(part_idx,1));
            xint2 = x.fd(part_idx,0) + (xint - x.fd(part_idx,0)) * delta;
            zint2 = x.fd(part_idx,2) + (zint - x.fd(part_idx,2)) * delta;

            x0[k] = x.fd(part_idx,0);
            y0[k] = x.fd(part_idx,1);
            z0[k] = x.fd(part_idx,2);

            x1[k] = xint2;
            y1[k] = yint2;
            z1[k] = zint2;

            qq[k] = q.fd(part_idx);

            ii[k] = ix.fd(part_idx,0) - 1;
            jj[k] = ix.fd(part_idx,1) - 1;
            kk[k] = ix.fd(part_idx,2) - 1;

            k=k+1;

            x0[k] = xint2;
            y0[k] = -yint2;
            z0[k] = zint2;

            x1[k] = xint;
            y1[k] = yint - dxi[1];
            z1[k] = zint;

            qq[k] = q.fd(part_idx);

            ii[k] = ix.fd(part_idx,0) - 1;
            jj[k] = ix.fd(part_idx,1) + dxi[1] - 1;
            kk[k] = ix.fd(part_idx,2) - 1;

            k=k+1;

            x0[k] = -xint;
            y0[k] = yint - dxi[1];
            z0[k] = zint;

            x1[k] = xbuf[0] - dxi[0];
            y1[k] = xbuf[1] - dxi[1];
            z1[k] = xbuf[2];

            qq[k] = q.fd(part_idx);

            ii[k] = ix.fd(part_idx,0) + dxi[0] - 1;
            jj[k] = ix.fd(part_idx,1) + dxi[1] - 1;
            kk[k] = ix.fd(part_idx,2) - 1;

            k=k+1;

          }

          // one of the 3 vp requires an additional z split
          zint = _0_5_p_k_fld * dxi[2];

          for ( int vpidx=0; vpidx<3; vpidx++ ) {
            if ( z1[vpidx] < -_0_5_p_k_fld or z1[vpidx] >= _0_5_p_k_fld ){
              delta = ( zint - z0[vpidx] ) / ( z1[vpidx] - z0[vpidx] );
              xint =  x0[vpidx] + (x1[vpidx] - x0[vpidx]) * delta;
              yint =  y0[vpidx] + (y1[vpidx] - y0[vpidx]) * delta;

              x0[k] = xint;
              y0[k] = yint;
              z0[k] = -zint;

              x1[k] = x1[vpidx];
              y1[k] = y1[vpidx];
              z1[k] = z1[vpidx] - dxi[2];

              qq[k] = qq[vpidx];

              ii[k] = ii[vpidx];
              jj[k] = jj[vpidx];
              kk[k] = kk[vpidx] + dxi[2];

              // correct old vp
              x1[vpidx] = xint;
              y1[vpidx] = yint;
              z1[vpidx] = zint;

              // correct remaining vp
              for ( int jjj=vpidx+1; jjj<k; jjj++ ) {
                z0[jjj] -= dxi[2];
                z1[jjj] -= dxi[2];
                kk[jjj] += dxi[2];
              }

              k=k+1;

              break;
            }
          }

          break;

      }
      __syncwarp();

      int nsplit = k;

      // --------------

      // now accumulate jay iooping through all virtual particles
      // and shifting grid indexes

      for (int i=0; i<nsplit; ++i){

        // Normalize charge
        qnx = qq[i] * jnorm1;
        qny = qq[i] * jnorm2;
        qnz = qq[i] * jnorm3;

        // get spline weights for x, y and z
        SPLINE( x0[i], S0x );
        SPLINE( x1[i], S1x );

        SPLINE( y0[i], S0y );
        SPLINE( y1[i], S1y );

        SPLINE( z0[i], S0z );
        SPLINE( z1[i], S1z );

        // get longitudinal motion weights
        // the last value is set to 0 so we can accumulate
        // all current components in a single pass
        WL( qnx, x0[i], x1[i], wl1 );
        wl1[P1] = 0;

        WL( qny, y0[i], y1[i], wl2 );
        wl2[P1] = 0;

        WL( qnz, z0[i], z1[i], wl3 );
        wl3[P1] = 0;

        // get perpendicular motion weights
        // (a division by 2 was moved to the jnorm vars. above)
        for (int k2=P0; k2<=P1; ++k2){
          for (int k1=P0; k1<=P1; ++k1){
            wp1[k1+k2*NB] = (S0y[k1]*S0z[k2] + S1y[k1]*S1z[k2]) +
                         _0_5_p_k_fld * (S0y[k1]*S1z[k2] + S1y[k1]*S0z[k2]);
            wp2[k1+k2*NB] = (S0x[k1]*S0z[k2] + S1x[k1]*S1z[k2]) +
                         _0_5_p_k_fld * (S0x[k1]*S1z[k2] + S1x[k1]*S0z[k2]);
            wp3[k1+k2*NB] = (S0x[k1]*S0y[k2] + S1x[k1]*S1y[k2]) +
                         _0_5_p_k_fld * (S0x[k1]*S1y[k2] + S1x[k1]*S0y[k2]);
          }
        }

        // accumulate j1, j2, j3 in a single pass
        for (int k3=P0; k3<=P1; ++k3){
          for (int k2=P0; k2<=P1; ++k2){
            for (int k1=P0; k1<=P1; ++k1){

              atomicAdd( jay_f3_til.fdp(0,ii[i]+k1,jj[i]+k2,kk[i]+k3), wl1[k1] * wp1[k2+k3*NB] );
              atomicAdd( jay_f3_til.fdp(1,ii[i]+k1,jj[i]+k2,kk[i]+k3), wl2[k2] * wp2[k1+k3*NB] );
              atomicAdd( jay_f3_til.fdp(2,ii[i]+k1,jj[i]+k2,kk[i]+k3), wl3[k3] * wp3[k1+k2*NB] );

            }
          }
        }

      }

      // --------------

      x.fd(part_idx,0) = xbuf[0] - dxi[0];
      x.fd(part_idx,1) = xbuf[1] - dxi[1];
      x.fd(part_idx,2) = xbuf[2] - dxi[2];
      ix.fd(part_idx,0) += dxi[0];
      ix.fd(part_idx,1) += dxi[1];
      ix.fd(part_idx,2) += dxi[2];

      // Note: this only works if chunk size divisible by block size (which we enforce at initialization)
      // See sort for how to do this more generally, if we decide to do that eventually
      part_idx += BLOCK_SIZE;

    } // end loop over chunk of particles ``while ( part_idx < num_par_proc_chunk )``

    num_par_proc_tot -= num_par_proc_chunk;
    chunk++;

  } // end loop over chunks  ``while ( num_par_remainin > 0 )``

  // one thread on each tile zero out counters to be used in sort
  if ( tx==0 and block_start_chunks[bx]==0 ) {
    n_hole[til] = 0;
    n_comm_ext[til] = 0;
    n_recv_loc[til] = 0;
  }

#undef NB

}




#undef DUDT_BORIS_1D
#undef DUDT_BORIS_2D
#undef DUDT_BORIS_3D
#undef DUDT_BORIS_1D_CS
#undef DUDT_BORIS_2D_CS
#undef DUDT_BORIS_3D_CS
#undef SPLINE
#undef SPLINEH
#undef LP
#undef UP
#undef WL
#undef P0
#undef P1
#undef ADVANCE_DEPOSIT_1D
#undef ADVANCE_DEPOSIT_2D
#undef ADVANCE_DEPOSIT_3D

#endif

#ifndef __TEMPLATE__

// ---------------------------------------------------------------------------------------
// Function stub to call kernel
// ---------------------------------------------------------------------------------------
void dudt_boris_wrapper( t_chunk ***chunks, int *block_start_tiles,
                         int *block_start_chunks, int *block_num_par,
                         p_k_part rqm, t_emf_c *emf, int interpolation, bool grid_center,
                         p_double *energy, bool report_energy, p_double dt,
                         int nblocks, int chunk_size, cudaStream_t stream ){

  if (grid_center) {
    // Interpolate fields centered at the cell corner
    switch( interpolation ) {
      case p_linear:
        switch( p_x_dim ) {
          case 1:
            dudt_boris_1d_cs1<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd1, emf->d_ed1, energy, report_energy, dt, chunk_size );
            break;
          case 2:
            dudt_boris_2d_cs1<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd2, emf->d_ed2, energy, report_energy, dt, chunk_size );
            break;
          case 3:
            dudt_boris_3d_cs1<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd3, emf->d_ed3, energy, report_energy, dt, chunk_size );
            break;
        }
        break;
      case p_quadratic:
        switch( p_x_dim ) {
          case 1:
            dudt_boris_1d_cs2<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd1, emf->d_ed1, energy, report_energy, dt, chunk_size );
            break;
          case 2:
            dudt_boris_2d_cs2<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd2, emf->d_ed2, energy, report_energy, dt, chunk_size );
            break;
          case 3:
            dudt_boris_3d_cs2<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd3, emf->d_ed3, energy, report_energy, dt, chunk_size );
            break;
        }
        break;
      case p_cubic:
        switch( p_x_dim ) {
          case 1:
            dudt_boris_1d_cs3<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd1, emf->d_ed1, energy, report_energy, dt, chunk_size );
            break;
          case 2:
            dudt_boris_2d_cs3<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd2, emf->d_ed2, energy, report_energy, dt, chunk_size );
            break;
          case 3:
            dudt_boris_3d_cs3<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd3, emf->d_ed3, energy, report_energy, dt, chunk_size );
            break;
        }
        break;
      case p_quartic:
        switch( p_x_dim ) {
          case 1:
            dudt_boris_1d_cs4<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd1, emf->d_ed1, energy, report_energy, dt, chunk_size );
            break;
          case 2:
            dudt_boris_2d_cs4<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd2, emf->d_ed2, energy, report_energy, dt, chunk_size );
            break;
          case 3:
            dudt_boris_3d_cs4<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd3, emf->d_ed3, energy, report_energy, dt, chunk_size );
            break;
        }
        break;
    }
  } else {
    // Interpolate fields on a staggered Yee mesh
    switch( interpolation ) {
      case p_linear:
        switch( p_x_dim ) {
          case 1:
            dudt_boris_1d_s1<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd1, emf->d_ed1, energy, report_energy, dt, chunk_size );
            break;
          case 2:
            dudt_boris_2d_s1<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd2, emf->d_ed2, energy, report_energy, dt, chunk_size );
            break;
          case 3:
            dudt_boris_3d_s1<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd3, emf->d_ed3, energy, report_energy, dt, chunk_size );
            break;
        }
        break;
      case p_quadratic:
        switch( p_x_dim ) {
          case 1:
            dudt_boris_1d_s2<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd1, emf->d_ed1, energy, report_energy, dt, chunk_size );
            break;
          case 2:
            dudt_boris_2d_s2<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd2, emf->d_ed2, energy, report_energy, dt, chunk_size );
            break;
          case 3:
            dudt_boris_3d_s2<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd3, emf->d_ed3, energy, report_energy, dt, chunk_size );
            break;
        }
        break;
      case p_cubic:
        switch( p_x_dim ) {
          case 1:
            dudt_boris_1d_s3<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd1, emf->d_ed1, energy, report_energy, dt, chunk_size );
            break;
          case 2:
            dudt_boris_2d_s3<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd2, emf->d_ed2, energy, report_energy, dt, chunk_size );
            break;
          case 3:
            dudt_boris_3d_s3<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd3, emf->d_ed3, energy, report_energy, dt, chunk_size );
            break;
        }
        break;
      case p_quartic:
        switch( p_x_dim ) {
          case 1:
            dudt_boris_1d_s4<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd1, emf->d_ed1, energy, report_energy, dt, chunk_size );
            break;
          case 2:
            dudt_boris_2d_s4<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd2, emf->d_ed2, energy, report_energy, dt, chunk_size );
            break;
          case 3:
            dudt_boris_3d_s4<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( chunks,
                        block_start_tiles, block_start_chunks, block_num_par, rqm,
                        emf->d_bd3, emf->d_ed3, energy, report_energy, dt, chunk_size );
            break;
        }
        break;
    }
  }

  chkLastErr("dudt_boris");
#ifdef __DEBUG__
  chkErr(cudaDeviceSynchronize());
#endif

}


// ---------------------------------------------------------------------------------------
// Function stub to call kernel
// ---------------------------------------------------------------------------------------
void advance_deposit_wrapper( t_chunk ***chunks, int *block_start_tiles,
                              int *block_start_chunks, int *block_num_par,
                              p_double *dx, t_current_c *jay, int interpolation,
                              p_double dt, int nblocks,
                              int chunk_size, int *n_hole, int *n_comm_ext,
                              int *n_recv_loc, cudaStream_t stream ){

  switch( interpolation ) {
    case p_linear:
      switch( p_x_dim ) {
        case 1:
          advance_deposit_1d_s1<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( chunks,
                                  block_start_tiles, block_start_chunks, block_num_par,
                                  dx[0], jay->d_jd1, dt, chunk_size,
                                  n_hole, n_comm_ext, n_recv_loc );
          break;
        case 2:
          advance_deposit_2d_s1<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( chunks,
                                  block_start_tiles, block_start_chunks, block_num_par,
                                  dx[0], dx[1], jay->d_jd2, dt, chunk_size,
                                  n_hole, n_comm_ext, n_recv_loc );
          break;
        case 3:
          advance_deposit_3d_s1<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( chunks,
                                  block_start_tiles, block_start_chunks, block_num_par,
                                  dx[0], dx[1], dx[2], jay->d_jd3, dt, chunk_size,
                                  n_hole, n_comm_ext, n_recv_loc );
          break;
      }
      break;
    case p_quadratic:
      switch( p_x_dim ) {
        case 1:
          advance_deposit_1d_s2<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( chunks,
                                  block_start_tiles, block_start_chunks, block_num_par,
                                  dx[0], jay->d_jd1, dt, chunk_size,
                                  n_hole, n_comm_ext, n_recv_loc );
          break;
        case 2:
          advance_deposit_2d_s2<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( chunks,
                                  block_start_tiles, block_start_chunks, block_num_par,
                                  dx[0], dx[1], jay->d_jd2, dt, chunk_size,
                                  n_hole, n_comm_ext, n_recv_loc );
          break;
        case 3:
          advance_deposit_3d_s2<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( chunks,
                                  block_start_tiles, block_start_chunks, block_num_par,
                                  dx[0], dx[1], dx[2], jay->d_jd3, dt, chunk_size,
                                  n_hole, n_comm_ext, n_recv_loc );
          break;
      }
      break;
    case p_cubic:
      switch( p_x_dim ) {
        case 1:
          advance_deposit_1d_s3<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( chunks,
                                  block_start_tiles, block_start_chunks, block_num_par,
                                  dx[0], jay->d_jd1, dt, chunk_size,
                                  n_hole, n_comm_ext, n_recv_loc );
          break;
        case 2:
          advance_deposit_2d_s3<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( chunks,
                                  block_start_tiles, block_start_chunks, block_num_par,
                                  dx[0], dx[1], jay->d_jd2, dt, chunk_size,
                                  n_hole, n_comm_ext, n_recv_loc );
          break;
        case 3:
          advance_deposit_3d_s3<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( chunks,
                                  block_start_tiles, block_start_chunks, block_num_par,
                                  dx[0], dx[1], dx[2], jay->d_jd3, dt, chunk_size,
                                  n_hole, n_comm_ext, n_recv_loc );
          break;
      }
      break;
    case p_quartic:
      switch( p_x_dim ) {
        case 1:
          advance_deposit_1d_s4<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( chunks,
                                  block_start_tiles, block_start_chunks, block_num_par,
                                  dx[0], jay->d_jd1, dt, chunk_size,
                                  n_hole, n_comm_ext, n_recv_loc );
          break;
        case 2:
          advance_deposit_2d_s4<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( chunks,
                                  block_start_tiles, block_start_chunks, block_num_par,
                                  dx[0], dx[1], jay->d_jd2, dt, chunk_size,
                                  n_hole, n_comm_ext, n_recv_loc );
          break;
        case 3:
          advance_deposit_3d_s4<__BLOCK_SIZE__> <<< nblocks, __BLOCK_SIZE__, 0, stream >>>( chunks,
                                  block_start_tiles, block_start_chunks, block_num_par,
                                  dx[0], dx[1], dx[2], jay->d_jd3, dt, chunk_size,
                                  n_hole, n_comm_ext, n_recv_loc );
          break;
      }
      break;
  }

  chkLastErr("advance_deposit");
#ifdef __DEBUG__
  chkErr(cudaDeviceSynchronize());
#endif

}

#endif
