# CUDA algorithm

## Overview

The CUDA module enables one to run OSIRIS on NVIDIA GPU accelerated architectures. It implements a GPU parallelized particle pusher and current deposit in CUDA, while relying on the existing CPU routines for the field solve, particle communication, particle boundary conditions, field boundary conditions, and particle initialization. This approach represents a compromise between computational performance and implementation complexity. Porting every feature of OSIRIS to the GPU would have been an onerous task tantamount to writing an entirely new code. Thus, we chose to only accelerate the most computationally expensive parts of the PIC algorithm, which are the particle pusher and current deposit. Particle sorting is also done on the GPU so that particles can live in device memory most of the time and expensive copy operations between device and host can be kept to a minimum.

CUDA is implemented as a simulation mode and is built on top of the tiles data structures/simulation mode. For this reason, some of the features and parameters used here will be new to a user of "standard" OSIRIS. See `source/tiles/README.md` for a full description. The optimal number of tiles in a simulation can depend on the physical setup as well as the computing hardware. In benchmarks on thermal plasmas, the optimal number of gridpoints per tile was found to be $\sim 32768 = 32^3 = 256 \times 128$ in both 2D and 3D; this would be a good starting point for most simulations. The default tile_topology is "snake", but the "hilbert" topology can offer better load balance at the cost of increased MPI communications. Thermal benchmarks also indicate that at least $2^{18}$ total gridpoints should be used per GPU for optimal push performance regardless of particles per cell.

To efficiently store particles in GPU memory, each MPI rank stores particles on its respective GPU in "chunks" of memory which store data for a specified number of particles. Each rank allocates most of the device memory at initialization as chunks, which can then be freely grabbed from and returned to the global chunk pool as necessary with little overhead. Each tile maintains a "pad" of enough extra chunks of memory such that buffer overflow during communication between tiles is unlikely, but few enough such that each tile nearly saturates its memory capacity (i.e. wastes little memory). Ideally, this results in maximal device memory utilization: each tile has just the right amount of memory, and there is nearly zero overhead to effectively resizing its memory by grabbing or returning chunks.

In practice, the chunk pool is controlled by two parameters: `chunk_size`, which is the number of particles each chunk can store, and `chunk_pool_size`, which is the total number of chunks allocated by each MPI rank. Using a small chunk size results in greater flexibility for memory management, but additional overhead for grabbing and returning chunks from the pool. We have found that the performance penalty from the latter is neglibigle; thus, we recommend using the smallest possible chunk size of 64. On the other hand, `chunk_pool_size` must be chosen such that `chunk_size` * `chunk_pool_size` $\geq$ N, where N is the number of particles on each GPU. Note that with dynamic load balancing, N can be assumed to be approximately constant. In practice, one can simply set the `chunk_size` and `chunk_pool_size` such that chunks take up most of the available memory on the GPU. The default is that the `chunk_pool_size`, for given `chunk_size`, takes up 75% of GPU memory, which should be appropriate for most applications.

As mentioned above, in addition to holding enough chunks to fit all its particles, each tile holds some additional chunks for buffering particles during the particle sort. Specifically, each tile holds `nchunks_pad` chunks to buffer particles coming from tiles located on the same MPI rank, and another `nchunks_pad` chunks to buffer particles going to or coming from tiles located on different MPI ranks. The user can specify this number directly, or (recommended for most cases) can specify this as a fraction of the number of chunks per tile via `nchunks_pad_frac`. See description of each below.

Related parameters are `ihole_size_frac` and `ihole_size`, which control the size of a buffer called "ihole" [1] used in the particle sorting algorithm. `ihole_size_frac` represents the number of particles which fit into ihole as a fraction of the total particle capacity of the chunk pool. For example, if `ihole_size_frac = 0.05`, then additional memory will be allocated to fit 5% of the particle capacity of the chunk pool; this memory is equally distributed among all tiles. The ihole itself does not use chunks, but because the required ihole size and `chunk_pool_size` both scale with the total number of particles, it is convenient to represent the size of the ihole as a fraction of the total chunk pool memory. The default value for `ihole_size_frac` of 0.05 is sufficient for most simulations. Alternatively, `ihole_size` can be specified. This parameter overrides `ihole_size_frac` to set the size of the ihole in absolute number of particles for each tile.

In the current implementation, if the ihole or chunks overflow during the particle sort, then the simulation crashes and prints an error informing the user to increase either `ihole_size_frac` or `nchunks_pad_frac` (or the absolute equivalents). In principle, overflows can be gracefully handled by using free chunks from the pool as additional memory within the particle sort. Normally, device memory cannot be reallocated with a CUDA kernel, but because the chunks are preallocated, they can be accessed and used arbitrarily. However, this form of dynamic error handling has not yet been implemented in OSIRIS.

In order to expedite communication for transfers between host and device, the code uses a "memory transfer manager" for batched data movements. This object buffers data and metadata for data movements in two fixed-size communication buffers allocated at initialization. One buffer resides on the host, the other equal in size on the device. The sending buffer is copied to the receiving buffer for movement between host and device whenever the sending buffer becomes full;  device-to-device copies between the device buffer and the original device memory addresses are performed via a kernel regardless of the movement direction. This construction greatly reduces the number of copies between host and device. The `com_buf_size` parameter controls the amount of memory allocated for the communication buffer in megabytes. The default
value of 1GB should be sufficient for most applications.


## Using CUDA

To use CUDA, OSIRIS must be compiled CUDA support by setting `ENABLE_CUDA = 1` in the configuration script. See config/osiris_sys.perlmutter for an example of this and additional configuration options required.

### Input parameters

* The cuda algorithm is currently implemented as a simulation mode. To activate it the following parameter must be added to the `simulation` section:
  * `algorithm=cuda,`

* The following parameters have been added to the `node_conf` section:
  * `chunk_size`, integer, default = 64. Size of a chunk in units of number of particles.
  
  * `chunk_pool_size`, integer, default = -1. Number of chunks allocated per MPI rank, in units of number of chunks. If not specified 75% of the available memory on the GPU will be allocated for chunks. This should be appropriate for most applications.
  
  * `com_buf_size`, integer, default = 1024. Size of buffer for CPU-GPU communication allocated for each MPI rank, in units of MB (so default is about 1GB). The default should be sufficient for most applications.
  
  * `ihole_size_frac`, real, default = 0.05. Amount of memory allocated for the ihole buffers for each tile specified as a fraction of the particle capacity of the chunk pool (`chunk_size`*`chunk_pool_size`). (Since this memory is distributed equally between all tiles, this parameter can also be thought of in terms of the approximate number of particles per tile.)
  
  * `ihole_size`, integer, default = -1. Optionally, one can specify the size of the ihole buffer for each tile in units of number of particles. If set, this overrides the value of `ihole_size_frac`.
  
  * `nchunks_pad_frac`, real, default = 0.05. Specifies, as a fraction of number of chunks per tile, (1) how many empty chunks each tile should have to buffer particles moving locally during the sort, and (2) how many chunks each tile should have to buffer particles moving externally (via MPI) during the sort.
  
  * `nchunks_pad`, integer, default = -1. Optionally, one can specify this in units of number of chunks. If set, this overrides the value of `nchunks_pad_frac`. 
  
  * `nblocks_requested`, integer, default = 10000. Number of GPU compute blocks requested. This default value should be appropriate for most applications.
  
  * `block_size`, integer, default = 32. Specifies the number of threads per block for GPU kernel launches. Must divide chunk_size. Currently only `block_size = 32,` has been implemented.

  * Note on setting `node_number`: the product of this number gives the total number of MPI ranks for the simulation. MPI ranks will be assigned to physical GPUs depending on how many are available. While it is possible to overload GPUs by setting the total number of MPI ranks to greater than the number of physical GPUs available, it is recommended that the user set the number of MPI ranks to less than or equal to the number of physical GPUs available (ie such that each MPI rank has it's own GPU). This will give best performance. Also, the default parameters above are set under this assumption.

See `decks/cuda/` for example input decks.

## References

[1] Decyk, Viktor K., and Tajendra V. Singh. "Particle-in-cell algorithms for emerging computer architectures." *Computer Physics Communications* 185.3 (2014): 708-719.
