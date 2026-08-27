# This is meant to be included by the osiris Makefile
ifdef ENABLE_CUDA

# Cuda algorithm source files
#src += os-simulation-cuda-define.f03 os-cuda-serialize-state.f03 os-simulation-cuda.f03
src += simple-cpp-profiler.cpp os-cuda-input-deck.cpp os-cuda-state.cpp os-cuda-device-runtime-interaction.cpp os-cuda-main.cpp

# Add local path to vpath
vpath %.f03 cuda
vpath %.cpp cuda
vpath %.h cuda
vpath %.cu cuda

# Specify dependencies and additional flags
# Needed for parsing of the Cuda part of the input deck.
INCPATH += -I$(builddir)/opencontainers_1_7_6/include
FPPF += -DENABLE_CUDA

# tell the final linking step to link in the CPP libraries
#	TODO: make this more robust.. I should just compile all the
#	Cuda/c++ stuff into a library will all Cuda/C++ dependenies resolved...
#   then use the library at the top level so as to not compilcate the linking process.
ifeq ("$(F90)","ifort")
	LDF += -cxxlib
endif

CPPF += $(FPPF)

# create the CUDA (nvcc) compile flags.
CUDA_CF  = $(CUDAFLAGS_$(COMPILATION_TYPE))
CUDA_CF += -DENABLE_CUDA
CUDA_CF += -D$(UNDERSCORE)
CUDA_CF += $(FPPF)

# Create/Fill-in needed CUDA configuration variables
# First specify the CUDA compilter... if not specfied, set to 'nvcc'
ifeq ($(strip $(NVCC)),)
	NVCC = nvcc
endif

# Find CUDA by using the linux 'which' command to see where the
# 	CUDA compiler 'nvcc' resides. 'which' returns the full path
# 	so use the 'dir' Make command to just substring out the path portion
# 	All the CUDA stuff is underneath that directory in ./include or ./lib or ./lib64
# 	directories.. Also some openGL and glut stuff is there under the cuda_samples directory.
ifeq ($(strip $(CUDABASEDIR)),)
PATH_TO_NVCC = $(shell which $(NVCC) 2> null )
endif

# now remove the bin/$(NVCC) (which usually evaluates to bin/nvcc ) at the end of $(PATH_TO_NVCC)
CUDABASEDIR = $(strip $(subst /bin/$(NVCC),,$(PATH_TO_NVCC) ) )

# Check to see that we actually found it...
ifeq ($(strip $(CUDABASEDIR)),)
$(error The CUDA compiler $(NVCC) as specified in the Makefile variable NVCC was not found in the path.)
endif

# Add include path to CUDA headers 
INCPATH += -I$(CUDABASEDIR)/include

# Library stuff.. TODO: make this more configurable
CUDA_LINKFLAGS = -L$(CUDABASEDIR)/lib64 -lcudart -lcuda -lcurand
LDF += $(CUDA_LINKFLAGS)

# Make a (hopefully) informative banner about the setting we will use.
$(info )
$(info ############################################################################## )
$(info -                            CUDA Enabled                                    )
$(info -  CUDA compiler: $(NVCC)                                                    )
$(info -  CUDA compiler found at: )
$(info -       $(PATH_TO_NVCC))
$(info -  CUDA library base path at:)
$(info -       $(CUDABASEDIR))
$(info -  CUDA library linker flags:)
$(info -       $(CUDA_LINKFLAGS))
$(info ############################################################################## )
$(info )

# add needed target rules for c++ and Cuda

# this rule untars the 'opencontainers' header-only library that reads python-style configuration
#	files.. it is used by the c++ cuda code to read any CUDA specific settings in the input deck.
$(builddir)/opencontainers_1_7_6: | $(builddir)
	@echo "Untaring supporting header libraries"
	$(shell tar xf ../tools/misc/cuda/opencontainers.tar -C $(builddir))

# rule to make sure the needed libraries are untared and ready before compiling certain files
os-cuda-input-deck.cpp : | $(builddir)/opencontainers_1_7_6

# rule for .cu (cuda kernal) files..
# .cu (cuda kernal and normal c files)
%.o : ../source/%.cu
	$(NVCC) $(CUDA_CF) $(INCPATH) -c $< 2>&1 | c++filt | tee $<.txt

# If not at top level go there and make
ifeq ("$(toplevel)","")

all:
	$(MAKE) -C ..

clean:
	$(MAKE) -C .. clean

endif

endif
