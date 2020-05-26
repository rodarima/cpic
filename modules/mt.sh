#!/bin/sh

module purge

unset PATH
unset LD_LIBRARY_PATH
unset C_INCLUDE_PATH
unset CXX_INCLUDE_PATH

# Sane PATH
export PATH=/bin:/usr/bin

# Load required modules
module load gcc/9.2.0 intel/2018.1 openmpi/3.0.0-cuda cuda/10.1 fftw/3.3.8

# Fix OpenMPI module configuration
export C_INCLUDE_PATH=/apps/OPENMPI/3.0.0/include/:$C_INCLUDE_PATH

BASE=/gpfs/projects/bsc15/bsc15557/mt

# Debug fails to compile in MinoTauro
#export PATH=$BASE/llvm-mono/build/bin:$PATH

# Load clang compiler with OmpSs-2 support
export PATH=$BASE/llvm-mono/build-release/bin:$PATH

# Load heFFTe
export LD_LIBRARY_PATH=$BASE/heffte/usr/lib:$LD_LIBRARY_PATH
export C_INCLUDE_PATH=$BASE/heffte/usr/include:$C_INCLUDE_PATH
export CXX_INCLUDE_PATH=$BASE/heffte/usr/include:$CXX_INCLUDE_PATH

# Load libconfig
export LD_LIBRARY_PATH=$BASE/libconfig/usr/lib:$LD_LIBRARY_PATH
export C_INCLUDE_PATH=$BASE/libconfig/usr/include:$C_INCLUDE_PATH

# Finally set LIBRARY_PATH as well
export LIBRARY_PATH=$LD_LIBRARY_PATH
