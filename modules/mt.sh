#!/bin/sh

module purge

unset PATH
unset LD_LIBRARY_PATH
unset C_INCLUDE_PATH
unset CXX_INCLUDE_PATH

# Sane PATH
export PATH=/bin:/usr/bin:/usr/sbin:/sbin:/opt/perf/bin

# Load required modules
#module load intel/2018.1 impi/2018.1
module load gcc/9.2.0 intel/2018.1 openmpi/3.0.0-cuda
module load cuda/10.1 fftw/3.3.8 ompss-2/git

# Load our installed software
export PATH=$HOME/mt/usr/bin:$PATH
export LIBRARY_PATH=$HOME/mt/usr/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=$HOME/mt/usr/lib:$LD_LIBRARY_PATH
export CPATH=$HOME/mt/usr/include:$CPATH

# Not used actually
export OMPI_CC=clang
export OMPI_CXX=clang++

# Fix OpenMPI module configuration
export C_INCLUDE_PATH=/apps/OPENMPI/3.0.0/include/:$C_INCLUDE_PATH

# Fix missing lib/ and include/ in OmpSs-2 module
export C_INCLUDE_PATH=/apps/PM/ompss-2/git/include:$C_INCLUDE_PATH
export CXX_INCLUDE_PATH=/apps/PM/ompss-2/git/include:$CXX_INCLUDE_PATH
export LD_LIBRARY_PATH=/apps/PM/ompss-2/git/lib:$LD_LIBRARY_PATH
export NANOS6_HOME=/apps/PM/ompss-2/git/


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
