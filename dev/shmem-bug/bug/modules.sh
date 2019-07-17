#!/bin/bash
USE_INTEL=1

module purge

export PATH="/usr/bin:/bin"
export CPATH=""
export LIBRARY_PATH=""
export LD_LIBRARY_PATH=""

module load gcc/7.2.0 EXTRAE/latest
module load intel impi
module load gcc/8.1.0

#export CPATH=$CPATH:$TAMPI_HOME/include:/apps/INTEL/2017.4/impi/2017.3.196/include64
#export LIBRARY_PATH=$LIBRARY_PATH:$TAMPI_HOME/lib:/apps/INTEL/2017.4/impi/2017.3.196/lib64:$HOME/root/usr/lib
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$TAMPI_HOME/lib:$HOME/root/usr/lib
