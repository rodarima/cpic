#!/bin/bash
USE_INTEL=0

module purge

export PATH="/usr/bin:/bin"
export CPATH=""
export LIBRARY_PATH=""
export LD_LIBRARY_PATH=""

if [ "$USE_INTEL" == "1" ]; then
	module load intel impi
	module load gcc/7.2.0 EXTRAE/latest
	module load gcc/8.1.0
else
	module load gcc/8.1.0
	module load openmpi/3.1.1
	module load gcc/7.2.0
	module load EXTRAE/latest
	module load gcc/8.1.0
fi

#export CPATH=$CPATH:$TAMPI_HOME/include:/apps/INTEL/2017.4/impi/2017.3.196/include64
#export LIBRARY_PATH=$LIBRARY_PATH:$TAMPI_HOME/lib:/apps/INTEL/2017.4/impi/2017.3.196/lib64:$HOME/root/usr/lib
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$TAMPI_HOME/lib:$HOME/root/usr/lib
