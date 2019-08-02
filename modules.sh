#!/bin/bash

USE_INTEL=1
USE_OPENMPI_4=0
USE_PATCHED_OPENMPI=0
USE_NEW_SCHEDULER=1

module purge

export PATH="/usr/bin:/bin"
export CPATH=""
export LIBRARY_PATH=""
export LD_LIBRARY_PATH=""

if [ $USE_INTEL == 0 ]; then

	module load intel mkl gsl
	echo "AFTER INTEL:"
	echo "LIBRARY_PATH=$LIBRARY_PATH"
	echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH"
	module load gcc/7.2.0 openmpi/3.1.1 extrae ompss-2 tampi
	#module load gcc/7.2.0 openmpi/3.1.1 extrae/3.7.0 ompss-2 tampi
	export CPATH=$CPATH:$TAMPI_HOME/include:/apps/OPENMPI/3.1.1/GCC/include:$HOME/root/usr/include:$HOME/cpic/include
	export LIBRARY_PATH=$LIBRARY_PATH:/apps/PM/TAMPI/1.0/openmpi/3.1.1/lib:/apps/OPENMPI/3.1.1/GCC/lib:$HOME/root/usr/lib
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/apps/PM/TAMPI/1.0/openmpi/3.1.1/lib:/apps/OPENMPI/3.1.1/GCC/lib:$HOME/root/usr/lib

	if [ $USE_OPENMPI_4 == 1 ]; then
		module unload openmpi/3.1.1
		module load gcc/8.1.0
		#openmpi/4.0.1

		if [ $USE_PATCHED_OPENMPI == 0 ]; then

			module load openmpi/4.0.1

		else
			# Load my patched openmpi/4.0.1
			export CPATH=$CPATH:$HOME/apps/openmpi-4.0.1/include
			export PATH=$HOME/apps/openmpi-4.0.1/bin:$PATH
			export LIBRARY_PATH=$HOME/root/usr/lib:$HOME/apps/openmpi-4.0.1/lib:$LIBRARY_PATH
			export LD_LIBRARY_PATH=$LIBRARY_PATH

			echo "----------- WARNING WARNING WARNING WARNING ----------"
			echo " MY FIXED OPENMPI/4.0.1 HAS BEEN LOADED, BE CAREFUL   "
			echo "----------- WARNING WARNING WARNING WARNING ----------"
		fi
	fi
else

	module load gcc/7.2.0 EXTRAE/latest
	module load intel mkl gsl impi fftw ompss-2 tampi
	module load gcc/8.1.0
	export CPATH=$CPATH:$TAMPI_HOME/include:/apps/INTEL/2017.4/impi/2017.3.196/include64:$HOME/cpic/include:$HOME/root/usr/include
	export LIBRARY_PATH=$LIBRARY_PATH:$TAMPI_HOME/lib:/apps/INTEL/2017.4/impi/2017.3.196/lib64:$HOME/root/usr/lib
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$TAMPI_HOME/lib:$HOME/root/usr/lib

fi


if [ $USE_NEW_SCHEDULER == 1 ]; then
	echo LOADING NEW SCHEDULER!
	export LD_LIBRARY_PATH=$HOME/apps/nanos6-scheduling-refactor/lib:$LD_LIBRARY_PATH
fi


module load hdf5
