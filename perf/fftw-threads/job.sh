#!/bin/bash

#interactive, debug, bsc_cs
#SBATCH --qos=debug
#SBATCH --job-name=cpic.@NT@
#SBATCH --workdir=/home/bsc15/bsc15557/cpic/perf/fftw-threads
#SBATCH --output=out/@NT@
#SBATCH --error=err/@NT@
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --exclusive
#SBATCH --time=00:30:00

. ../../modules

ulimit -s 8192

#export NANOS6_REPORT_PREFIX="#"

srun --cpu-bind=verbose,cores ../../cpic conf/@NT@
