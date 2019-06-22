#!/bin/bash

#interactive, debug, bsc_cs
#SBATCH --qos=debug
#SBATCH --job-name=cpic.@NCPUS@
#SBATCH --workdir=/home/bsc15/bsc15557/cpic/perf/constant-cpus
#SBATCH --output=out/@NCPUS@
#SBATCH --error=err/@NCPUS@
#SBATCH --ntasks=@NPROCS@
#SBATCH --nodes=1
#SBATCH --cpus-per-task=@NCPUS@
#SBATCH --exclusive
#SBATCH --time=00:30:00

. ../../modules

ulimit -s 8192

#export NANOS6_REPORT_PREFIX="#"

srun --cpu-bind=verbose,cores ../../cpic base.conf
