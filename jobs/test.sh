#!/bin/bash

#SBATCH --qos=debug
#SBATCH --job-name=cpic
#SBATCH --workdir=.
#SBATCH --output=log/%j.out
#SBATCH --error=log/%j.err
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=24
#SBATCH --tasks-per-node=1
#SBATCH --time=00:02:00

. hard-modules

extrae/trace.sh ./cpic conf/mpi.conf

mpi2prv -f TRACE.mpits -o trace/cpic.prv
