#!/bin/sh
set -x
mpirun -n 1 ./example >> log.txt
mpirun -n 2 ./example >> log.txt
mpirun -n 4 ./example >> log.txt
mpirun -n 8 ./example >> log.txt
mpirun -n 16 ./example >> log.txt
