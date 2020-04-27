#!/bin/sh
set -x
mpirun --bind-to core -n  1 ./fft
mpirun --bind-to core -n  2 ./fft
mpirun --bind-to core -n  4 ./fft
mpirun --bind-to core -n  8 ./fft
mpirun --bind-to core -n 16 ./fft
