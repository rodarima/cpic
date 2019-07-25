#!/bin/bash

echo "It's working!"

#. modules.sh
#env

#ompi_info --all > log.info

export LD_PRELOAD=""
valgrind --log-file=${PMI_RANK}.vg ./mft_worker
./mft_worker
