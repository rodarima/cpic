#!/bin/bash

#export EXTRAE_HOME=/usr

source ${EXTRAE_HOME}/etc/extrae.sh

export EXTRAE_CONFIG_FILE=extrae.xml
export LD_PRELOAD=${EXTRAE_HOME}/lib/libmpitrace.so # For C apps
#export LD_PRELOAD=${EXTRAE_HOME}/lib/libnanosmpitrace.so
#export LD_PRELOAD=${EXTRAE_HOME}/lib/libnanostrace.so

#export LD_PRELOAD="/apps/GCC/8.1.0/lib/gcc/x86_64-suse-linux/8.1.0/libasan.so.5 $LD_PRELOAD"
#export LD_PRELOAD="$LD_PRELOAD /usr/lib64/libasan.so.0"


export EXTRAE_ON=1
#export NANOS6_REPORT_PREFIX="@@@"
#export NANOS6=debug

# For extrae
export NANOS6=extrae

#For heavy debug
#export NANOS6=verbose-debug
#export NANOS6_VERBOSE=all
export NANOS6_EXTRAE_AS_THREADS=0

## Run the desired program
#NANOS6=extrae
$*
