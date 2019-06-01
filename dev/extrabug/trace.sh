#!/bin/bash

#export EXTRAE_HOME=/usr

source ${EXTRAE_HOME}/etc/extrae.sh

export EXTRAE_CONFIG_FILE=extrae.xml
export LD_PRELOAD=${EXTRAE_HOME}/lib/libmpitrace.so # For C apps

export EXTRAE_ON=1
export NANOS6_REPORT_PREFIX="@@@"
export NANOS6=extrae
#export NANOS6=debug
#export NANOS6=vaerbose-debug

## Run the desired program
#NANOS6=extrae
$*
