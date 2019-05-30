#!/bin/bash

export EXTRAE_HOME=/usr

source /usr/etc/extrae.sh

export EXTRAE_CONFIG_FILE=extrae/extrae2.xml
export LD_PRELOAD=${EXTRAE_HOME}/lib/libmpitrace.so # For C apps

## Run the desired program
$*
