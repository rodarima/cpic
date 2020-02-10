#!/bin/bash

CC=clang
CPU=1
VER=aslr0

N=20

DIR="log/${VER}/"

mkdir -p $DIR

for i in $(seq $N $((2*$N))); do
	echo "CC $CC, CPU $CPU, $i"
	sudo taskset -c $CPU ./test.$CC > ${DIR}/$i
done
