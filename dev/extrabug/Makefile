#CC=gcc
#CC=clang
CC=mcc --v --ompss-2 --line-markers --cpp=mpicc --ld=mpicc
#OCC=mcc
LDLIBS:=
CFLAGS:=-g -Wall
#CFLAGS:=-g -Wall -Werror
#LDFLAGS:=-L. -Wl,-rpath,.

# Enable profiling with gprof
#CFLAGS+=-pg

# Use debug messages
#CFLAGS+=-DGLOBAL_DEBUG

# Instrument functions so Extrae can get some information
CFLAGS+=-finstrument-functions

# Avoid optimized instructions, so we can still use Valgrind
#CFLAGS+=-march=x86-64 -mtune=generic

# Stack protector
#CFLAGS+=-fstack-protector-all

# Optimization
#CFLAGS+=-O2

all: extrabug

function.all: extrabug
	nm -g $^ > $@

function.list: function.all filter.list
	grep -wf filter.list $< | sed -e 's/ . /#/g' > $@


trace: extrabug extrae.xml trace.sh function.list Makefile
	rm -rf set-0/ TRACE.sym TRACE.mpits
	mpirun -n 2 --map-by NUMA:PE=4 trace.sh ./extrabug
	#mpirun -n 2 ./trace.sh ./extrabug
	#./trace.sh ./extrabug
	mpi2prv -f TRACE.mpits -o bug.prv
	grep 'User function' bug.pcf

.PHONY: trace

clean:
	rm -f extrabug bug.*
