#MODULES:=src
MODULES:=src test
#MODULES+=test user


#GCC=gcc
#CC=clang
#CC=mcc --ompss-2 --line-markers
#I_MPI_CC=mcc
#CC=I_MPI_CC=mcc mpiicc --ompss-2 --line-markers
MPICC=mpicc
#CC=OMPI_CC=mcc $(MPICC) -k --ompss-2
#CC=OMPI_CC=mcc $(MPICC) --line-markers

#CC=I_MPI_CC=mcc mpiicc --cc=clang --ompss-2 --line-markers

#CC=OMPI_CC=mcc $(MPICC) --ompss-2 --line-markers

CC=clang
CPP=clang

#CC=gcc
#CPP=gcc

#CC=icc
#CPP=icpc

#MCC=mcc --ompss-2 --cc=$(CC) --pp=$(CC) --v --cpp=$(CC)

#CFLAGS+=-qopt-zmm-usage=high

#OCC=mcc
LDLIBS:=

# Add libtasio
#LDLIBS+=-ltasio

CFLAGS+=-g -Wall -Wextra -Wpedantic
CFLAGS+=-Wsign-conversion
#CFLAGS+=-Wpadded
CFLAGS+=-Wmissing-prototypes
CFLAGS+=-Werror

CFLAGS+=-std=c11

# Use the new clang with ompss2 support
#CFLAGS+=-fompss-2

#Disable ompss pragma warnings if we are not using OmpSs-2
CFLAGS+=-Wno-unknown-pragmas

#Extra warnings
#CFLAGS+=-Wstrict-prototypes -Wshadow -Wconversion

# Optimization enabled
#CFLAGS+=-O3
#CFLAGS+=-ffp-contract=fast
CFLAGS+=-O0
CFLAGS+=-g
#CFLAGS+=-Wall -Werror
#LDFLAGS:=-L. -Wl,-rpath,.

# Enable debug in Mercurium
#CC+=--v

# Enable profiling with gprof
#CFLAGS+=-pg

# Ignore pragma warnings
#CFLAGS+=-Wno-unknown-pragmas

# Use debug messages
CFLAGS+=-DGLOBAL_DEBUG

# No asserts
#CFLAGS+=-DNDEBUG

# No extra assers
#CFLAGS+=-DNO_EXTRA_ASSERTS

# For intel compiler
CFLAGS+=-fPIE

# Use 256 bits for vector operations
CFLAGS+=-DUSE_VECTOR_256
CFLAGS+=-march=core-avx2
#CFLAGS+=-xHost

# Use TAMPI
USE_TAMPI?=0

ifeq ($(USE_TAMPI), 1)
CFLAGS+=-DWITH_TAMPI
endif

# Instrument functions so Extrae can get some information
#CFLAGS+=-finstrument-functions

# Debug
#CFLAGS+=-fsanitize=address

# For perf
#CFLAGS+=-ggdb -g3
CFLAGS+=-fno-omit-frame-pointer

# Debug race conditions
#CFLAGS+=-fsanitize=thread

# Avoid optimized instructions, so we can still use Valgrind
#CFLAGS+=-march=x86-64 -mtune=generic

# Stack protector
CFLAGS+=-fstack-protector-all

# Optimization
#CFLAGS+=-O2

#Include all modules for headers
CFLAGS+=$(patsubst %,-I%,$(MODULES))

BIN:=
SRC:=
GEN:=
OBJ:=

HOSTNAME=$(shell hostname)

#NANOS6_HEADER=NANOS6=verbose NANOS6_VERBOSE=all
#NANOS6_HEADER=NANOS6=verbose

#NANOS6_DEBUG=NANOS6=verbose NANOS6_VERBOSE=all
NANOS6_DEBUG=NANOS6=verbose

#LSAN_OPTIONS=detect_leaks=0

# MPI configuration based on the computing device
ifeq ($(HOSTNAME), mio)
 NPROCS?=2
 NCORES?=2
 MPIRUN=mpirun -n $(NPROCS) --map-by NUMA:PE=$(NCORES) --oversubscribe
 ENV_RANK=PMIX_RANK
else
 NNODES?=1
 PROCS_PER_NODE?=1
 CPUS_PER_TASK?=32
# MPIRUN=mpirun --bind-to core -n $(PROCS_PER_NODE) --map-by NUMA:PE=$(CPUS_PER_TASK)
 #MPIRUN=mpirun -n $(PROCS_PER_NODE) --map-by node:pe=$(CPUS_PER_TASK)
 #ENV_RANK=PMIX_RANK

 # Doesn't start the simulator main()
 MPIRUN=srun -N $(NNODES) --cpu-bind=verbose,cores --cpus-per-task $(CPUS_PER_TASK) --tasks-per-node $(PROCS_PER_NODE)
 ENV_RANK=PMI_RANK
endif

#CPIC_CONF=perf/constant-cpus/base.conf
CPIC_CONF=conf/mpi.conf


all:


# include the description for each module
INC_MAKEFILES:=$(patsubst %,%/build.mk,$(MODULES))

include $(INC_MAKEFILES)

# include the C include dependencies
#$(info Including ${OBJ:.o=.d})
#include $(OBJ:.o=.d)
#$(info $$OBJ is [${OBJ}])
DEP:=$(patsubst %.c,%.d, $(filter %.c,$(SRC)))

include $(DEP)

# determine the object files
#OBJ += \
#$(patsubst %.c,%.o, $(filter %.c,$(SRC))) \
#$(patsubst %.y,%.o, $(filter %.y,$(SRC)))
#
#$(info $$SRC is [${SRC}])
#$(info $$OBJ is [${OBJ}])
#OBJ=$(SRC:.c=.o)
#DEP:=$(SRC:.c=.d)
#$(info Including ${OBJ})

#include $(DEP)

# rule to generate a dep file by using the C preprocessor
# (see man cpp for details on the -MM and -MT options)

#%.mcc.d: %.mcc.c
#	$(CPP) $(CFLAGS) $< -MM -MT $(@:.d=.o) >$@

%.d: %.c
	@$(CC) $(CFLAGS) $< -MM -MT $(@:.d=.o) >$@


#%.o: %.c
#	@echo CC $@
#	$(COMPILE.c) $(OUTPUT_OPTION) $^


#@echo "CC $<"
#	$(MCC) $(CFLAGS) -c -o $@ $<
%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

all: $(BIN)

#test/cyclotron: $(CPIC_OBJ)

#nm -g cpic | sed -e 's/ . / /g' -e '/GLIBC/d' -e '/ _/d' -e '/^ /d' -e 's/ /#/g' -e '/interpol/d' > function.list

extrae/function.all: cpic
	nm -g cpic > $@

extrae/function.list: extrae/function.all extrae/filter.list
	grep -wf extrae/filter.list $< | sed -e 's/ . /#/g' > $@

clean:
	rm -f $(OBJ) $(BIN) $(DEP) $(GEN)

#%.mcc.c: %.c
#	$(OCC) $(CFLAGS) $(OCFLAGS) -y -o $@ $<
#
#load:
#	module load gcc/7.2.0 extrae ompss-2
#
trace/cpic.prv: cpic extrae/extrae2.xml extrae/trace.sh extrae/function.list conf/mpi.conf Makefile
	rm -rf set-0/ TRACE.sym TRACE.mpits
	$(MPIRUN) extrae/trace.sh ./cpic conf/mpi.conf
	#mpirun -n 2 --cpus-per-proc 4 extrae/trace.sh ./cpic conf/mpi.conf
	#mpirun -n 2 taskset -c 0-15 extrae/trace.sh ./cpic conf/mpi.conf
	mpi2prv -f TRACE.mpits -o trace/cpic.prv
	grep 'User function' trace/cpic.pcf

valgrind:
	$(MPIRUN) bash -c 'valgrind ./cpic conf/mpi.conf 2> log/$$$(ENV_RANK).log'

trace: trace/cpic.prv

run: cpic
	ulimit -s $$((8*1024))
	rm -f log/*
	$(MPIRUN) bash -c '$(NANOS6_HEADER) ./cpic $(CPIC_CONF) 2> log/$$$(ENV_RANK).log'

debug: cpic
	rm -f log/*
	$(MPIRUN) bash -c '$(NANOS6_DEBUG) ./cpic $(CPIC_CONF) 2> log/$$$(ENV_RANK).log'

gdb: cpic
	#$(MPIRUN) xterm -e bash -c 'gdb --args ./cpic conf/mpi.conf 2> log/$$$(ENV_RANK).log'
	$(MPIRUN) xterm -e gdb --args ./cpic conf/mpi.conf

gprof:
	GMON_OUT_PREFIX=gmon taskset -c 0-15 mpirun --oversubscribe -n 16 bash -c './cpic conf/mpi.conf 2> log/$$$(ENV_RANK).log'

.PHONY: run valgrind trace gprof

#runmn:
#	LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/apps/PM/ompss-2/2018.11/lib taskset -c 0-20 ./cpic
#	#NANOS6=extrae taskset -c 0-25 ./cpic
#	#NANOS6=extrae taskset -c 0-20 ./cpic
#	${EXTRAE_HOME}/bin/mpi2prv -f TRACE.mpits -o output.prv
#
#vg: cpic
#	valgrind --fair-sched=yes ./cpic
#
#.PRECIOUS: %.mcc.c
