MODULES:=src test
#MODULES+=test user


#CC=gcc
#CC=clang
CC=mcc --ompss-2 --line-markers
#OCC=mcc
LDLIBS:=
CFLAGS:=-g -Wall
#CFLAGS:=-g -Wall -Werror
#LDFLAGS:=-L. -Wl,-rpath,.

# Enable debug in Mercurium
#CC+=--v

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

#Include all modules for headers
CFLAGS+=$(patsubst %,-I%,$(MODULES))

BIN:=
SRC:=
GEN:=
OBJ:=

HOSTNAME=$(shell hostname)

#NANOS6_HEADER=NANOS6=verbose NANOS6_VERBOSE=all
#NANOS6_HEADER=NANOS6=verbose

# MPI configuration based on the computing device
ifeq ($(HOSTNAME), mio)
 NPROCS?=2
 NCORES?=2
 MPIRUN=mpirun -n $(NPROCS) --map-by NUMA:PE=$(NCORES) --oversubscribe
else
 NPROCS?=4
 NCORES?=4
 MPIRUN=mpirun -n $(NPROCS) --map-by NUMA:PE=$(NCORES)
endif


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
OBJ += \
$(patsubst %.c,%.o, $(filter %.c,$(SRC))) \
$(patsubst %.y,%.o, $(filter %.y,$(SRC)))
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
	@$(CPP) $(CFLAGS) $< -MM -MT $(@:.d=.o) >$@


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
	$(MPIRUN) bash -c 'valgrind ./cpic conf/mpi.conf 2> log/$$PMIX_RANK.log'

trace: trace/cpic.prv

run:
	rm -f log/*
	$(MPIRUN) bash -c '$(NANOS6_HEADER) ./cpic conf/mpi.conf 2> log/$$PMIX_RANK.log'

gprof:
	GMON_OUT_PREFIX=gmon taskset -c 0-15 mpirun --oversubscribe -n 16 bash -c './cpic conf/mpi.conf 2> log/$$PMIX_RANK.log'

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
