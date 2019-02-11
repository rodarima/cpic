CC=gcc
OCC=mcc
LDLIBS=-lm 
CFLAGS=-g -I./include/ -I/apps/PM/ompss-2/2018.11/include/ -L /apps/PM/ompss-2/2018.11/lib

USE_OMPSS=no

ifeq ($(USE_OMPSS), yes)
#OMPSS_CFLAGS=-k --ompss-2 --instrumentation
CFLAGS+=--ompss-2
LDLIBS+=-lnanos6-extrae
LOADER+=loader.o
endif

all: cpic plot test

test: test.mcc.c
	$(CC) $(CFLAGS) $(LDLIBS) $^ -o $@

cpic: cpic.mcc.o $(LOADER) specie.mcc.o mat.mcc.o block.mcc.o
	$(CC) $(CFLAGS) $(LDLIBS) $^ -o $@

%.mcc.c: %.c
	$(OCC) $(CFLAGS) $(OCFLAGS) -y -o $@ $<

plot: plot.c
	$(CC) -lGL -lGLU -lglut -lm $< -o $@

clean:
	rm -rf *.o *.mcc.c cpic

load:
	module load gcc/7.2.0 extrae ompss-2

run:
	./cpic
	mpi2prv -f TRACE.mpits -o trace/output.prv
runmn:
	LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/apps/PM/ompss-2/2018.11/lib taskset -c 0-20 ./cpic
	#NANOS6=extrae taskset -c 0-25 ./cpic
	#NANOS6=extrae taskset -c 0-20 ./cpic
	${EXTRAE_HOME}/bin/mpi2prv -f TRACE.mpits -o output.prv

.PRECIOUS: %.mcc.c
