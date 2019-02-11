CC=gcc
LDLIBS=-lm -lnanos6-extrae
CFLAGS=-g -Iinclude/

USE_OMPSS=no

OCC=mcc
OCFLAGS=--ompss-2
#OMPSS_CFLAGS=-k --ompss-2 --instrumentation

ifeq ($(USE_OMPSS),yes)
	CC=$(OMPSS_CC)
	CFLAGS+=$(OMPSS_CFLAGS)
endif

all: cpic plot test

test: test.mcc.c
	$(CC) $(CFLAGS) $(LDLIBS) $^ -o $@

cpic: cpic.mcc.o loader.o specie.mcc.o mat.mcc.o block.mcc.o
	$(CC) $(CFLAGS) $(LDLIBS) $^ -o $@

%.mcc.c: %.c
	$(OCC) $(OCFLAGS) -y -o $@ $<

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
	NANOS6=extrae taskset -c 0-25 ./cpic
	#NANOS6=extrae taskset -c 0-20 ./cpic
	${EXTRAE_HOME}/bin/mpi2prv -f TRACE.mpits -o output.prv

#.PRECIOUS: %.mcc.c
