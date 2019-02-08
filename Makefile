CC=gcc
LDLIBS=-lm
CFLAGS=-g -Iinclude/

USE_OMPSS=no

OMPSS_CC=mcc
OMPSS_CFLAGS=-k --ompss-2
#OMPSS_CFLAGS=-k --ompss-2 --instrumentation

ifeq ($(USE_OMPSS),yes)
	CC=$(OMPSS_CC)
	CFLAGS+=$(OMPSS_CFLAGS)
endif

all: cpic plot

cpic: cpic.c specie.o mat.o block.o

plot: plot.c
	$(CC) -lGL -lGLU -lglut -lm $< -o $@

clean:
	rm -rf *.o cpic

load:
	module load gcc/7.2.0 extrae ompss-2

run:
	NANOS6=extrae taskset -c 0-25 ./cpic
	#NANOS6=extrae taskset -c 0-20 ./cpic
	${EXTRAE_HOME}/bin/mpi2prv -f TRACE.mpits -o output.prv
