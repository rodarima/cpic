CC=clang
OCC=mcc
LDLIBS=-lm -lconfig -lfftw3 -lgsl -lgslcblas
CFLAGS=-g -I./include/ -I/apps/PM/ompss-2/2018.11/include/ -L /apps/PM/ompss-2/2018.11/lib

USE_OMPSS=yes

CPIC_SRC=specie.c particle.c block.c mat.c block.c sim.c \
	 field.c cpic.c


ifeq ($(USE_OMPSS), no)
#OMPSS_CFLAGS=-k --ompss-2 --instrumentation
OCFLAGS=--ompss-2
LDLIBS+=-lnanos6-optimized
CPIC_SRC:=$(CPIC_SRC:.c=.mcc.c)
CPIC_SRC+=loader.c
endif

CPIC_OBJ=$(CPIC_SRC:.c=.o)

SRC=$(CPIC_SRC)
OBJ=$(SRC:.c=.o)

BIN=cpic eplot pplot plot config fft solver

all: $(BIN)


test: test.mcc.c
	$(CC) $(CFLAGS) $(LDLIBS) $^ -o $@

test2: test2.mcc.c loader.c
	$(CC) $(CFLAGS) $(LDLIBS) $^ -o $@

cpic: $(CPIC_OBJ)

%.mcc.c: %.c
	$(OCC) $(CFLAGS) $(OCFLAGS) -y -o $@ $<

solver: mat.o solver.o

plot: plot.c
	$(CC) $(CFLAGS) $(LDLIBS) -lGL -lGLU -lglut -lm $< -o $@

pplot: pplot.c
	$(CC) $(CFLAGS) -lGL -lGLU -lglut -lm $< -o $@

eplot: eplot.c
	$(CC) $(CFLAGS) -lGL -lGLU -lglut -lm $< -o $@

config: config.o

clean:
	rm -rf *.o *.mcc.c $(BIN)

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

vg: cpic
	valgrind --fair-sched=yes ./cpic

doc: $(SRC)
	doxygen .doxygen

.PHONY: doc

.PRECIOUS: %.mcc.c
