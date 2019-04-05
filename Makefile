CC=clang
OCC=mcc
LDLIBS=-lm -lconfig -lfftw3 -lgsl -lgslcblas -lGL -lGLU -lglut -lmgl2
CFLAGS=-g -pthread -Wall
LDFLAGS=-pthread

CFLAGS_GLFW3=`pkg-config --cflags glfw3`
LDLIBS_GLFW3=`pkg-config --libs glfw3`

CFLAGS+=$(CFLAGS_GLFW3)
LDLIBS+=$(LDLIBS_GLFW3)


USE_OMPSS=yes

CPIC_SRC=specie.c particle.c block.c mat.c block.c sim.c \
	 field.c cpic.c solver.c config.c plot.c interpolate.c


ifeq ($(USE_OMPSS), no)
#OMPSS_CFLAGS=-k --ompss-2 --instrumentation
OCFLAGS=--ompss-2
CFLAGS+=-I./include/ -I/apps/PM/ompss-2/2018.11/include/ -L /apps/PM/ompss-2/2018.11/lib
LDLIBS+=-lnanos6-optimized
CPIC_SRC:=$(CPIC_SRC:.c=.mcc.c)
CPIC_SRC+=loader.c
endif

CPIC_OBJ=$(CPIC_SRC:.c=.o)

SRC=$(CPIC_SRC)
OBJ=$(SRC:.c=.o)
DEP=$(SRC:.c=.d)

BIN=cpic interpolate.test



all: $(BIN)

include $(DEP)

# rule to generate a dep file by using the C preprocessor
# (see man cpp for details on the -MM and -MT options)
%.d: %.c
	@$(CPP) $(CFLAGS) $< -MM -MT $(@:.d=.o) >$@

test: test.mcc.c
	$(CC) $(CFLAGS) $(LDLIBS) $^ -o $@

test2: test2.mcc.c loader.c
	$(CC) $(CFLAGS) $(LDLIBS) $^ -o $@

interpolate.test: interpolate.test.c interpolate.o mat.o


cpic: $(CPIC_OBJ)

%.mcc.c: %.c
	$(OCC) $(CFLAGS) $(OCFLAGS) -y -o $@ $<

clean:
	rm -f $(OBJ) $(BIN) $(DEP)

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
