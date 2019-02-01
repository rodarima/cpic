CC=clang
LDLIBS=-lm
CFLAGS=-g

USE_OMPSS=no

OMPSS_CC=mcc
OMPSS_CFLAGS=--ompss-2 --instrumentation

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
