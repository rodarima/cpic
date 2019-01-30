LDLIBS=-lm
CFLAGS=-g

USE_OMPSS=yes

OMPSS_CC=mcc
OMPSS_CFLAGS=--ompss-2 --instrumentation

ifeq ($(USE_OMPSS),yes)
	CC=$(OMPSS_CC)
	CFLAGS+=$(OMPSS_CFLAGS)
endif

all: cpic plot

cpic: cpic.c specie.c mat.c

plot: plot.c
	$(CC) -lGL -lGLU -lglut -lm $< -o $@
