LDLIBS=-lm
CFLAGS=-g

all: cpic plot

cpic: cpic.c specie.c mat.c

plot: plot.c
	$(CC) -lGL -lGLU -lglut -lm $< -o $@
