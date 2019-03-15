CC=gcc
LDLIBS=-lm -lstdc++ -lGL -lGLU -lglut -lmgl2-glut -lmgl2 -lplplot

BIN=mgl

all: $(BIN) plplot

mgl: mgl-glue.o

clean:
	rm -rf *.o $(BIN)
