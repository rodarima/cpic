CC=gcc
CFLAGS_GLFW3=`pkg-config --cflags glfw3`
LDLIBS_GLFW3=`pkg-config --libs glfw3`

LDLIBS=-lm -lstdc++ -lGL -lmgl2

CFLAGS+=$(CFLAGS_GLFW3)
LDLIBS+=$(LDLIBS_GLFW3)

BIN=mgl glvid

all: $(BIN)

clean:
	rm -rf *.o $(BIN)
