CC=mpicc

CFLAGS=-g -O0

BIN=main worker

all: $(BIN)

run: $(BIN)
	mpiexec.hydra -n 1 ./main
	#mpirun --oversubscribe -n 1 ./main

trace: $(BIN)
	mpiexec.hydra -n 1 ./trace.sh ./main
	#mpirun --oversubscribe -n 1 ./trace.sh ./main

clean:
	rm -f $(BIN) *.o
