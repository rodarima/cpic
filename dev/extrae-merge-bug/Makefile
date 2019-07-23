CC=mpicc

CFLAGS=-g -O0

BIN=main worker

all: $(BIN)

run: $(BIN)
	#mpiexec.hydra -n 1 ./main
	mpirun -n 1 --oversubscribe ./main

trace: $(BIN)
	#mpiexec.hydra -n 1 ./trace.sh ./main
	mpirun -n 1 --oversubscribe ./trace.sh ./main

clean:
	rm -f $(BIN) *.o
