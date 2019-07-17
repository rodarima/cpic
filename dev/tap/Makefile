#CC=OMPI_CC=mcc mpicc
#MCC=I_MPI_CC="mcc -k" mpiicc
MPICC=mpiicc
MCC=$(MPICC)
#OMPI_CC=clang
CC=$(MPICC)

CFLAGS=-g -O0
#CFLAGS+=--ompss-2

#CFLAGS+=$(shell mpicc --showme:compile)
#LDFLAGS+=$(shell mpicc --showme:link)

#CFLAGS+=-fsanitize=address -fno-omit-frame-pointer

#LDLIBS+=-lm
#LDLIBS+=-ltampi-c
#LDLIBS+=-lompitrace

TASKS=4
CPUS_PER_TASK=4

BIN=main worker

all: $(BIN)

worker: tap.o

main: main.c tap.o
	$(MCC) $(CFLAGS) -o $@ $^

#run: test
	#mpirun -n 2 --map-by NUMA:PE=2 ./test
#	srun -w $$SLURM_NODELIST --cpu_bind=cores --ntasks=${TASKS} --cpus-per-task=${CPUS_PER_TASK} ./test

clean:
	rm -f $(BIN) *.o
