CC=OMPI_CC=mcc mpicc
#CC=I_MPI_CC=mcc mpiicc

CFLAGS=-g --ompss-2

#CFLAGS+=$(shell mpicc --showme:compile)
#LDFLAGS+=$(shell mpicc --showme:link)

CFLAGS+=-fsanitize=address -fno-omit-frame-pointer

LDLIBS+=-lm
LDLIBS+=-ltampi-c
#LDLIBS+=-lompitrace

TASKS=4
CPUS_PER_TASK=4

all: test

test: perf.o

run: test
	#mpirun -n 2 --map-by NUMA:PE=2 ./test
	srun -w $$SLURM_NODELIST --cpu_bind=cores --ntasks=${TASKS} --cpus-per-task=${CPUS_PER_TASK} ./test

clean:
	rm -f test *.o
