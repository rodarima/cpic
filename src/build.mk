module:=src

#bin:=cpic mft_worker cpic.a
bin:=cpic cpic.a

all_src:=$(wildcard $(module)/*.c)
all_src:=$(filter-out $(wildcard $(module)/*.mcc.c),$(all_src))
all_src:=$(filter-out $(module)/plot.c,$(all_src))
all_src:=$(filter-out $(module)/video.c,$(all_src))
src:=$(filter-out $(module)/mft_worker.c,$(all_src))

USE_MCC?=0
obj:=

ifeq ($(USE_MCC), 1)
src:=$(subst .c,.mcc.c,$(src))
CFLAGS+=-D_MCC -D_MERCURIUM -D_OMPSS_2=1
#CFLAGS+=-D_MCC -D_MERCURIUM -D_OMPSS_2=1 -include nanos6.h
#OBJ+=/usr/lib/nanos6-main-wrapper.o
#obj+=/usr/lib/nanos6-main-wrapper.o
LDFLAGS+=-L/usr/lib -Wl,-z,lazy -Xlinker -rpath -Xlinker /usr/lib
LDLIBS+=-lnanos6 -ldl /usr/lib/nanos6-main-wrapper.o
GEN+=$(src)
endif


obj_cpic=$(subst .c,.o,$(src))
obj_worker=src/mft_worker.o src/tap.o src/utils.o src/mat.o


src_lib:=$(filter-out $(module)/cpic.c,$(src))
obj_lib:=$(subst .c,.o,$(src))

src_cflags:=
src_ldlibs:=

ifeq ($(USE_TAMPI), 1)
# Add TAMPI BEFORE MPI
TAMPI_HOME?=/usr
#src_ldlibs=-ltampi-c
src_ldlibs=$(TAMPI_HOME)/lib/libtampi-c.a
src_ldlibs+=-lmpi_cxx
endif

src_cflags+=$(shell mpicc --showme:compile)
src_ldlibs+=$(shell mpicc --showme:link)
src_ldlibs+=-lmpi

src_ldlibs+=-lm -lconfig
#src_ldlibs+=-lgsl -lgslcblas
src_cflags+=-g -pthread -Wall

#src_ldlibs+=-lfftw3_omp
src_ldlibs+=-lfftw3_threads
src_ldlibs+=-lfftw3_mpi -lfftw3

# HDF5
src_ldlibs+=-lhdf5

# Extrae API (not needed anymore)
#src_ldlibs+=-lmpitrace

#src_ldlibs+=-lGL -lmgl2
#src_cflags+=`pkg-config --cflags glfw3`
#src_ldlibs+=`pkg-config --libs glfw3`

MCC_CFLAGS:=

# Enable OmpSs 2
#MCC_CFLAGS+=--ompss-2

# Show line markers in the generated file
MCC_CFLAGS+=--line-markers

# For intel
MCC_CFLAGS+=--cc=$(CC)
MCC_CFLAGS+=--cxx=$(CPP)
MCC_CFLAGS+=--ld=$(CPP)
MCC_CFLAGS+=--v

%.mcc.c: %.c
	mcc $(MCC_CFLAGS) -y $^ -o $@

%.mcc.o: %.mcc.c
	$(COMPILE.c) $(OUTPUT_OPTION) $^

.PRECIOUS: %.mcc.c


#echo "LD $@"

#cpic: $(obj_cpic)
#	mcxx --ompss-2 --line-markers $(MCC_CFLAGS) $(CFLAGS) $(src_cflags) $^ $(src_ldlibs) $(LDFLAGS) $(LDLIBS) -o $@

cpic: $(obj_cpic)
	$(CXX) $(CFLAGS) $(src_cflags) $^ $(src_ldlibs) $(LDFLAGS) $(LDLIBS) -o $@

#mcxx --ld=clang++ --ompss-2 --line-markers $(CFLAGS) $(src_cflags) $^ $(src_ldlibs) $(LDFLAGS) $(LDLIBS) -o $@

#WORKERS_CFLAGS=-O0 -static-libasan -fsanitize=address -fno-omit-frame-pointer -no-pie
#WORKERS_CFLAGS+=-DGLOBAL_DEBUG
#WORKERS_CFLAGS=

src/mft_worker.o: src/mft_worker.c
	$(CC) $(WORKERS_CFLAGS) -c -o $@ $^

mft_worker: $(obj_worker)
	@echo "LD $@"
	@I_MPI_CC=gcc $(MPICC) $(WORKERS_CFLAGS) -lm -lfftw3_mpi -lfftw3 $(CFLAGS) -o $@ $^

#$(GCC) $(CFLAGS) $(src_cflags) $(WORKERS_CFLAGS) $(src_ldlibs) $(LDFLAGS) $(LDLIBS) -o $@ $^

cpic.a: $(obj_lib)
	ar rcs $@ $^

# Add to main rules
SRC += $(src)
OBJ += $(obj_cpic) $(obj_worker)
BIN += $(bin)
