module:=src

bin:=cpic cpic.a

src:=$(wildcard $(module)/*.c)
src:=$(filter-out $(wildcard $(module)/*.mcc.c),$(src))
src:=$(filter-out $(module)/plot.c,$(src))
src:=$(filter-out $(module)/video.c,$(src))

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


obj+=$(subst .c,.o,$(src))


src_lib:=$(filter-out $(module)/cpic.c,$(src))
obj_lib:=$(subst .c,.o,$(src))

src_cflags:=
src_ldlibs:=

ifeq ($(USE_TAMPI), 1)
# Add TAMPI BEFORE MPI
TAMPI_HOME?=/usr
#src_ldlibs=-ltampi-c
src_ldlibs=$(TAMPI_HOME)/lib/libtampi.a
endif

#src_cflags+=$(shell mpicc --showme:compile)
#src_ldlibs+=$(shell mpicc --showme:link)
src_ldlibs+=-lmpi

src_ldlibs+=-lm -lconfig -lgsl -lgslcblas
src_cflags+=-g -pthread -Wall

#src_ldlibs+=-lfftw3_omp
src_ldlibs+=-lfftw3_threads
src_ldlibs+=-lfftw3_mpi -lfftw3

# Extrae API (not needed anymore)
#src_ldlibs+=-lmpitrace

#src_ldlibs+=-lGL -lmgl2
#src_cflags+=`pkg-config --cflags glfw3`
#src_ldlibs+=`pkg-config --libs glfw3`

MCC_CFLAGS:=

# Enable OmpSs 2
MCC_CFLAGS+=--ompss-2

# Show line markers in the generated file
MCC_CFLAGS:=--line-markers

%.mcc.c: %.c
	mcc $(MCC_CFLAGS) -y $^ -o $@

%.mcc.o: %.mcc.c
	$(COMPILE.c) $(OUTPUT_OPTION) $<

.PRECIOUS: %.mcc.c


cpic: $(obj)
	mcxx --ompss-2 --line-markers $(CFLAGS) $(src_cflags) $^ $(src_ldlibs) $(LDFLAGS) $(LDLIBS) -o $@

cpic.a: $(obj_lib)
	ar rcs $@ $^

# Add to main rules
SRC += $(src)
BIN += $(bin)
