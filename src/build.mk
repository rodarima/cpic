module:=src

bin:=cpic cpic.a

src:=$(wildcard $(module)/*.c)
src:=$(filter-out $(wildcard $(module)/*.mcc.c),$(src))
src:=$(filter-out $(module)/plot.c,$(src))
src:=$(filter-out $(module)/video.c,$(src))

USE_MCC=0

ifeq ($(USE_MCC), 1)
	src:=$(subst .c,.mcc.c,$(src))
endif


obj:=$(subst .c,.o,$(src))


src_lib:=$(filter-out $(module)/cpic.c,$(src))
obj_lib:=$(subst .c,.o,$(src))

src_ldlibs:=-lm -lconfig -lgsl -lgslcblas
src_cflags:=-g -pthread -Wall

src_ldlibs+=-lfftw3_mpi -lfftw3

# Extrae API (not needed anymore)
#src_ldlibs+=-lmpitrace

#src_ldlibs+=-lGL -lmgl2
#src_cflags+=`pkg-config --cflags glfw3`
#src_ldlibs+=`pkg-config --libs glfw3`

src_cflags+=$(shell mpicc --showme:compile)
src_ldlibs+=$(shell mpicc --showme:link)

%.mcc.c: %.c
	mcc --ompss-2 -y $^ -o $@

%.mcc.o: %.mcc.c
	$(COMPILE.c) $(OUTPUT_OPTION) $<

.PRECIOUS: %.mcc.c


cpic: $(obj)
	$(CC) $(CFLAGS) $(src_cflags) $(LDFLAGS) $(LDLIBS) $(src_ldlibs) $^ -o $@

cpic.a: $(obj_lib)
	ar rcs $@ $^

# Add to main rules
SRC += $(src)
BIN += $(bin)
