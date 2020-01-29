CFLAGS=

VECTOR_BIT_WIDTH=256
#VECTOR_BIT_WIDTH=512

CFLAGS+=-DUSE_VECTOR_$(VECTOR_BIT_WIDTH)

# For DEBUG
#CFLAGS+=-fsanitize=address
#CFLAGS+=-fno-omit-frame-pointer
#CFLAGS+=-O0

# Common
#CFLAGS+=-O0
#CFLAGS+=-O2
CFLAGS+=-O3
CFLAGS+=-g -Wall
#CFLAGS+=-masm=intel

CFLAGS+=-std=c11

CFLAGS+=-I../../src/
#CFLAGS+=-fno-omit-frame-pointer

# Default compiler
CC=clang


## GCC
GCC=gcc
GCFLAGS=$(CFLAGS)
#GCFLAGS+=-fopenmp
GCFLAGS+=-fopenmp-simd
#GCFLAGS+=-fopt-info-vec-all
GCFLAGS+=-march=native
GCFLAGS+=-mtune=native
#GCFLAGS+=-march=skylake-avx512
#GCFLAGS+=-mtune=skylake-avx512
#GCFLAGS+=-mprefer-vector-width=512
GCFLAGS+=-fverbose-asm
GCFLAGS+=-fdump-tree-optimized
GCFLAGS+=-ftree-vectorize
#GCFLAGS+=-ftree-vectorizer-verbose=1
#GCFLAGS+=-fopt-info-vec-all
#GCFLAGS+=-fopt-info-vec-all-internals
#GCFLAGS+=-fira-verbose=50
#GCFLAGS+=-fsched-pressure

# LLVM/Clang
LCC=clang
LCFLAGS=$(CFLAGS)
#LCFLAGS+=--cc=clang --ompss-2 --v
#LCFLAGS+=-march=native
LCFLAGS+=-march=core-avx2
#LCFLAGS+=-march=skylake-avx512

# Intel
ICC=icc
ICFLAGS=$(CFLAGS)
##CFLAGS+=-axCORE-AVX512
ICFLAGS+=-fsource-asm
ICFLAGS+=-xHost
#ICFLAGS+=-xCOMMON-AVX512
#ICFLAGS+=-xCOMMON-AVX512
#ICFLAGS+=-qopt-report=5
#ICFLAGS+=-xCORE-AVX512 -qopt-report=5
ICFLAGS+=-qopt-zmm-usage=high
#ICFLAGS+=-xCOMMON-AVX512 -qopt-report=5
#ICFLAGS+=-fopenmp-simd

MCC=mcc
MCFLAGS=$(CFLAGS)
MCFLAGS+=--cc=clang --cpp=clang
MCFLAGS+=--ompss-2 -k --v
#MCFLAGS+=-fsanitize=address -fno-omit-frame-pointer

LDLIBS+=-L/apps/PAPI/5.6.0/lib -lm
#LDLIBS+=-lpapi

COMPILERS=$(MCC)
#COMPILERS=$(LCC)
#COMPILERS=$(LCC) $(GCC)
#COMPILERS=$(ICC) $(GCC) $(LCC)
BIN=$(addprefix test.,$(COMPILERS)) #my-kernel
ASM=$(addsuffix .s, $(BIN) $(OBJ))
OBJ=$(addprefix kernel.,$(COMPILERS)) #my-kernel
VINST=$(addsuffix .v, $(OBJ))

#.SECONDARY:
.PRECIOUS: $(ASM)

all: $(BIN) #$(ASM) $(VINST)

%.icc.s: %.c
	$(ICC) $(filter-out -g,$(ICFLAGS)) -S $^ -o $@

%.icc.o: %.icc.s
	$(ICC) -c $(ICFLAGS) $< -o $@

test.icc: test.c perf.c mem.icc.o kernel.icc.o ../../src/utils.o
	$(ICC) $(ICFLAGS) $^ $(LDLIBS) -o $@

################################################################

%.clang.s: %.c
	$(LCC) $(filter-out -g,$(LCFLAGS)) -S $^ -o $@

#%.clang.o: %.c
%.clang.o: %.clang.s
	$(LCC) -c $(LCFLAGS) $< -o $@

test.clang: test.c perf.c mem.o kernel.clang.o ../../src/utils.o
	$(LCC) $(LCFLAGS) $^ $(LDLIBS) -o $@

################################################################

%.gcc.s: %.c
	$(GCC) $(filter-out -g,$(GCFLAGS)) -S $^ -o $@

%.gcc.o: %.gcc.s
	$(GCC) -c $(GCFLAGS) $< -o $@

test.gcc: test.c perf.c mem.o kernel.gcc.o ../../src/utils.o
	$(GCC) $(GCFLAGS) $(LDLIBS) $^ -o $@

################################################################

test.mcc: test.c perf.c mem.o kernel.clang.s ../../src/utils.o
	$(MCC) $(MCFLAGS) $^ $(LDLIBS) -o $@

################################################################

my-kernel.o: my-kernel.s
	$(GCC) -c $(GCFLAGS) $< -o $@

my-test: test.c perf.c my-kernel.o
	$(GCC) $(GCFLAGS) $(LDLIBS) $^ -o $@

################################################################

%.v: %.s
	grep '^\s*v' $^ | sort | awk '{print $$1}' | uniq -c > $@

################################################################

#icc-intrinsic: intrinsic.c
#	$(ICC) $(ICFLAGS) $^ $(LDLIBS) -o $@
#
#gcc-intrinsic: intrinsic.c
#	$(GCC) $(GCFLAGS) $^ $(LDLIBS) -o $@

clean:
	rm -f *.o $(BIN) $(ASM) $(VINST)


perf: test.clang
	sudo taskset -c 1 perf stat -a -C 1 --topdown -I 100 ./test.clang

sync:
	rsync -av test.c kernel.c test.clang mn:
