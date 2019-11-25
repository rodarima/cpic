CFLAGS=

# Common
CFLAGS+=-Ofast
#CFLAGS+=-O2
CFLAGS+=-g -Wall
#CFLAGS+=-masm=intel

## GCC
GCC=gcc
GCFLAGS=$(CFLAGS)
#GCFLAGS+=-fopenmp
GCFLAGS+=-fopenmp-simd
#GCFLAGS+=-fopt-info-vec-all
#GCFLAGS+=-march=native
#GCFLAGS+=-mtune=native
GCFLAGS+=-march=skylake-avx512
GCFLAGS+=-mtune=skylake-avx512
GCFLAGS+=-mprefer-vector-width=512
GCFLAGS+=-fverbose-asm
GCFLAGS+=-fdump-tree-optimized
GCFLAGS+=-ftree-vectorize
#GCFLAGS+=-ftree-vectorizer-verbose=1
#GCFLAGS+=-fopt-info-vec-all
#GCFLAGS+=-fopt-info-vec-all-internals
#GCFLAGS+=-fira-verbose=50
#GCFLAGS+=-fsched-pressure

# Intel
ICC=icc
ICFLAGS=$(CFLAGS)
##CFLAGS+=-axCORE-AVX512
ICFLAGS+=-fsource-asm
ICFLAGS+=-xCOMMON-AVX512
#ICFLAGS+=-qopt-report=5
#ICFLAGS+=-xCORE-AVX512 -qopt-report=5
ICFLAGS+=-qopt-zmm-usage=high
#ICFLAGS+=-xCOMMON-AVX512 -qopt-report=5
#ICFLAGS+=-fopenmp-simd

#CC=mcc
#CFLAGS+=--ompss-2 -k --v
#CFLAGS+=-fsanitize=address -fno-omit-frame-pointer

LDLIBS+=-lm

BIN=icc-test gcc-test my-test
ASM=icc-kernel.s gcc-kernel.s

all: $(BIN) $(ASM)

icc-kernel.s: kernel.c
	$(ICC) $(filter-out -g,$(ICFLAGS)) -S $^ -o $@

icc-kernel.o: icc-kernel.s
	$(ICC) -c $(ICFLAGS) $< -o $@

icc-test: test.c perf.c icc-kernel.o
	$(ICC) $(ICFLAGS) $^ $(LDLIBS) -o $@

################################################################

gcc-kernel.s: kernel.c
	$(GCC) $(filter-out -g,$(GCFLAGS)) -S $^ -o $@

gcc-kernel.o: gcc-kernel.s
	$(GCC) -c $(GCFLAGS) $< -o $@

gcc-test: test.c perf.c gcc-kernel.o
	$(GCC) $(GCFLAGS) $(LDLIBS) $^ -o $@

my-kernel.o: my-kernel.s
	$(GCC) -c $(GCFLAGS) $< -o $@

my-test: test.c perf.c my-kernel.o
	$(GCC) $(GCFLAGS) $(LDLIBS) $^ -o $@


################################################################

#icc-intrinsic: intrinsic.c
#	$(ICC) $(ICFLAGS) $^ $(LDLIBS) -o $@
#
#gcc-intrinsic: intrinsic.c
#	$(GCC) $(GCFLAGS) $^ $(LDLIBS) -o $@

clean:
	rm -f *.o $(BIN) $(ASM)

