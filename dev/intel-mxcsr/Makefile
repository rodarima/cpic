BIN=main.gcc main.icc main.clang main.mcc

all: $(BIN)

main.gcc: main.c
	gcc -Ofast -ffast-math $^ -o $@

main.icc: main.c
	icc -Ofast -fp-model fast $^ -o $@

main.clang: main.c
	clang -Ofast -ffast-math $^ -o $@

main.mcc: main.c
	mcc -ffast-math $^ -o $@

run: $(BIN)
	@for EXE in $(BIN) ; do \
		echo -en "$$EXE\\t";\
		./$$EXE ; \
	done

clean:
	rm -f $(BIN)
