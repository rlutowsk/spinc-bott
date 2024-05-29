CFLAGS := -Wall -O3 -fopenmp -ffast-math -march=native

all: bott

debug: CFLAGS := -Wall -g -fopenmp
debug: bott

bott: bott.c
	gcc -o bott bott.c $(CFLAGS)

clean:
	@rm bott
