CFLAGS := -Wall -O3 -fopenmp -ffast-math -march=native -I.

all: count

debug: CFLAGS := -Wall -g -fopenmp
debug: all

count: count.c bott.o
	gcc -o count count.c $(CFLAGS) bott.o

bott.o: bott.c bott.h
	gcc -c -o bott.o bott.c $(CFLAGS)

clean:
	@rm count bott.o
