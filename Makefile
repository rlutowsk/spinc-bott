CFLAGS := -Wall -O3 -fopenmp -ffast-math -march=native -I.

all: count gap.so

debug: CFLAGS := -Wall -g -fopenmp
debug: all

count: count.c bott.o
	gcc -o count count.c $(CFLAGS) bott.o

bott.o: bott.c bott.h
	gcc -c -o bott.o bott.c $(CFLAGS)

bott.lo: bott.c
	gac -p "$(CFLAGS)" -c bott.c -o bott.lo

gap.so: gap.c bott.lo
	gac -d gap.c -o gap.so -p "$(CFLAGS)" -L "bott.lo"

clean:
	@rm count bott.o gap.so
