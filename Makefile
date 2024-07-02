CFLAGS := -Wall -O3 -fopenmp -ffast-math -march=native

all: count

debug: CFLAGS := -Wall -g -fopenmp
debug: all

count: count.c
	gcc -o count count.c $(CFLAGS)

clean:
	@rm count
