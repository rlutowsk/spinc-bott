CFLAGS := -Wall -O3 -fopenmp -ffast-math -march=native -I.

TARGETS := count backtrack gap.so
OBJ  := bott.o
LOBJ := bott.lo

all: $(TARGETS)

debug: CFLAGS := -Wall -g -fopenmp
debug: all

count: count.c $(OBJ)
	gcc -o count count.c $(CFLAGS) bott.o

backtrack: backtrack.c $(OBJ)
	gcc -o backtrack backtrack.c $(CFLAGS) bott.o

bott.o: bott.c bott.h
	gcc -c -o bott.o bott.c $(CFLAGS)

bott.lo: bott.c
	gac -p "$(CFLAGS)" -c bott.c -o bott.lo

gap.so: gap.c $(LOBJ)
	gac -d gap.c -o gap.so -p "$(CFLAGS)" -L "bott.lo"

clean:
	@rm $(wildcard $(TARGETS) $(OBJ) $(LOBJ))
