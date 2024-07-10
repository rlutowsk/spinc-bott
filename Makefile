CFLAGS := -Wall -O3 -fopenmp -ffast-math -march=native -I.

SRC := count.c backtrack.c backtrack1.c
APP := $(patsubst %.c, %, $(SRC)) 
LIB := gap.so
TARGETS := $(APP) $(LIB)
OBJ  := bott.o
LOBJ := bott.lo

CC = gcc
GAC = gac

all: $(TARGETS)

debug: CFLAGS := -Wall -g -fopenmp -O0
debug: all

$(APP): $(SRC) $(OBJ)
	$(CC) -o $@ $@.c $(CFLAGS) $(OBJ)

#count: count.c $(OBJ)
#	gcc -o count count.c $(CFLAGS) bott.o
#
#backtrack: backtrack.c $(OBJ)
#	gcc -o backtrack backtrack.c $(CFLAGS) bott.o

%.o: %.c %.h
	$(CC) -o $@ $< -c $(CFLAGS)

#bott.o: bott.c bott.h
#	gcc -c -o bott.o bott.c $(CFLAGS)

#bott.lo: bott.c
#	gac -p "$(CFLAGS)" -c bott.c -o bott.lo

%.lo: %.c %.h
	$(GAC) -o $@ $< -c -p "$(CFLAGS)"

%.so: %.c $(LOBJ)
	$(GAC) -d $< -o $@ -p "$(CFLAGS)" -L "$(LOBJ)"

clean:
	@rm $(wildcard $(TARGETS) $(OBJ) $(LOBJ))
