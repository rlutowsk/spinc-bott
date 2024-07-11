CFLAGS := -Wall -Ofast -fopenmp -march=native -I.

SRC := count.c backtrack.c
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

profile: CFLAGS := -Wall -pg -fopenmp -Ofast -fno-inline
profile: $(APP)

$(APP): $(SRC) $(OBJ)
	$(CC) -o $@ $@.c $(CFLAGS) $(OBJ)

%.o: %.c %.h
	$(CC) -o $@ $< -c $(CFLAGS)

%.lo: %.c %.h
	$(GAC) -o $@ $< -c -p "$(CFLAGS)"

%.so: %.c $(LOBJ)
	$(GAC) -d $< -o $@ -p "$(CFLAGS)" -L "$(LOBJ)"

clean:
	@rm $(wildcard $(TARGETS) $(OBJ) $(LOBJ))
