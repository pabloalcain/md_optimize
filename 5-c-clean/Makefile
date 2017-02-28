.PHONY: default help object executable all clean
CC = gcc

CC_FLAGS = -g -std=gnu99 -O3 -pg
LD_FLAGS = -lm -pg

LD = $(CC)

SOURCE_C = $(wildcard *.c)
OBJECTS_C = $(patsubst %.c, %_c.o, $(SOURCE_C))

EXECUTABLE = md_c.e

default: all

objects: $(OBJECTS_C)

executable: $(EXECUTABLE)

all: objects executable

%_c.o: %.c
	$(CC) $(CC_FLAGS) -c $^ -o $@

md_c.e: $(OBJECTS_C)
	$(LD) $^ $(LD_FLAGS) -o $@

clean:
	rm -rfv $(OBJECTS_C) $(EXECUTABLE)
