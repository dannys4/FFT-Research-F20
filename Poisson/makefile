INCDIR=/home/dsharp/Documents/FFT-Research-F20
CC=g++
SIMDARGS= -mfma -mavx2 -fopenmp
STD=11
CFLAGS=-I$(INCDIR) $(SIMDARGS) -Ofast -ffast-math -Wall -std=c++$(STD) -D_P_IDX=$(P)
SRCDIR=./src
OBJDIR=$(SRCDIR)/obj

LIBS=

# This will only work on a Unix-like environment
dir_guard=@mkdir -p $(@D)

_DEPS=$(wildcard $(INCDIR)/*.hpp) 
DEPS=$(patsubst $(INCDIR)/%,%,$(_DEPS))

_OBJ=$(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(wildcard $(SRCDIR)/*.cpp)) 
OBJ=$(patsubst %,%,$(_OBJ))

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(dir_guard)
	$(CC) $(CFLAGS) -c $< -o $@

run: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

.PHONY: clean all default

clean:
	rm -f $(OBJDIR)/*.o *~ $(INCDIR)/*~ run
