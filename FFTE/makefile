INCDIR=./include
CC=g++
SIMDARGS= -mfma -mavx2
CFLAGS=-I$(INCDIR) $(SIMDARGS) -ggdb3 -Wall -O3
SRCDIR=./src
OBJDIR=$(SRCDIR)/obj

LIBS=

_DEPS=$(wildcard $(INCDIR)/*.hpp) 
DEPS=$(patsubst $(INCDIR)/%,%,$(_DEPS))

_OBJ=$(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(wildcard $(SRCDIR)/*.cpp)) 
OBJ=$(patsubst %,%,$(_OBJ))

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

run: $(OBJ)
	$(CC) $(LINKARGS) $(CFLAGS) -o $@ $^ $(LIBS)

.PHONY: clean all default

clean:
	rm -f $(OBJDIR)/*.o *~ core $(INCDIR)/*~ run
