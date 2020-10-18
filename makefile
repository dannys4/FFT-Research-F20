INCDIR=.
CC=g++-8
SIMDARGS= -mavx2
CFLAGS=-I$(INCDIR) $(SIMDARGS) --std=c++14 -ggdb3 -Wall
LINKARGS=-lfmt -pthread
SRCDIR=.
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
