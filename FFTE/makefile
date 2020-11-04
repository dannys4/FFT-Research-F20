INCDIR=./include
CC=g++-7
SIMDARGS= -mfma -mavx2
CFLAGS=-I$(INCDIR) $(SIMDARGS) -ggdb3 -Wall
SRCDIR=./src
OBJDIR=$(SRCDIR)/obj

LIBS=

_DEPS=$(wildcard $(INCDIR)/*.hpp) 
DEPS=$(patsubst $(INCDIR)/%,%,$(_DEPS))

_OBJ=$(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(wildcard $(SRCDIR)/*.cpp)) 
OBJ=$(patsubst %,%,$(_OBJ))

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(INCDIR)/%.hpp
	$(CC) $(CFLAGS) -c $< -o $@

run: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

.PHONY: clean all default

clean:
	rm -f $(OBJDIR)/*.o *~ $(INCDIR)/*~ run
