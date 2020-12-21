INCDIR=./include
CC=g++-9
SIMDARGS= -mfma -mavx2
STD=11
CFLAGS=-I$(INCDIR) $(SIMDARGS) -ggdb3 -Wall -std=c++$(STD)
SRCDIR=./src
OBJDIR=$(SRCDIR)/obj

LIBS=

# This will only work on a Unix-like environment
dir_guard=@mkdir -p $(@D)

_DEPS=$(wildcard $(INCDIR)/*.hpp) 
DEPS=$(patsubst $(INCDIR)/%,%,$(_DEPS))

_OBJ=$(patsubst $(SRCDIR)/%.cpp, $(OBJDIR)/%.o, $(wildcard $(SRCDIR)/*.cpp)) 
OBJ=$(patsubst %,%,$(_OBJ))

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(INCDIR)/%.hpp
	$(dir_guard)
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

run: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

.PHONY: clean all default

clean:
	rm -f $(OBJDIR)/*.o *~ $(INCDIR)/*~ run
