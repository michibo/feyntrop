# Customize compiler version here if desired:
CC = gcc
CXX = g++
# This might be necessary on macOS where the default C++ compiler does not support OpenMP.

#
RM = rm -f
LN = ln -f

# On macOS it might help to remove the following line. Doing so unfortunately, makes feyntrop much slower. A better way is to install a C++ compiler that supports OpenMP.
CXXFLAGS+= -fopenmp
# Configuration of C++ compiler
CXXFLAGS+= -Iextern/eigen -Iextern
CXXFLAGS+= -std=c++11 
CXXFLAGS+= -ffast-math -funsafe-math-optimizations -fno-finite-math-only -O3
#CXXFLAGS+= -DNDEBUG
CXXFLAGS+= -Wno-unused-variable -Wno-maybe-uninitialized

MAIN=feyntrop examples/feyntrop tests/feyntrop

.PHONY: depend clean

all:    $(MAIN)

.depend :

examples/feyntrop : feyntrop
	$(LN) $^ $@

tests/feyntrop : feyntrop
	$(LN) $^ $@

feyntrop : feyntrop.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

feyntrop.o : srcs/feyntrop.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


depend: .depend

.depend : srcs/feyntrop.cpp
	$(RM) ./.depend
	$(CXX) $(CXXFLAGS) -MM $^>>./.depend;

clean:
	$(RM) $(MAIN) feyntrop.o .depend

include .depend

