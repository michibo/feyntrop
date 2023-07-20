# Customize compiler version here if desired:
CC = gcc
CXX = g++

#
RM = rm -f
LN = ln -f

CXXFLAGS+= -Iextern/eigen -Iextern
CXXFLAGS+= -std=c++11 -fopenmp
CXXFLAGS+= -ffast-math -funsafe-math-optimizations -fno-finite-math-only -O3
CXXFLAGS+= -DNDEBUG
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

