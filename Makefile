# Customize compiler and python version here if desired:
CC = gcc
CXX = g++
PY = python3

#
RM = rm -f

#CXXFLAGS+= $(shell $(PY)-config --cflags)
CXXFLAGS+= -Iextern/eigen -Iextern
CXXFLAGS+= -std=c++11 -fopenmp
CXXFLAGS+= -ffast-math -funsafe-math-optimizations -fno-finite-math-only -O3
CXXFLAGS+= -DNDEBUG
CXXFLAGS+= -Wno-unused-variable -Wno-maybe-uninitialized
#CXXFLAGS+= -fPIC -flto -fsized-deallocation

#LDFLAGS+= -shared
#LDFLAGS+= $(shell $(PY)-config --ldflags)

MAIN=feyntrop examples/feyntrop tests/feyntrop

.PHONY: depend clean

all:    $(MAIN)

.depend :

examples/feyntrop : feyntrop
	ln $^ $@

tests/feyntrop : feyntrop
	ln $^ $@

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

