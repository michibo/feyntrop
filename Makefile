# Customize compiler and python version here if desired:
CC = gcc
CXX = g++
PY = python3

#
RM = rm -f

CXXFLAGS+= $(shell $(PY)-config --cflags)
CXXFLAGS+= -Iextern/eigen -Iextern/pybind11/include
CXXFLAGS+= -std=c++14 -fopenmp
CXXFLAGS+= -ffast-math -funsafe-math-optimizations -fno-finite-math-only -O3
CXXFLAGS+= -DNDEBUG
CXXFLAGS+= -Wno-unused-variable -Wno-maybe-uninitialized
CXXFLAGS+= -fPIC -flto -fsized-deallocation

LDFLAGS+= -shared
LDFLAGS+= $(shell $(PY)-config --ldflags)

MAIN=feyntrop.so

.PHONY: depend clean

all:    $(MAIN)

.depend :

feyntrop.so : feyntrop.o
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

