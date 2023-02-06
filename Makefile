CC = gcc
CXX = g++
RM = rm -f

CXXFLAGS= 
CXXFLAGS+= `python3-config --cflags`
CXXFLAGS+= -Iextern/eigen -Iextern/pybind11/include
CXXFLAGS+= -std=c++14 -fopenmp
CXXFLAGS+= -ffast-math -funsafe-math-optimizations -fno-finite-math-only -O3
CXXFLAGS+= -DNDEBUG 
CXXFLAGS+= -Wno-unused-variable -Wno-maybe-uninitialized
CXXFLAGS+= -fPIC -flto -fsized-deallocation

LDFLAGS+= -shared `python3-config --ldflags`

MAIN=feyntrop.so

.PHONY: depend clean

all:    $(MAIN)

.depend :

feyntrop.so : feyntrop.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

feyntrop.o : srcs/feyntrop.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


depend: .depend

.depend: $(SRCS)
	$(RM) ./.depend
	$(CXX) $(CXXFLAGS) -MM $^>>./.depend;

clean:
	$(RM) $(MAIN) feyntrop.o .depend 

include .depend

