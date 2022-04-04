.PHONY: all clean

OMPFLAGS=-fopenmp
libdir=$(shell pwd)/../lib/

all: openmpvec

clean:
	rm -f openmpvec

openmpvec: openmpvec.cpp
	$(CXX) $(OMPFLAGS) -L$(libdir) -I../src/ -Wl,-rpath=$(libdir) -lprofile_util_omp openmpvec.cpp -o openmpvec
