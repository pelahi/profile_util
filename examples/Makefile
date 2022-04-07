.PHONY: all clean

OMPFLAGS=-fopenmp
libdir=$(shell pwd)/../lib/

all: openmpvec_cpp openmpvec_c openmpvec_fort

clean:
	rm -f openmpvec_cpp openmpvec_c openmpvec_fort

openmpvec_cpp: openmpvec.cpp
	$(CXX) $(OMPFLAGS) -L$(libdir) -I../src/ -Wl,-rpath=$(libdir) -lprofile_util_omp openmpvec.cpp -o openmpvec_cpp

openmpvec_c: openmpvec.c
	$(CC) $(OMPFLAGS) -L$(libdir) -I../src/ -Wl,-rpath=$(libdir) -lprofile_util_omp openmpvec.c -o openmpvec_c

openmpvec_fort: openmpvec.f90
	$(F90) $(FOMPFLAGS) -L$(libdir) -I../src/ -Wl,-rpath=$(libdir) -lprofile_util_omp openmpvec.f90 -o openmpvec_fort
	
