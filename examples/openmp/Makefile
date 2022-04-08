.PHONY: all clean

F90=gfortran 
OMPFLAGS=-fopenmp -DUSEOPENMP
FOMPFLAGS=-fopenmp 
libdir=$(shell pwd)/../../lib/

all: bin/openmpvec_cpp bin/openmpvec_c bin/openmpvec_fort

clean:
	rm -f bin/openmpvec_cpp bin/openmpvec_c bin/openmpvec_fort

bin/openmpvec_cpp: openmpvec.cpp
	$(CXX) $(OMPFLAGS) -L$(libdir) -I../src/ -Wl,-rpath=$(libdir) -lprofile_util_omp openmpvec.cpp -o bin/openmpvec_cpp

bin/openmpvec_c: openmpvec.c
	$(CC) $(OMPFLAGS) -L$(libdir) -I../src/ -Wl,-rpath=$(libdir) -lprofile_util_omp openmpvec.c -o bin/openmpvec_c

bin/openmpvec_fort: openmpvec.f90
	$(F90) $(FOMPFLAGS) -L$(libdir) -I../src/ -Wl,-rpath=$(libdir) -lprofile_util_omp openmpvec.f90 -o bin/openmpvec_fort
	
