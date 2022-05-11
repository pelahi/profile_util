.PHONY: all clean

F90=gfortran 
OMPFLAGS=-fopenmp -DUSEOPENMP
FOMPFLAGS=-fopenmp 
libdir=$(shell pwd)/../../lib/
incdir=$(shell pwd)/../../src/

all: bin/openmpvec_cpp bin/openmpvec_c bin/openmpvec_fort

clean:
	rm -f bin/openmpvec_cpp bin/openmpvec_c bin/openmpvec_fort

bin/openmpvec_cpp: openmpvec.cpp
	$(CXX) $(OMPFLAGS) -L$(libdir) -I$(incdir) -Wl,-rpath=$(libdir) openmpvec.cpp -o bin/openmpvec_cpp -lprofile_util_omp 

bin/openmpvec_c: openmpvec.c
	$(CC) $(OMPFLAGS) -L$(libdir) -I$(incdir) -Wl,-rpath=$(libdir) openmpvec.c -o bin/openmpvec_c -lprofile_util_omp 

bin/openmpvec_fort: openmpvec.f90
	$(F90) $(FOMPFLAGS) -L$(libdir) -I$(incdir) -Wl,-rpath=$(libdir) openmpvec.f90 -o bin/openmpvec_fort -lprofile_util_omp

	
