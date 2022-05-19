.PHONY: all clean

LIBDIR=$(shell pwd)/../../lib/
LIB=-lprofile_util_mpi 
INCDIR=$(shell pwd)/../../include/
FLAGS=-L$(LIBDIR) -I$(INCDIR) -Wl,-rpath=$(LIBDIR) -D_MPI

all: bin/mpivec_cpp 

clean:
	rm -f bin/mpivec_cpp

bin/mpivec_cpp: mpivec.cpp
	$(MPICXX) $(FLAGS) mpivec.cpp -o bin/mpivec_cpp $(LIB)

	
