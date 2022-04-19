.PHONY: all clean

libdir=$(shell pwd)/../../lib/
incdir=$(shell pwd)/../../src/

all: bin/mpivec_cpp 

clean:
	rm -f bin/mpivec_cpp

bin/mpivec_cpp: mpivec.cpp
	$(MPICXX) -L$(libdir) -I$(incdir) -Wl,-rpath=$(libdir) -lprofile_util_mpi mpivec.cpp -o bin/mpivec_cpp

	
