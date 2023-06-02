.PHONY: all clean

GPU=hip
MPICXX=hipcc
LIBDIR=$(shell pwd)/../../lib/
LIB=-lprofile_util_${GPU}_mpi_omp -lmpi -lmpi_gtl_hsa
INCDIR=$(shell pwd)/../../include/ 
FLAGS=-L$(LIBDIR) -I$(INCDIR) -Wl,-rpath=$(LIBDIR) -D_HIP -D_MPI -std=c++17 -fopenmp 
CXXFLAGS?=-I${CRAY_MPICH_DIR}/include/ -L${CRAY_MPICH_DIR}/lib/ -L${CRAY_MPICH_ROOTDIR}/gtl/lib/

all: bin/gpu-mpi-comm

clean:
	rm -f bin/gpu-mpi-comm


bin/gpu-mpi-comm: gpu-mpi-comm.cpp
	$(MPICXX) $(FLAGS) $(CXXFLAGS) gpu-mpi-comm.cpp -o bin/gpu-mpi-comm $(LIB)

