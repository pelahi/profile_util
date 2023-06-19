.PHONY: all clean

CXX=hipcc
OMPFLAGS=-DUSEOPENMP -fopenmp --offload-arch=gfx90a -fopenmp-targets=amdgcn-amd-amdhsa -Xopenmp-target=amdgcn-amd-amdhsa -march=gfx90a  
GPUFLAGS=-D_HIP -D__HIP_PLATFORM_AMD__ 
LIBDIR=$(shell pwd)/../../lib/
LIB=-lprofile_util_hip_omp
INCDIR=$(shell pwd)/../../include/
FLAGS=-std=c++14 -fPIC -I$(INCDIR) 
LFLAGS=-L$(LIBDIR)  -Wl,-rpath=$(LIBDIR)

all: bin/gpu-openmp

clean:
	rm -f bin/gpu-openmp

bin/gpu-openmp: src/gpu-openmp.cpp src/kernels.h src/kernels.cu src/common.h
	$(CXX) $(GPUFLAGS) $(FLAGS) -c src/kernels.cu -o obj/kernels.o
	$(CXX) $(GPUFLAGS) $(FLAGS) $(OMPFLAGS) -c src/gpu-openmp.cpp -o obj/gpu-openmp.o
	$(CXX) $(GPUFLAGS) $(FLAGS) $(OMPFLAGS) obj/*.o -o bin/gpu-openmp $(LFLAGS) $(LIB)

