.PHONY: all clean serial openmp mpi mpi_openmp 
.PHONY: cuda_set cuda cuda_mpi cuda_openmp cuda_mpi_openmp
.PHONY: hip_set hip hip_mpi hip_openmp hip_mpi_openmp

OUTPUTFILE=libprofile_util
OMPFLAGS=-fopenmp
MPICXX?=$(CXX)

#all: lib/$(OUTPUTFILE).so lib/$(OUTPUTFILE)_omp.so lib/$(OUTPUTFILE)_mpi.so lib/$(OUTPUTFILE)_mpi_omp.so
all: non_gpu_set cuda_set hip_set
#lib/$(OUTPUTFILE).a

clean:
	rm -f obj/*o
	rm -f lib/*

non_gpu_set: serial openmp mpi mpi_openmp

serial: lib/$(OUTPUTFILE).so
	echo "Making serial library"
openmp: lib/$(OUTPUTFILE)_omp.so
	echo "Making OpenMP library"
mpi: lib/$(OUTPUTFILE)_mpi.so
	echo "Making MPI library"
mpi_openmp: lib/$(OUTPUTFILE)_mpi_omp.so
	echo "Making MPI+OpenMP library"

cuda_set: cuda cuda_mpi cuda_openmp cuda_mpi_openmp 

cuda: lib/$(OUTPUTFILE)_cuda.so
	echo "Making CUDA"
cuda_mpi: lib/$(OUTPUTFILE)_cuda_mpi.so
	echo "Making CUDA MPI library"
cuda_openmp: lib/$(OUTPUTFILE)_cuda_omp.so
	echo "Making CUDA MPI library"
cuda_mpi_openmp: lib/$(OUTPUTFILE)_cuda_mpi_omp.so
	echo "Making CUDA MPI + OpenMP library"

hip_set: hip hip_mpi hip_openmp hip_mpi_openmp 

hip: lib/$(OUTPUTFILE)_hip.so
	echo "Making HIP"
hip_mpi: lib/$(OUTPUTFILE)_hip_mpi.so
	echo "Making HIP MPI library"
hip_openmp: lib/$(OUTPUTFILE)_hip_omp.so
	echo "Making HIP MPI library"
hip_mpi_openmp: lib/$(OUTPUTFILE)_hip_mpi_omp.so
	echo "Making HIP MPI + OpenMP library"

lib/$(OUTPUTFILE).so: obj/mem_util.o obj/time_util.o obj/thread_affinity_util.o
	mkdir -p lib; 
	$(CXX) -shared obj/mem_util.o obj/time_util.o obj/thread_affinity_util.o -o lib/$(OUTPUTFILE).so

obj/mem_util.o: src/profile_util.h src/mem_util.cpp
	mkdir -p obj;
	$(CXX) $(CXXFLAGS) -fPIC -std=c++14 -c src/mem_util.cpp -o obj/mem_util.o

obj/time_util.o: src/profile_util.h src/time_util.cpp
	mkdir -p obj;
	$(CXX) $(CXXFLAGS)  -fPIC -std=c++14 -c src/time_util.cpp -o obj/time_util.o

obj/thread_affinity_util.o: src/profile_util.h src/thread_affinity_util.cpp
	mkdir -p obj;
	$(CXX) $(CXXFLAGS) -fPIC -std=c++14 -c src/thread_affinity_util.cpp -o obj/thread_affinity_util.o


lib/$(OUTPUTFILE)_omp.so: obj/mem_util_omp.o obj/time_util_omp.o obj/thread_affinity_util_omp.o
	mkdir -p lib; 
	$(CXX) $(OMPFLAGS) -shared obj/mem_util_omp.o obj/time_util_omp.o obj/thread_affinity_util_omp.o -o lib/$(OUTPUTFILE)_omp.so

obj/mem_util_omp.o: src/profile_util.h src/mem_util.cpp
	mkdir -p obj;
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) -fPIC -std=c++14 -c src/mem_util.cpp -o obj/mem_util_omp.o

obj/time_util_omp.o: src/profile_util.h src/time_util.cpp
	mkdir -p obj;
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) -fPIC -std=c++14 -c src/time_util.cpp -o obj/time_util_omp.o

obj/thread_affinity_util_omp.o: src/profile_util.h src/thread_affinity_util.cpp
	mkdir -p obj;
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) -fPIC -std=c++14 -c src/thread_affinity_util.cpp -o obj/thread_affinity_util_omp.o


lib/$(OUTPUTFILE)_mpi.so: obj/mem_util_mpi.o obj/time_util_mpi.o obj/thread_affinity_util_mpi.o
	mkdir -p lib; 
	$(MPICXX) $(MPIFLAGS) -shared obj/mem_util_mpi.o obj/time_util_mpi.o obj/thread_affinity_util_mpi.o -o lib/$(OUTPUTFILE)_mpi.so

obj/mem_util_mpi.o: src/profile_util.h src/mem_util.cpp
	mkdir -p obj;
	$(MPICXX) $(MPIFLAGS) $(CXXFLAGS) -fPIC -std=c++14 -c src/mem_util.cpp -o obj/mem_util_mpi.o

obj/time_util_mpi.o: src/profile_util.h src/time_util.cpp
	mkdir -p obj;
	$(MPICXX) $(MPIFLAGS) $(CXXFLAGS)  -fPIC -std=c++14 -c src/time_util.cpp -o obj/time_util_mpi.o

obj/thread_affinity_util_mpi.o: src/profile_util.h src/thread_affinity_util.cpp
	mkdir -p obj;
	$(MPICXX) $(MPIFLAGS) $(CXXFLAGS) -fPIC -std=c++14 -c src/thread_affinity_util.cpp -o obj/thread_affinity_util_mpi.o

lib/$(OUTPUTFILE)_mpi_omp.so: obj/mem_util_mpi_omp.o obj/time_util_mpi_omp.o obj/thread_affinity_util_mpi_omp.o
	mkdir -p lib; 
	$(MPICXX) $(MPIFLAGS) $(OMPFLAGS) -shared obj/mem_util_mpi_omp.o obj/time_util_mpi_omp.o obj/thread_affinity_util_mpi_omp.o -o lib/$(OUTPUTFILE)_mpi_omp.so

obj/mem_util_mpi_omp.o: src/profile_util.h src/mem_util.cpp
	mkdir -p obj;
	$(MPICXX) $(CXXFLAGS) $(MPIFLAGS) $(OMPFLAGS) -fPIC -std=c++14 -c src/mem_util.cpp -o obj/mem_util_mpi_omp.o

obj/time_util_mpi_omp.o: src/profile_util.h src/time_util.cpp
	mkdir -p obj;
	$(MPICXX) $(CXXFLAGS) $(MPIFLAGS) $(OMPFLAGS) -fPIC -std=c++14 -c src/time_util.cpp -o obj/time_util_mpi_omp.o

obj/thread_affinity_util_mpi_omp.o: src/profile_util.h src/thread_affinity_util.cpp
	mkdir -p obj;
	$(MPICXX) $(CXXFLAGS) $(MPIFLAGS) $(OMPFLAGS) -fPIC -std=c++14 -c src/thread_affinity_util.cpp -o obj/thread_affinity_util_mpi_omp.o

# CUDA 
lib/$(OUTPUTFILE)_cuda.so: obj/mem_util.o obj/time_util.o obj/thread_affinity_util.o
	mkdir -p lib; 
	$(NVCXX) -shared obj/mem_util.o obj/time_util.o obj/thread_affinity_util.o -o lib/$(OUTPUTFILE)_cuda.so

obj/mem_util.o: src/profile_util.h src/mem_util.cpp
	mkdir -p obj;
	$(NVCXX) $(CXXFLAGS) -fPIC -std=c++14 -c src/mem_util.cpp -o obj/mem_util.o

obj/time_util.o: src/profile_util.h src/time_util.cpp
	mkdir -p obj;
	$(NVCXX) $(CXXFLAGS)  -fPIC -std=c++14 -c src/time_util.cpp -o obj/time_util.o

obj/thread_affinity_util.o: src/profile_util.h src/thread_affinity_util.cpp
	mkdir -p obj;
	$(NVCXX) $(CXXFLAGS) -fPIC -std=c++14 -c src/thread_affinity_util.cpp -o obj/thread_affinity_util.o


lib/$(OUTPUTFILE)_cuda_omp.so: obj/mem_util_omp.o obj/time_util_omp.o obj/thread_affinity_util_omp.o
	mkdir -p lib; 
	$(NVCXX) $(OMPFLAGS) -shared obj/mem_util_omp.o obj/time_util_omp.o obj/thread_affinity_util_omp.o -o lib/$(OUTPUTFILE)_cuda_omp.so

obj/mem_util_omp.o: src/profile_util.h src/mem_util.cpp
	mkdir -p obj;
	$(NVCXX) $(CXXFLAGS) $(OMPFLAGS) -fPIC -std=c++14 -c src/mem_util.cpp -o obj/mem_util_omp.o

obj/time_util_omp.o: src/profile_util.h src/time_util.cpp
	mkdir -p obj;
	$(NVCXX) $(CXXFLAGS) $(OMPFLAGS) -fPIC -std=c++14 -c src/time_util.cpp -o obj/time_util_omp.o

obj/thread_affinity_util_omp.o: src/profile_util.h src/thread_affinity_util.cpp
	mkdir -p obj;
	$(NVCXX) $(CXXFLAGS) $(OMPFLAGS) -fPIC -std=c++14 -c src/thread_affinity_util.cpp -o obj/thread_affinity_util_omp.o


lib/$(OUTPUTFILE)_cuda_mpi.so: obj/mem_util_mpi.o obj/time_util_mpi.o obj/thread_affinity_util_mpi.o
	mkdir -p lib; 
	$(MPINVCXX) $(MPIFLAGS) -shared obj/mem_util_mpi.o obj/time_util_mpi.o obj/thread_affinity_util_mpi.o -o lib/$(OUTPUTFILE)_cuda_mpi.so

obj/mem_util_mpi.o: src/profile_util.h src/mem_util.cpp
	mkdir -p obj;
	$(MPINVCXX) $(MPIFLAGS) $(CXXFLAGS) -fPIC -std=c++14 -c src/mem_util.cpp -o obj/mem_util_mpi.o

obj/time_util_mpi.o: src/profile_util.h src/time_util.cpp
	mkdir -p obj;
	$(MPINVCXX) $(MPIFLAGS) $(CXXFLAGS)  -fPIC -std=c++14 -c src/time_util.cpp -o obj/time_util_mpi.o

obj/thread_affinity_util_mpi.o: src/profile_util.h src/thread_affinity_util.cpp
	mkdir -p obj;
	$(MPINVCXX) $(MPIFLAGS) $(CXXFLAGS) -fPIC -std=c++14 -c src/thread_affinity_util.cpp -o obj/thread_affinity_util_mpi.o

lib/$(OUTPUTFILE)_cuda_mpi_omp.so: obj/mem_util_mpi_omp.o obj/time_util_mpi_omp.o obj/thread_affinity_util_mpi_omp.o
	mkdir -p lib; 
	$(MPINVCXX) $(MPIFLAGS) $(OMPFLAGS) -shared obj/mem_util_mpi_omp.o obj/time_util_mpi_omp.o obj/thread_affinity_util_mpi_omp.o -o lib/$(OUTPUTFILE)_cuda_mpi_omp.so

obj/mem_util_mpi_omp.o: src/profile_util.h src/mem_util.cpp
	mkdir -p obj;
	$(MPINVCXX) $(CXXFLAGS) $(MPIFLAGS) $(OMPFLAGS) -fPIC -std=c++14 -c src/mem_util.cpp -o obj/mem_util_mpi_omp.o

obj/time_util_mpi_omp.o: src/profile_util.h src/time_util.cpp
	mkdir -p obj;
	$(MPINVCXX) $(CXXFLAGS) $(MPIFLAGS) $(OMPFLAGS) -fPIC -std=c++14 -c src/time_util.cpp -o obj/time_util_mpi_omp.o

obj/thread_affinity_util_mpi_omp.o: src/profile_util.h src/thread_affinity_util.cpp
	mkdir -p obj;
	$(MPINVCXX) $(CXXFLAGS) $(MPIFLAGS) $(OMPFLAGS) -fPIC -std=c++14 -c src/thread_affinity_util.cpp -o obj/thread_affinity_util_mpi_omp.o

# HIP 
lib/$(OUTPUTFILE)_hip.so: obj/mem_util.o obj/time_util.o obj/thread_affinity_util.o
	mkdir -p lib; 
	$(HIPCXX) -shared obj/mem_util.o obj/time_util.o obj/thread_affinity_util.o -o lib/$(OUTPUTFILE)_hip.so

obj/mem_util.o: src/profile_util.h src/mem_util.cpp
	mkdir -p obj;
	$(HIPCXX) $(CXXFLAGS) -fPIC -std=c++14 -c src/mem_util.cpp -o obj/mem_util.o

obj/time_util.o: src/profile_util.h src/time_util.cpp
	mkdir -p obj;
	$(HIPCXX) $(CXXFLAGS)  -fPIC -std=c++14 -c src/time_util.cpp -o obj/time_util.o

obj/thread_affinity_util.o: src/profile_util.h src/thread_affinity_util.cpp
	mkdir -p obj;
	$(HIPCXX) $(CXXFLAGS) -fPIC -std=c++14 -c src/thread_affinity_util.cpp -o obj/thread_affinity_util.o


lib/$(OUTPUTFILE)_hip_omp.so: obj/mem_util_omp.o obj/time_util_omp.o obj/thread_affinity_util_omp.o
	mkdir -p lib; 
	$(HIPCXX) $(OMPFLAGS) -shared obj/mem_util_omp.o obj/time_util_omp.o obj/thread_affinity_util_omp.o -o lib/$(OUTPUTFILE)_hip_omp.so

obj/mem_util_omp.o: src/profile_util.h src/mem_util.cpp
	mkdir -p obj;
	$(HIPCXX) $(CXXFLAGS) $(OMPFLAGS) -fPIC -std=c++14 -c src/mem_util.cpp -o obj/mem_util_omp.o

obj/time_util_omp.o: src/profile_util.h src/time_util.cpp
	mkdir -p obj;
	$(HIPCXX) $(CXXFLAGS) $(OMPFLAGS) -fPIC -std=c++14 -c src/time_util.cpp -o obj/time_util_omp.o

obj/thread_affinity_util_omp.o: src/profile_util.h src/thread_affinity_util.cpp
	mkdir -p obj;
	$(HIPCXX) $(CXXFLAGS) $(OMPFLAGS) -fPIC -std=c++14 -c src/thread_affinity_util.cpp -o obj/thread_affinity_util_omp.o


lib/$(OUTPUTFILE)_hip_mpi.so: obj/mem_util_mpi.o obj/time_util_mpi.o obj/thread_affinity_util_mpi.o
	mkdir -p lib; 
	$(MPIHIPCXX) $(MPIFLAGS) -shared obj/mem_util_mpi.o obj/time_util_mpi.o obj/thread_affinity_util_mpi.o -o lib/$(OUTPUTFILE)_hip_mpi.so

obj/mem_util_mpi.o: src/profile_util.h src/mem_util.cpp
	mkdir -p obj;
	$(MPIHIPCXX) $(MPIFLAGS) $(CXXFLAGS) -fPIC -std=c++14 -c src/mem_util.cpp -o obj/mem_util_mpi.o

obj/time_util_mpi.o: src/profile_util.h src/time_util.cpp
	mkdir -p obj;
	$(MPIHIPCXX) $(MPIFLAGS) $(CXXFLAGS)  -fPIC -std=c++14 -c src/time_util.cpp -o obj/time_util_mpi.o

obj/thread_affinity_util_mpi.o: src/profile_util.h src/thread_affinity_util.cpp
	mkdir -p obj;
	$(MPIHIPCXX) $(MPIFLAGS) $(CXXFLAGS) -fPIC -std=c++14 -c src/thread_affinity_util.cpp -o obj/thread_affinity_util_mpi.o

lib/$(OUTPUTFILE)_hip_mpi_omp.so: obj/mem_util_mpi_omp.o obj/time_util_mpi_omp.o obj/thread_affinity_util_mpi_omp.o
	mkdir -p lib; 
	$(MPIHIPCXX) $(MPIFLAGS) $(OMPFLAGS) -shared obj/mem_util_mpi_omp.o obj/time_util_mpi_omp.o obj/thread_affinity_util_mpi_omp.o -o lib/$(OUTPUTFILE)_hip_mpi_omp.so

obj/mem_util_mpi_omp.o: src/profile_util.h src/mem_util.cpp
	mkdir -p obj;
	$(MPIHIPCXX) $(CXXFLAGS) $(MPIFLAGS) $(OMPFLAGS) -fPIC -std=c++14 -c src/mem_util.cpp -o obj/mem_util_mpi_omp.o

obj/time_util_mpi_omp.o: src/profile_util.h src/time_util.cpp
	mkdir -p obj;
	$(MPIHIPCXX) $(CXXFLAGS) $(MPIFLAGS) $(OMPFLAGS) -fPIC -std=c++14 -c src/time_util.cpp -o obj/time_util_mpi_omp.o

obj/thread_affinity_util_mpi_omp.o: src/profile_util.h src/thread_affinity_util.cpp
	mkdir -p obj;
	$(MPIHIPCXX) $(CXXFLAGS) $(MPIFLAGS) $(OMPFLAGS) -fPIC -std=c++14 -c src/thread_affinity_util.cpp -o obj/thread_affinity_util_mpi_omp.o
