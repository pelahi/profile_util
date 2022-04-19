.PHONY: all clean serial openmp mpi mpi_openmp 

OUTPUTFILE=libprofile_util
OMPFLAGS=-fopenmp
MPICXX?=$(CXX)

#all: lib/$(OUTPUTFILE).so lib/$(OUTPUTFILE)_omp.so lib/$(OUTPUTFILE)_mpi.so lib/$(OUTPUTFILE)_mpi_omp.so
all: serial openmp mpi mpi_openmp
#lib/$(OUTPUTFILE).a

clean:
	rm -f obj/*o
	rm -f lib/*

serial: lib/$(OUTPUTFILE).so
	echo "Making serial library"
openmp: lib/$(OUTPUTFILE)_omp.so
	echo "Making OpenMP library"
mpi: lib/$(OUTPUTFILE)_mpi.so
	echo "Making MPI library"
mpi_openmp: lib/$(OUTPUTFILE)_mpi_omp.so
	echo "Making MPI+OpenMP library"

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
