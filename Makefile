# Simple Makefile

OUTPUTFILEBASE=libprofile_util
OMPFLAGS=-fopenmp
MPICXX?=$(CXX)
CXXFLAGS=-fPIC -std=c++14 -O2

# different sets of devices to build for
deviceset = cpu cuda hip
# different builds 
buildset = serial openmp mpi mpi_openmp

.PHONY: all clean 
.PHONY: cpu_set serial openmp mpi mpi_openmp 
.PHONY: cuda_set cuda_serial cuda_mpi cuda_openmp cuda_mpi_openmp
.PHONY: hip_set hip_serial hip_mpi hip_openmp hip_mpi_openmp

all: cpu_set

cpu_set: serial openmp mpi mpi_openmp 
cuda_set: serial openmp mpi mpi_openmp 
hip_set: serial openmp mpi mpi_openmp 

buildname=
COMPILER=$(CXX)
COMPILEFLAGS=$(CXXFLAGS)
OBJS = obj/mem_util$(buildname).o obj/time_util$(buildname).o obj/thread_affinity_util$(buildname).o 

serial: lib/$(OUTPUTFILEBASE)$(buildname).so
	echo "Making serial library"

lib/$(OUTPUTFILEBASE)$(buildname).so: $(OBJS)
	$(CXX) -shared $(OBJS) -o lib/$(OUTPUTFILE).so

$(OBJS): obj/%.o : src/%.cpp include/profile_util.h

obj/%.o: src/%.cpp include/profile_util.h
	$(COMPILER) $(COMPILERFLAGS) -I../include/ -c $< -o $@ 


# define make_target = 
# OBJS = obj/mem_util$(DT)$(BT).o obj/time_util$(DT)$(BT).o obj/thread_affinity_util$(DT)$(BT).o 
# $(OBJS): obj/%.o : src/%.cpp
# obj/%$(DT)$(BT).o: src/%.cpp
# 	$(CXX) $(CXXFLAGS) $(COMPILEFLAGS) -c $< -o $@ 

# $(BT): lib/$(OUTPUTFILEBASE)$(DT)$(BT).so
# 	echo "Making serial library"

# lib/$(OUTPUTFILEBASE)$(DT)$(BT).so: $(OBJS)
# 	$(CXX) -shared $(OBJS) -o lib/$(OUTPUTFILE).so

# endef

# $(eval $(foreach BT,$(buildset),$(call make_target,$(BT))))

clean:
	rm -f obj/*o
	rm -f lib/*
