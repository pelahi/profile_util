.PHONY: all clean

OUTPUTFILE=libprofile_util.so
all: lib/$(OUTPUTFILE)

clean:
	rm -f obj/*o
	rm -f lib/*

lib/$(OUTPUTFILE): obj/mem_util.o obj/time_util.o obj/thread_affinity_util.o
	mkdir -p lib; 
	$(CXX) -shared obj/mem_util.o obj/time_util.o obj/thread_affinity_util.o -o lib/$(OUTPUTFILE)

obj/mem_util.o: src/profile_util.h src/mem_util.cpp
	mkdir -p obj;
	$(CXX) $(CXXFLAGS) -fPIC -std=c++14 -c src/mem_util.cpp -o obj/mem_util.o

obj/time_util.o: src/profile_util.h src/time_util.cpp
	mkdir -p obj;
	$(CXX) $(CXXFLAGS) -fPIC -std=c++14 -c src/time_util.cpp -o obj/time_util.o

obj/thread_affinity_util.o: src/profile_util.h src/thread_affinity_util.cpp
	mkdir -p obj;
	$(CXX) $(CXXFLAGS) -fPIC -std=c++14 -c src/thread_affinity_util.cpp -o obj/thread_affinity_util.o
