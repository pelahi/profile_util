#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <random>
#include <tuple>
#include <thread>
#include <profile_util.h>
#include <mpi.h>
#include <getopt.h>


int ThisTask, NProcs;
#define Rank0Log() if (ThisTask==0) Log()

void WriteCollective(MPI_Comm &comm, std::string &fnamebase, size_t numints=100) {
    size_t bytes_written = NProcs * numints * sizeof(int);
    std::string fname=fnamebase+std::string(".collective.example.txt");
    Rank0Log()<<" Starting collective binary write to "<<bytes_written<<" bytes to "<<fname<<std::endl;
    MPI_File file;
    MPI_Offset offset;
    MPI_Status status;
    std::vector<int> buffer(numints);
    for (auto i=0u;i<buffer.size();i++) buffer[i] = ThisTask*NProcs + i;

    // Open the file in parallel
    auto t1 = NewTimer();
    MPI_File_open(comm, fname.c_str(), 
        MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, 
        &file);
    // Calculate the offset based on the rank
    offset = ThisTask * numints*sizeof(int);
    // Write to the file using collective I/O of binary data
    MPI_File_write_at_all(file, offset, buffer.data(), numints, MPI_INT, &status);
    // Close the file
    MPI_File_close(&file);
    auto deltat = static_cast<double>((t1.get()))*1e-9; // convert to seconds
    // LogTimeTaken(t1);
    Rank0Log()<<" Finished collective binary write of "<<bytes_written<<" bytes at rate of "<<bytes_written/deltat<<" bytes/second"<<std::endl; 
}

void WriteNonCollective(MPI_Comm &comm, std::string &fnamebase, size_t numints=100) {
    size_t bytes_written = NProcs * numints * sizeof(int);
    std::string fname=fnamebase+std::string(".non-collective.example.txt");
    Rank0Log()<<" Starting non-collective binary write to "<<bytes_written<<" bytes to "<<fname<<std::endl;
    MPI_File file;
    MPI_Offset offset;
    MPI_Status status;
    std::vector<int> buffer(numints);
    for (auto i=0u;i<buffer.size();i++) buffer[i] = ThisTask*NProcs + i;

    // Open the file in parallel
    auto t1 = NewTimer();
    MPI_File_open(comm, fname.c_str(), 
        MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, 
        &file);
    // Calculate the offset based on the rank
    offset = ThisTask * numints*sizeof(int);
    // Write to the file using non-collective I/O
    MPI_File_write_at(file, offset, buffer.data(), 50, MPI_INT, &status);
    // Write to the file again using non-collective I/O
    MPI_File_write_at(file, offset + sizeof(int)*50, buffer.data() + 50, 50, MPI_INT, &status);

    // Close the file
    MPI_File_close(&file);
    auto deltat = static_cast<double>((t1.get()))*1e-9; // convert to seconds
    // LogTimeTaken(t1);
    Rank0Log()<<" Finished non-collective binary write of "<<bytes_written<<" bytes at rate of "<<bytes_written/deltat<<" bytes/second"<<std::endl; 
}


int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    auto comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &ThisTask);
    MPI_Comm_size(comm, &NProcs);
    MPISetLoggingComm(comm);
    std::string fname = "mpi-io";
    if (argc >= 2) fname = std::string(argv[1]);
    size_t numints = 100;
    if (argc >= 3) numints = std::stoul(argv[2]);

    Rank0Log()<<"Starting job "<<std::endl;
    MPILog0Version();
    MPILog0ParallelAPI();
    MPILog0Binding();
    MPILog0NodeSystemMem();
    MPI_Barrier(comm);

    // trial some writes 
    WriteCollective(comm,fname, numints);
    WriteNonCollective(comm,fname, numints);

    // finish job
    Rank0Log()<<"Ending job"<<std::endl;
    MPI_Finalize();
    return 0;
}
