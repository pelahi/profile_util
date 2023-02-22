#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <algorithm>
#include <random>
#include <tuple>
#include <thread>
#include <mpi.h>
#include <profile_util.h>

int ThisTask, NProcs;
std::chrono::system_clock::time_point logtime;
std::time_t log_time;
char wherebuff[1000];
std::string whenbuff;

#define Where() sprintf(wherebuff,"[%04d] @%sL%d ", ThisTask,__func__, __LINE__);
#define When() logtime = std::chrono::system_clock::now(); log_time = std::chrono::system_clock::to_time_t(logtime);whenbuff=std::ctime(&log_time);whenbuff.erase(std::find(whenbuff.begin(), whenbuff.end(), '\n'), whenbuff.end());
#define LocalLogger() Where();std::cout<<wherebuff<<" : " 
#define Rank0LocalLogger() Where();if (ThisTask==0) std::cout<<wherebuff<<" : " 
#define LocalLoggerWithTime() Where();When(); std::cout<<wherebuff<<" ("<<whenbuff<<") : "
#define Rank0LocalLoggerWithTime() Where();When(); if (ThisTask==0) std::cout<<wherebuff<<" ("<<whenbuff<<") : "
#define LogMPITest() Rank0LocalLoggerWithTime()<<" running "<<mpifunc<< " test"<<std::endl;
#define LogMPIBroadcaster() if (ThisTask == itask) LocalLoggerWithTime()<<" running "<<mpifunc<<" broadcasting "<<sendsize<<" GB"<<std::endl;
#define LogMPISender() LocalLoggerWithTime()<<" Running "<<mpifunc<<" sending "<<sendsize<<" GB"<<std::endl;
#define LogMPIReceiver() if (ThisTask == itask) LocalLoggerWithTime()<<" running "<<mpifunc<<std::endl;
#define LogMPIAllComm() Rank0LocalLoggerWithTime()<<" running "<<mpifunc<<" all "<<sendsize<<" GB"<<std::endl;
#define Rank0ReportMem() if (ThisTask==0) {Where();When();std::cout<<wherebuff<<" ("<<whenbuff<<") : ";LogMemUsage();std::cout<<wherebuff<<" ("<<whenbuff<<") : ";LogSystemMem();}

void WriteCollective(MPI_Comm &comm) {
    Rank0LocalLoggerWithTime()<<" Starting collective binary write "<<std::endl;    
    MPI_File file;
    MPI_Offset offset;
    MPI_Status status;
    std::vector<int> buffer(100);
    for (auto i=0;i<buffer.size();i++) buffer[i] = ThisTask*NProcs + i;

    // Open the file in parallel
    MPI_File_open(comm, "mpi-io.collective.example.txt", 
        MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, 
        &file);
    // Calculate the offset based on the rank
    offset = ThisTask * buffer.size()*sizeof(int);
    // Write to the file using collective I/O of binary data
    MPI_File_write_at_all(file, offset, buffer.data(), buffer.size(), MPI_INT, &status);
    // Close the file
    MPI_File_close(&file);
    Rank0LocalLoggerWithTime()<<" Finished collective binary write "<<std::endl;    
}

void WriteNonCollective(MPI_Comm &comm) {
    Rank0LocalLoggerWithTime()<<" Starting non-collective binary write "<<std::endl;    
    MPI_File file;
    MPI_Offset offset;
    MPI_Status status;
    std::vector<int> buffer(100);
    for (auto i=0;i<buffer.size();i++) buffer[i] = ThisTask*NProcs + i;

    // Open the file in parallel
    MPI_File_open(comm, "mpi-io.non-collective.example.txt", 
        MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, 
        &file);
    // Calculate the offset based on the rank
    offset = ThisTask * buffer.size()*sizeof(int);
    // Write to the file using non-collective I/O
    MPI_File_write_at(file, offset, buffer.data(), 50, MPI_INT, &status);
    // Write to the file again using non-collective I/O
    MPI_File_write_at(file, offset + sizeof(int)*50, buffer.data() + 50, 50, MPI_INT, &status);

    // Close the file
    MPI_File_close(&file);
    Rank0LocalLoggerWithTime()<<" Finished non-collective binary write "<<std::endl;    
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    auto comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &ThisTask);
    MPI_Comm_size(comm, &NProcs);

    // init logger time
    logtime = std::chrono::system_clock::now();
    log_time = std::chrono::system_clock::to_time_t(logtime);
    auto start = std::chrono::system_clock::now();
    std::time_t start_time = std::chrono::system_clock::to_time_t(start);
    Rank0LocalLoggerWithTime()<<"Starting job "<<std::endl;
    Rank0ReportMem();
    MPILog0NodeMemUsage(comm);
    MPILog0NodeSystemMem(comm);
    MPILog0ParallelAPI();
    MPILog0Binding();

    // trial some writes 
    WriteCollective(comm);
    WriteNonCollective(comm);

    // finish job
    Rank0LocalLoggerWithTime()<<"Ending job "<<std::endl;
    MPI_Finalize();
    return 0;
}
