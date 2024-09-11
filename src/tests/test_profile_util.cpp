#include <array>
#include <iostream>
#include <numeric>
#include <vector>
#include <profile_util.h>

int main(int argc, char *argv[])
{
    int NProcs = 1;
    int ThisTask = 0;
#ifdef _MPI
    auto comm = MPI_COMM_WORLD;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(comm, &ThisTask);
    MPI_Comm_size(comm, &NProcs);
    MPISetLoggingComm(comm);
#endif 
#ifdef _MPI
    MPILog0Version();
    MPILog0ParallelAPI();
    MPILog0NodeSystemMem();
    MPI_Barrier(comm);
#else 
    LogVersion();
    LogParallelAPI();
    LogSystemMem();
#endif
    LogMemUsage();

#ifdef _MPI
    MPI_Finalize();
#endif 

}
