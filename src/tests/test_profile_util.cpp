#include <array>
#include <iostream>
#include <numeric>
#include <vector>
#include <profile_util.h>

int main(int argc, char *argv[])
{
#ifdef _MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
    MPI_Comm_size(MPI_COMM_WORLD, &NProcs);
#else
    int NProcs = 1;
    int ThisTask = 0;
#endif 
#ifdef _MPI
    MPILog0Version();
    MPILog0ParallelAPI();
#else 
    LogVersion();
    LogParallelAPI();
#endif
    LogSystemMem();
    LogMemUsage();

#ifdef _MPI
    MPI_Finalize();
#endif 

}
