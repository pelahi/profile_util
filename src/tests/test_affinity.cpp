#include <profile_util.h>

int main(int argc, char *argv[])
{
    int ThisTask = 0;
#ifdef _MPI
    auto comm = MPI_COMM_WORLD;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(comm, &ThisTask);
    MPISetLoggingComm(comm);
#endif 
#ifdef _MPI
    MPILog0ParallelAPI();
    MPILog0Binding();
#else 
    LogParallelAPI();
    LogBinding();
#endif

#ifdef _MPI
    MPI_Finalize();
#endif 

}
