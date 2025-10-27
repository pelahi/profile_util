#include <profile_util_api_c.h>
#include <stdio.h>
#include<unistd.h>

#ifdef _MPI
#include <mpi.h>
#endif

int main(int *argc, char ***argv) {
#ifdef _MPI
    MPI_Init(argc, argv);
#endif
    #ifdef _MPI
    MPI_Comm comm = MPI_COMM_WORLD;
    pu_mpi_set_logging_comm(comm);
    #endif
    LogParallelAPI_c("main");
    LogBinding_c("main");
    LogThreadAffinity_c("main");
    LogMemUsage_c("main");
    LogSystemMem_c("main");
    #ifdef _MPI
    LogNodeSystemMem_c("main");
    #endif
    struct Timer_c *timer = NewTimer_c("main");
    struct ComputeSampler_c *sampler = NewComputeSampler_c("main", 0.5);
    ComputeSampler_c_keepfiles(sampler, 1);
    sleep(1); // Simulate some work
    int i,j;
    double sum = 0;
    for (j=0; j<10000; j++){
        for (i=0; i<1000000; i++) {
            double x = 3.14159 * 2.71828; // Dummy computation
            sum += x;
        }
    }
    printf("%f\n", sum);
    LogTimeTaken_c(timer, "main");
    LogCPUUsage_c(sampler, "main");
    Timer_c_destroy(timer);
    ComputeSampler_c_destroy(sampler);
#ifdef _MPI
    MPI_Finalize();
#endif
    return 0;
}
