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
    pu_get_version();
    pu_report_parallel_api();
    pu_report_binding();
    pu_report_thread_affinity("main", __FILE__, __LINE__);
    pu_report_mem_usage("main", __FILE__, __LINE__);
    pu_report_system_mem("main", __FILE__, __LINE__);
    struct Timer_c *timer = Timer_c_create("main", __FILE__, __LINE__, 1);
    struct ComputeSampler_c *sampler = ComputeSampler_c_create("main", __FILE__, __LINE__, 0.1, 1);
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
    pu_report_time_taken(timer, "main", __FILE__, __LINE__);
    pu_report_cpu_usage(sampler, "main", __FILE__, __LINE__);
    Timer_c_destroy(timer);
    ComputeSampler_c_destroy(sampler);
#ifdef _MPI
    MPI_Finalize();
#endif
    return 0;
}
