#include <profile_util_api_c.h>
#include <stdio.h>

#ifdef _MPI
#include <mpi.h>
#endif

int main(int *argc, char ***argv) {
#ifdef _MPI
  MPI_Init(argc, argv);
#endif
  pu_get_version();
  pu_report_parallel_api();
  pu_report_thread_affinity("main", __FILE__, __LINE__);
  pu_report_mem_usage("main", __FILE__, __LINE__);
  pu_report_system_mem("main", __FILE__, __LINE__);
#ifdef _MPI
  MPI_Finalize();
#endif
  return 0;
}
