/*! \file thread_affinity_util.cpp
 *  \brief Get thread to core affinity
 */

#include "profile_util.h"

namespace profiling_util {

    /*
    Code to facilitate core binding reporting
    Borrowed from VELOCIraptor, which itself
    borrowed from util-linux-2.13-pre7/schedutils/taskset.c
    */
    #ifdef __APPLE__

    static inline void
    CPU_ZERO(cpu_set_t *cs) { cs->count = 0; }

    static inline void
    CPU_SET(int num, cpu_set_t *cs) { cs->count |= (1 << num); }

    static inline int
    CPU_ISSET(int num, cpu_set_t *cs) { return (cs->count & (1 << num)); }

    int sched_getaffinity(pid_t pid, size_t cpu_size, cpu_set_t *cpu_set)
    {
        int32_t core_count = 0;
        size_t  len = sizeof(core_count);
        int ret = sysctlbyname(SYSCTL_CORE_COUNT, &core_count, &len, 0, 0);
        if (ret) {
            printf("error while get core count %d\n", ret);
            return -1;
        }
        cpu_set->count = 0;
        for (int i = 0; i < core_count; i++) cpu_set->count |= (1 << i);
        return 0;
    }
    #endif

    void cpuset_to_cstr(cpu_set_t *mask, char *str)
    {
        char *ptr = str;
        int i, j, entry_made = 0;
        for (i = 0; i < CPU_SETSIZE; i++) {
            if (CPU_ISSET(i, mask)) {
                int run = 0;
                entry_made = 1;
                for (j = i + 1; j < CPU_SETSIZE; j++) {
                    if (CPU_ISSET(j, mask)) run++;
                    else break;
                }
                if (!run) {
                    sprintf(ptr, "%d ", i);
                }
                else if (run == 1) {
                    sprintf(ptr, "%d,%d ", i, i + 1);
                    i++;
                } else {
                    sprintf(ptr, "%d-%d ", i, i + run);
                    i += run;
                }
                while (*ptr != 0) ptr++;
            }
        }
        ptr -= entry_made;
        ptr = nullptr;
    }

    std::string ReportThreadAffinity()
    {
        std::string binding_report;
        // if there is no MPI and no OMP do not report any binding
#if !defined(_MPI) && !defined(_OPENMP)
        binding_report = "Serial code, binding uninformative \n ";
        return binding_report;
#endif

        int ThisTask=0, NProcs=1;
#ifdef _MPI
        MPI_Comm_size(MPI_COMM_WORLD, &NProcs);
        MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
#endif
        binding_report = "Core binding \n ";
        cpu_set_t coremask;
        char clbuf[7 * CPU_SETSIZE], hnbuf[64];
        memset(clbuf, 0, sizeof(clbuf));
        memset(hnbuf, 0, sizeof(hnbuf));
        (void)gethostname(hnbuf, sizeof(hnbuf));
    #ifdef _OPENMP
        #pragma omp parallel shared (binding_report) private(coremask, clbuf)
    #endif
        {
            std::string result;
            (void)sched_getaffinity(0, sizeof(coremask), &coremask);
            cpuset_to_cstr(&coremask, clbuf);
            result = "\t On node " + std::string(hnbuf) + " : ";
    #ifdef _MPI
            result += "MPI Rank " + std::to_string(ThisTask + " : ");
    #endif
    #ifdef _OPENMP
            auto thread = omp_get_thread_num();
            result +=" OMP Thread " + std::to_string(thread) + " : ";
    #endif
            result += " Core affinity = " + std::string(clbuf) + " \n ";
    #ifdef _OPENMP
            #pragma omp critical
    #endif
            {
                binding_report += result;
            }
        }
        return binding_report;
    }
}
