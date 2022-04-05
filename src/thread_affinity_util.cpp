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

    std::string MPICallingRank(int task){
        char s[20];
        sprintf(s,"MPI [%04d]: ",task);
        return std::string(s);
    }

    std::string ReportParallelAPI() 
    {
        std::string s;
        s = "Parallel API's \n ======== \n";
#ifdef _MPI
        int rank, size;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        s += "MPI Comm world size " + std::to_string(size);
        s += "\n";
#endif 
#ifdef _OPENMP 
        s += "OpenMP version " + std::to_string(_OPENMP);
        s += " with total number of threads = " + std::to_string(omp_get_max_threads());
        s += "\n";
#endif
        return s;
    }

    std::string ReportBinding()
    {
        std::string binding_report;
        int ThisTask=0, NProcs=1;
#ifdef _MPI
        MPI_Comm_size(MPI_COMM_WORLD, &NProcs);
        MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
#endif
        binding_report = "Core Binding \n ======== \n";
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
            result += "MPI Rank " + std::to_string(ThisTask) + " : ";
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
    /// return binding as called within openmp region 
    std::string ReportThreadAffinity(std::string func, std::string line)
    {
        std::string result;
        cpu_set_t coremask;
        char clbuf[7 * CPU_SETSIZE], hnbuf[64];
        memset(clbuf, 0, sizeof(clbuf));
        memset(hnbuf, 0, sizeof(hnbuf));
        (void)gethostname(hnbuf, sizeof(hnbuf));
        result = "Thread affinity report @ " + func + " L" + line + " : ";
        (void)sched_getaffinity(0, sizeof(coremask), &coremask);
        cpuset_to_cstr(&coremask, clbuf);
        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif
        result += " Thread " + std::to_string(thread) + " : ";
        result += " Core affinity = " + std::string(clbuf) + " ";
        result += " Core placement = " + std::to_string(sched_getcpu()) + " ";
        result += "\n";

        return result;
    }

    /// return binding as called within openmp region, MPI aware 
#ifdef _MPI 
    std::string MPIReportThreadAffinity(std::string func, std::string line, MPI_Comm &comm)
    {
        std::string result;
        int ThisTask=0, NProcs=1;
        cpu_set_t coremask;
        char clbuf[7 * CPU_SETSIZE], hnbuf[64];

        MPI_Comm_size(comm, &NProcs);
        MPI_Comm_rank(comm, &ThisTask);
        memset(hnbuf, 0, sizeof(hnbuf));
        memset(clbuf, 0, sizeof(clbuf));
        (void)gethostname(hnbuf, sizeof(hnbuf));
        result = "Thread affinity report @ " + func + " L" + line + " : ";
        result += "::\t On node " + std::string(hnbuf) + " : ";
        result += "MPI Rank " + std::to_string(ThisTask) + " : ";
        (void)sched_getaffinity(0, sizeof(coremask), &coremask);
        cpuset_to_cstr(&coremask, clbuf);
        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif
        result += " Thread " + std::to_string(thread) + " : ";
        result += " Core affinity = " + std::string(clbuf) + " \n ";
        result += " Core placement = " + std::to_string(sched_getcpu()) + " ";
        result += "\n";

        return result;
    }
#endif

}

extern "C" {
    void report_binding(char c[])
    {        
        strcpy(c,profiling_util::ReportBinding().c_str());
    }
    void report_thread_affinity(char c[], char *f, int l)
    {        
        strcpy(c,profiling_util::ReportThreadAffinity(std::string(f), std::to_string(l)).c_str());
    }
}
