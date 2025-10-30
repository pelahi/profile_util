/*! \file mem_util.cpp
 *  \brief Get memory 
 */


#ifdef ENABLE_C_API
#include "profile_util_api_c.h"
#include "profile_util.h"
#include <cstdio>
// the C API
#define MYLIB_C_API __attribute__((visibility("default")))

#define _where_calling_from_c(func,file,line) std::string("@")+std::string(func)+std::string(" ")+profiling_util::__extract_filename(std::string(file))+":L"+std::to_string(line)+" "
#define _when_calling_from_c "("+profiling_util::__when()+") : "
#ifdef _MPI 
#define _MPI_calling_rank_c "["+std::to_string(profiling_util::__comm_rank)+"] "
#define _log_header_c(func,file,line) _MPI_calling_rank_c+_where_calling_from_c(func,file,line)+_when_calling_from_c
#else 
#define _log_header_c(func,file,line) _where_calling_from_c(func,file,line)+_when_calling_from_c
#endif


// extern "C" wrapper function
extern "C" {
    MYLIB_C_API void pu_get_version(){
        printf("%s\n", profiling_util::__version().c_str());
    };
#ifdef _MPI
    MYLIB_C_API void pu_mpi_set_logging_comm(MPI_Comm comm){
        MPISetLoggingComm(comm);
    };
#endif
    MYLIB_C_API void pu_report_parallel_api(const char *func, const char *file, int line){
        std::string header = _log_header_c(func,file,line);
#ifdef _MPI
        if (profiling_util::__comm_rank == 0)
#endif
        printf("%s %s\n", header.c_str(), profiling_util::ReportParallelAPI().c_str());
    };
    MYLIB_C_API void pu_report_binding(const char *func, const char *file, int line){
        std::string header = _log_header_c(func,file,line);
#ifdef _MPI
        if (profiling_util::__comm_rank == 0)
#endif
        printf("%s %s\n", header.c_str(), profiling_util::ReportBinding().c_str());
    };

    MYLIB_C_API void pu_report_thread_affinity(const char *func, const char *file, int line){
        std::string header = _log_header_c(func,file,line);
        printf("%s %s\n", header.c_str(), profiling_util::ReportThreadAffinity(func, file, std::to_string(line)).c_str());
    };
    MYLIB_C_API void pu_report_mem_usage(const char *func, const char *file, int line){
        std::string header = _log_header_c(func,file,line);
        printf("%s %s\n", header.c_str(), profiling_util::ReportMemUsage(func, file, std::to_string(line)).c_str());
    };
    MYLIB_C_API void pu_report_system_mem(const char *func, const char *file, int line){
        std::string header = _log_header_c(func,file,line);
        printf("%s %s\n", header.c_str(), profiling_util::ReportSystemMem(func, file, std::to_string(line)).c_str());
    };
#ifdef _MPI
    MYLIB_C_API void pu_report_node_system_mem(const char *func, const char *file, int line){
        std::string header = _log_header_c(func,file,line);
        std::string __s = profiling_util::MPIReportNodeSystemMem(profiling_util::__comm, func, file, std::to_string(line));
        if (profiling_util::__comm_rank == 0)
            printf("%s %s\n", header.c_str(), __s.c_str());
    };
#endif

    struct Timer_c * Timer_c_create(const char *func, const char *file, int line, int _use_device)
    {
        return reinterpret_cast<struct Timer_c*>(new profiling_util::Timer(func, file, std::to_string(line), _use_device));
    }
    void Timer_c_destroy(struct Timer_c * v )
    {
        delete reinterpret_cast<profiling_util::Timer*>(v);
    }
    
    MYLIB_C_API void pu_report_time_taken(struct Timer_c *timer, const char *func, const char *file, int line){
        std::string header = _log_header_c(func,file,line);
        printf("%s %s\n", header.c_str(), profiling_util::ReportTimeTaken(*(reinterpret_cast<profiling_util::Timer*>(timer)), func, file, std::to_string(line)).c_str());
    };
#ifdef _GPU
    MYLIB_C_API void pu_report_time_taken_on_device(struct Tiner_c *timer, const char *func, const char *file, int line){
        std::string header = _log_header_c(func,file,line);
        printf("%s %s\n", header.c_str(), profiling_util::ReportTimeTakenOnDevice(*(reinterpret_cast<profiling_util::Timer*>(timer)), func, file, std::to_string(line)).c_str());
    };
#endif

    struct ComputeSampler_c * ComputeSampler_c_create(const char *func, const char *file, int line, float samples_per_sec, int _use_device)
    {
        return reinterpret_cast<struct ComputeSampler_c*>(new profiling_util::ComputeSampler(func, file, std::to_string(line), samples_per_sec, _use_device));
    }
    void ComputeSampler_c_destroy(struct ComputeSampler_c * v )
    {
        delete reinterpret_cast<profiling_util::ComputeSampler*>(v);
    }
    void ComputeSampler_c_keepfiles(struct ComputeSampler_c * v, int _keep_files)
    {
        reinterpret_cast<profiling_util::ComputeSampler*>(v)->SetKeepFiles(_keep_files);
    }

    MYLIB_C_API void pu_report_cpu_usage(struct ComputeSampler_c *sampler, const char *func, const char *file, int line){
        std::string header = _log_header_c(func,file,line);
        printf("%s %s\n", header.c_str(), profiling_util::ReportCPUUsage(*(reinterpret_cast<profiling_util::ComputeSampler*>(sampler)), func, file, std::to_string(line)).c_str());
    };
#ifdef _GPU
    MYLIB_C_API void pu_report_gpu_usage(struct ComputeSampler_c *sampler, const char *func, const char *file, int line){
        std::string header = _log_header_c(func,file,line);
        printf("%s %s\n", header.c_str(), profiling_util::ReportGPUUsage(*(reinterpret_cast<profiling_util::ComputeSampler*>(sampler)), func, file, std::to_string(line)).c_str());
    };
    MYLIB_C_API void pu_report_gpu_statistics(struct ComputeSampler_c *sampler, const char *func, const char *file, int line){
        std::string header = _log_header_c(func,file,line);
        printf("%s %s\n", header.c_str(), profiling_util::ReportGPUStatistics(*(reinterpret_cast<profiling_util::ComputeSampler*>(sampler)), func, file, std::to_string(line)).c_str());
    };
#endif

}

#endif