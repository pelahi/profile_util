/*! \file mem_util.cpp
 *  \brief Get memory 
 */


#include "profile_util_api_c.h"
#include "profile_util.h"
#include <cstdio>
#ifdef ENABLE_C_API
// the C API
#define MYLIB_C_API __attribute__((visibility("default")))

// extern "C" wrapper function
extern "C" {
    MYLIB_C_API void pu_get_version(){
        printf("%s\n", profiling_util::__version().c_str());
    };
    MYLIB_C_API void pu_report_parallel_api(){
        printf("%s\n", profiling_util::ReportParallelAPI().c_str());
    };
    MYLIB_C_API void pu_report_binding(){
        printf("%s\n", profiling_util::ReportBinding().c_str());
    };
    MYLIB_C_API void pu_report_thread_affinity(const char *func, const char *file, int line){
        printf("%s\n", profiling_util::ReportThreadAffinity(func, file, std::to_string(line)).c_str());
    };
    MYLIB_C_API void pu_report_mem_usage(const char *func, const char *file, int line){
        printf("%s\n", profiling_util::ReportMemUsage(func, file, std::to_string(line)).c_str());
    };
    MYLIB_C_API void pu_report_system_mem(const char *func, const char *file, int line){
        printf("%s\n", profiling_util::ReportSystemMem(func, file, std::to_string(line)).c_str());
    };

    struct Timer_c * Timer_c_create(const char *func, const char *file, int line, int _use_device)
    {
        return reinterpret_cast<struct Timer_c*>(new profiling_util::Timer(func, file, std::to_string(line), _use_device));
    }
    void Timer_c_destroy(struct Timer_c * v )
    {
        delete reinterpret_cast<profiling_util::Timer*>(v);
    }
    
    MYLIB_C_API void pu_report_time_taken(struct Timer_c *timer, const char *func, const char *file, int line){
        printf("%s\n", profiling_util::ReportTimeTaken(*(reinterpret_cast<profiling_util::Timer*>(timer)), func, file, std::to_string(line)).c_str());
    };
#ifdef _GPU
    MYLIB_C_API void pu_report_time_taken_on_device(struct Tiner_c *timer, const char *func, const char *file, int line){
        printf("%s\n", profiling_util::ReportTimeTakenOnDevice(*(reinterpret_cast<profiling_util::Timer*>(timer)), func, file, std::to_string(line)).c_str());
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
        printf("%s\n", profiling_util::ReportCPUUsage(*(reinterpret_cast<profiling_util::ComputeSampler*>(sampler)), func, file, std::to_string(line)).c_str());
    };
// extern "C" MYLIB_C_API void pu_report_cpu_usage(const char *sampler, const char *func, const char *file, int line){
//     printf("%s\n", profiling_util::ReportCPUUsage(sampler, func, file, std::to_string(line)).c_str());
// };
// extern "C" MYLIB_C_API void pu_report_gpu_usage(const char *sampler, const char *func, const char *file, int line){
//     printf("%s\n", profiling_util::ReportGPUUsage(sampler, func, file, std::to_string(line)).c_str());
// };

}

#endif