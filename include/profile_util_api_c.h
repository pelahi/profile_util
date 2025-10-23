/*! \file profile_util_api_c.h
 *  \brief this file contains all definitions that provide useful C API to library calls.
 */


#ifndef _PROFILE_UTIL_C_API
#define _PROFILE_UTIL_C_API

#ifdef ENABLE_C_API
#define MYLIB_C_API __attribute__((visibility("default")))

#ifdef __cplusplus
extern "C" {
#endif

    // Declare the function prototype with extern "C"
    MYLIB_C_API void pu_get_version();
    MYLIB_C_API void pu_report_parallel_api(const char *func, const char *file, int line);
    MYLIB_C_API void pu_report_binding(const char *func, const char *file, int line);
    MYLIB_C_API void pu_report_thread_affinity(const char *func, const char *file, int line);
    MYLIB_C_API void pu_report_mem_usage(const char *func, const char *file, int line);
    MYLIB_C_API void pu_report_system_mem(const char *func, const char *file, int line);
#ifdef _MPI
    MYLIB_C_API void pu_report_node_system_mem(const char *func, const char *file, int line);
#endif
    struct Timer_c {
        unsigned long int t0;
        unsigned long int tref;
        char *ref;
        int use_device;
#if defined(_GPU)
        void* t0_event;
        int device_id;
        int other_device_id;
        bool swap_device;
#endif
    };
    struct Timer_c * Timer_c_create(const char *func, const char *file, int line, int _use_device);
    void Timer_c_destroy(struct Timer_c * v );

    MYLIB_C_API void pu_report_time_taken(struct Timer_c *timer, const char *func, const char *file, int line);
#ifdef _GPU
    MYLIB_C_API void pu_report_time_taken_on_device(struct Timer_c *timer, const char *func, const char *file, int line);
#endif

    struct ComputeSampler_c {
        int nDevices;
        char *cpu_energy_fname, *cpu_usage_fname, *cpu_freq_fname;
#ifdef _GPU
        char *gpu_energy_fname, gpu_usage_fname, gpu_mem_fname, gpu_memusage_fname;
#endif
    };
    struct ComputeSampler_c * ComputeSampler_c_create(const char *func, const char *file, int line, float samples_per_sec, int _use_device);
    void ComputeSampler_c_destroy(struct ComputeSampler_c * v );
    void ComputeSampler_c_keepfiles(struct ComputeSampler_c * v , int _keep_files);

    MYLIB_C_API void pu_report_cpu_usage(struct ComputeSampler_c *sampler, const char *func, const char *file, int line);
#ifdef _GPU
    MYLIB_C_API void pu_report_gpu_usage(struct ComputeSampler_c *sampler, const char *func, const char *file, int line);
    MYLIB_C_API void pu_report_gpu_statistics(struct ComputeSampler_c *sampler, const char *func, const char *file, int line);
#endif
#ifdef __cplusplus
}
#endif 
#endif

#define LogParallelAPI_c(func) pu_report_parallel_api(func, __FILE__, __LINE__);
#define LogBinding_c(func) pu_report_binding(func, __FILE__, __LINE__);
#define LogThreadAffinity_c(func) pu_report_thread_affinity(func, __FILE__, __LINE__);
#define LogMemUsage_c(func) pu_report_mem_usage(func, __FILE__, __LINE__);
#define LogSystemMem_c(func) pu_report_system_mem(func, __FILE__, __LINE__);
#ifdef _MPI
#define LogNodeSystemMem_c(func) pu_report_node_system_mem(func, __FILE__, __LINE__); 
#endif

#define NewTimer_c(func) Timer_c_create(func, __FILE__, __LINE__, 1);
#define NewComputeSampler_c(func, sampling) ComputeSampler_c_create(func, __FILE__, __LINE__, sampling, 1);

#define LogTimeTaken_c(timer, func) pu_report_time_taken(timer, func, __FILE__, __LINE__);
#ifdef _GPU
#define LogTimeTakenOnDevice_c(timer, func) pu_report_time_taken_on_device(timer, func, __FILE__, __LINE__);
#endif
#define LogCPUUsage_c(sampler, func) pu_report_cpu_usage(sampler, func, __FILE__, __LINE__);
#ifdef _GPU
#define LogGPUUsage_c(sampler, func) pu_report_gpu_usage(sampler, func, __FILE__, __LINE__); 
#define LogGPUStatistics_c(sampler, func) pu_report_gpu_statistics(sampler, func, __FILE__, __LINE__); 
#endif 



#endif
