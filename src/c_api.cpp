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
extern "C" MYLIB_C_API void pu_get_version(){
    printf("%s\n", profiling_util::__version().c_str());
};
extern "C" MYLIB_C_API void pu_report_parallel_api(){
    printf("%s\n", profiling_util::ReportParallelAPI().c_str());
};
extern "C" MYLIB_C_API void pu_report_thread_affinity(const char *func, const char *file, int line){
    printf("%s\n", profiling_util::ReportThreadAffinity(func, file, std::to_string(line)).c_str());
};
extern "C" MYLIB_C_API void pu_report_mem_usage(const char *func, const char *file, int line){
    printf("%s\n", profiling_util::ReportMemUsage(func, file, std::to_string(line)).c_str());
};
extern "C" MYLIB_C_API void pu_report_system_mem(const char *func, const char *file, int line){
    printf("%s\n", profiling_util::ReportSystemMem(func, file, std::to_string(line)).c_str());
};


#endif