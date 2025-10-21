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
    MYLIB_C_API void pu_report_parallel_api();
    MYLIB_C_API void pu_report_thread_affinity(const char *func, const char *file, int line);
    MYLIB_C_API void pu_report_mem_usage(const char *func, const char *file, int line);
    MYLIB_C_API void pu_report_system_mem(const char *func, const char *file, int line);
#ifdef __cplusplus
}
#endif 
#endif

#endif
