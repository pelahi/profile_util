/*! \file profile_util_api_c.h
 *  \brief this file contains all definitions that provide useful C API to library calls.
 */


#ifndef _PROFILE_UTIL_C_API
#define _PROFILE_UTIL_C_API

#ifdef PU_ENABLE_C_API
#define MYLIB_C_API __attribute__((visibility("default")))


#ifdef __cplusplus
extern "C" {
#endif
    // Declare the function prototype with extern "C"
    MYLIB_C_API void pu_get_version();
#ifdef __cplusplus
}
#endif 
#endif

#endif
