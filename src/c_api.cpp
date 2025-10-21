/*! \file mem_util.cpp
 *  \brief Get memory 
 */


#include "profile_util_api_c.h"
#ifdef PU_ENABLE_C_API
// the C API
#define MYLIB_C_API __attribute__((visibility("default")))

// extern "C" wrapper function
extern "C" MYLIB_C_API void pu_get_version(){
    printf("%s\n", profiling_util::__version());
};
#endif