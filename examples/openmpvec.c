#include <string.h>
#include <profile_util.h>

#ifdef USEOPENMP
#include <omp.h>
#endif

int main() {
    log_parallel_api();
}
