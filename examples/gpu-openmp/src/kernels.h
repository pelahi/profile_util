/*! \file kernels.h
 *  \brief kernels of for warming up and running on GPU
 */

#ifndef KERNELS_H
#define KERNELS_H

#include <profile_util.h>

/// GPU kernels types definitions
//@{
#define KERNEL_DEFAULT 0
//@}

/// \defgroup kernels
/// GPU kernels
//@{
template <typename T> __global__ void vector_square(const T *A_d, T *C_d, size_t N);
void silly_template_instansiation();

//@}
#endif
