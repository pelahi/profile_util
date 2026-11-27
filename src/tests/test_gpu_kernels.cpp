#include <cmath>
#include <profile_util_gpu.h>
#include "test_gpu_kernels.h"

__global__ void vec_add(const int N, const float * __restrict__ din, float * __restrict__ dout){
    // Calculate the unique global thread index
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    if (index >= N) return;
    dout[index] = din[index]+din[index]*index;
}

__global__ void vec_add_silly(const int N, const int Ni, const float * __restrict__ din, float * __restrict__ dout){
    // Calculate the unique global thread index
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    if (index >= N) return;
    for (int i=0;i<Ni;i++) {
        dout[index] = din[index]+din[index]*index;
    }
}

__global__ void vec_silly(const int N, const int Ni, const float * __restrict__ din, float * __restrict__ dout){
    // Calculate the unique global thread index
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    if (index >= N) return;
    for (int i=0;i<Ni;i++) {
        double delta = 0;
        if (i==0 || i==N-1) delta = din[index] ;
        else delta = din[index]-din[index+1];
        dout[index] = din[index]+din[index]*index / pow(1.1,i) + sin(din[index])+din[index]*delta;
    }
}

void launch_gpu_kernel(int kernel_type, int blocksize, int threadsperblock,
    int Nentries, int Nloop, const float *din, float *dout)
{
    size_t shsize = 0;
    pu_gpuStream_t stream = 0;
    if (kernel_type == 1) {
        pu_gpuLaunchKernel(vec_add, dim3(blocksize), dim3(threadsperblock),
            shsize, stream, Nentries, din, dout);
    }
    else if (kernel_type == 2) {
        pu_gpuLaunchKernel(vec_add_silly, dim3(blocksize), dim3(threadsperblock),
            shsize, stream, Nentries, Nloop, din, dout);
    }
    else if (kernel_type == 3) {
        pu_gpuLaunchKernel(vec_silly, dim3(blocksize), dim3(threadsperblock),
            shsize, stream, Nentries, Nloop, din, dout);
    }
}
