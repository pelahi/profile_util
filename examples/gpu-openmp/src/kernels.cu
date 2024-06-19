#include "kernels.h"
#include "common.h"

template <typename T> __global__ 
void vector_square(const T *A_d, T *C_d, size_t N)
{
    size_t offset = (blockIdx.x * blockDim.x + threadIdx.x);
    size_t stride = blockDim.x * gridDim.x;

    for (size_t i=offset; i<N; i+=stride) {
        C_d[i] = A_d[i] * A_d[i];
    }
}

void silly_template_instansiation()
{
    int ix,iy;
    float fx,fy;
    double dx,dy;
    size_t N;
    size_t blocksize = 256;
    size_t threadsperblock = 1024;
    pu_gpuLaunchKernel(vector_square, 
        dim3(blocksize), dim3(threadsperblock), 
        0, 0,
        &ix, &iy, N);
    pu_gpuLaunchKernel(vector_square, 
        dim3(blocksize), dim3(threadsperblock), 
        0, 0,
        &fx, &fy, N);
    pu_gpuLaunchKernel(vector_square, 
        dim3(blocksize), dim3(threadsperblock), 
        0, 0,
        &dx, &dy, N);
}

void compute_kernel1(size_t N, 
    std::vector<int*> &x_int_gpu, 
    std::vector<int*> &y_int_gpu, 
    std::vector<float*> &x_float_gpu, 
    std::vector<float*> &y_float_gpu, 
    std::vector<double*> &x_double_gpu, 
    std::vector<double*> &y_double_gpu,
    size_t blocksize,
    size_t threadsperblock,
    int Niter
    ) 
{
    int nDevices;
    size_t dynsharedsize = 0;
    pu_gpuStream_t stream = 0;
    pu_gpuErrorCheck(pu_gpuGetDeviceCount(&nDevices));
    for (auto idev=0;idev<nDevices;idev++) {
        Log()<<" at device "<<idev<<" and doing stuff "<<std::endl;
        pu_gpuErrorCheck(pu_gpuSetDevice(idev));
        auto time_kernel = NewTimer();
        for (auto i=0; i<Niter;i++) {
            pu_gpuLaunchKernel(vector_square, 
                dim3(blocksize), dim3(threadsperblock), 
                dynsharedsize, stream,
                x_int_gpu[idev], y_int_gpu[idev], N);
            pu_gpuLaunchKernel(vector_square, 
                dim3(blocksize), dim3(threadsperblock), 
                dynsharedsize, stream,
                x_float_gpu[idev], y_float_gpu[idev], N);
            pu_gpuLaunchKernel(vector_square, 
                dim3(blocksize), dim3(threadsperblock), 
                dynsharedsize, stream,
                x_double_gpu[idev], y_double_gpu[idev], N);
        }
        LogTimeTakenOnDevice(time_kernel);
    }
}

// __global__ void vector_square(const int *, int *, size_t);
// __global__ void vector_square(const float *, float *, size_t);
// __global__ void vector_square(const double *, double *, size_t);
