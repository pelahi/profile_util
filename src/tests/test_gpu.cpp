/*!
    \file test_gpu.cpp
    \brief Test GPU profiling and computation using the profiling utility library.
    \details This test initializes the profiling utility, performs vector addition on the GPU while logging various metrics, and verifies the results.
*/

#include "profile_util_gpu.h"

#include <vector>
#include <random>
#include <profile_util.h>

__global__ void vec_add(const int N, const float * __restrict__ din, float * __restrict__ dout){
    // Calculate the unique global thread index
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    if (index >= N) return;
    dout[index] = din[index]+din[index]*index;
}

int main(int argc, char *argv[])
{
    int ThisTask = 0;
#ifdef _MPI
    auto comm = MPI_COMM_WORLD;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(comm, &ThisTask);
    MPISetLoggingComm(comm);
#endif 
#ifdef _MPI
    MPILog0Version();
    MPILog0ParallelAPI();
    MPILog0Binding();
    MPILog0NodeSystemMem();
    MPI_Barrier(comm);
#else 
    LogVersion();
    LogParallelAPI();
    LogBinding();
    LogSystemMem();
#endif
    LogMemUsage();
    auto t1 = NewTimer();
    auto s1 = NewComputeSampler(0.01);
#ifdef _GPU
    int Nentries = 100000000;
    int blocksize = 1024;
    int threadsperblock = 256;
    size_t shsize = 0;
    pu_gpuStream_t stream = 0;
    std::vector<float> xvec(Nentries), yvec(Nentries);
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 gen(rd()); // seed the generator
    std::normal_distribution<> distr(0,1);
    float sum = 0;
    for (auto &x:xvec) {
        x = distr(gen);
        sum += x;
    }
    Log()<<Nentries<<" entries with initial sum "<<sum<<std::endl;
    int nDevices;
    std::vector<float*> xdev, ydev;
    pu_gpuErrorCheck(pu_gpuGetDeviceCount(&nDevices));
    xdev.resize(nDevices);
    ydev.resize(nDevices);
    auto tgpumem = NewTimerHostOnly();
    for (auto idev=0;idev<nDevices;idev++) {
        pu_gpuErrorCheck(pu_gpuSetDevice(idev));
        auto time_local = NewTimer();
        auto sgpu = NewComputeSampler(0.01);
        auto nbytes = Nentries * sizeof(float);
        pu_gpuErrorCheck(pu_gpuMalloc(&xdev[idev], nbytes));
        pu_gpuErrorCheck(pu_gpuMalloc(&ydev[idev], nbytes));
        pu_gpuErrorCheck(pu_gpuMemcpy(xdev[idev], xvec.data(), nbytes, pu_gpuMemcpyHostToDevice));
        pu_gpuLaunchKernel(vec_add,
            dim3(blocksize), dim3(threadsperblock), 
            shsize, stream,
            Nentries, xdev[idev], ydev[idev]);
        pu_gpuErrorCheck(pu_gpuMemcpy(yvec.data(), ydev[idev], nbytes, pu_gpuMemcpyDeviceToHost));
        pu_gpuErrorCheck(pu_gpuDeviceSynchronize());
        LogTimeTakenOnDevice(time_local);
        LogGPUStatistics(sgpu);
        sum = 0;
        for (auto &x:yvec) sum += x;
        Log()<<Nentries<<" entries with final sum "<<sum<<" on gpu "<<idev<<std::endl;
    }
    LogTimeTaken(tgpumem);
#endif

#ifdef _MPI
    MPI_Finalize();
#endif 

}
