/*!
    \file test_gpu.cpp
    \brief Test GPU profiling and computation using the profiling utility library.
    \details This test initializes the profiling utility, performs vector addition on the GPU while logging various metrics, and verifies the results.
*/

#include "profile_util_gpu.h"

#include <vector>
#include <random>
#include <string>
#include <sstream>
#include <profile_util.h>
#include <getopt.h>


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

int ThisTask=0;
#define Rank0Log() if (ThisTask==0) Log()


struct Options
{
    int Nentries = 100000000;
    int Niter = 100;
    int blocksize = 1024;
    int threadsperblock = 256;
    int kernel_type = 1;
    std::vector<int> device_list = {-1};
};

/// usage
void usage()
{
    Rank0Log()<<"Options: "<<std::endl;
    Rank0Log()<<"  -n <number of entries> (default 100000000) "<<std::endl;
    Rank0Log()<<"  -i <number of iterations for each test> (default 100) "<<std::endl;
    Rank0Log()<<"  -b <gpu block size> (default 1024) "<<std::endl;
    Rank0Log()<<"  -t <threads per block> (default 256) "<<std::endl;
    Rank0Log()<<"  -k <kernel type> (default 1 [vec_add]) "<<std::endl;
    Rank0Log()<<"     kernel type 1: simple vec add "<<std::endl;
    Rank0Log()<<"     kernel type 2: vec add with inner loop to increase runtime"<<std::endl;
    Rank0Log()<<"  -d <comma separated list of gpu devices to use> (default -1 [all devices]) "<<std::endl;

#ifdef _MPI
    MPI_Finalize();
#endif
    exit(0);
}

///routine to get arguments from command line
void GetArgs(int argc, char *argv[], Options &opt)
{
    int option;
    while ((option = getopt(argc, argv, ":n:i:b:t:k:d:")) != EOF)
    {
        switch(option)
        {
            case 'n':
                opt.Nentries = atoll(optarg);
                break;
            case 'i':
                opt.Niter = atoi(optarg);
                break;
            case 'b':
                opt.blocksize = atoi(optarg);
                break;
            case 't':
                opt.threadsperblock = atoi(optarg);
                break;
            case 'k':
                opt.kernel_type = atoi(optarg);
                break;
            case 'd':
            {
                opt.device_list.clear();
                // Split by comma delimiter
                std::stringstream ss(optarg);
                std::string token;
                while (std::getline(ss, token, ',')) {
                    opt.device_list.push_back(std::stoi(token));
                }
                break;
            }
            case '?':
                usage();
                break;
        }
    }
}


int main(int argc, char *argv[])
{
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
    Options opt;
    GetArgs(argc, argv, opt);
    int nDevices;
    std::vector<float*> xdev, ydev;
    pu_gpuErrorCheck(pu_gpuGetDeviceCount(&nDevices));
    xdev.resize(nDevices);
    ydev.resize(nDevices);
    Rank0Log() << "Running GPU tests"<< std::endl;
    Log()<<"Found "<<nDevices<<" gpu devices "<<std::endl;
    Log()<<"Running "<<opt.Nentries<<" computations for "<<opt.Niter<<" iterations on gpu devices "<<std::endl;
    if (opt.device_list.size()==1 && opt.device_list[0]==-1) {
        opt.device_list.clear();
        for (int i=0;i<nDevices;i++) opt.device_list.push_back(i);
    }
    else {
        for (auto &idev:opt.device_list) {
            if (idev>=nDevices || idev<0) {
                Rank0Log()<<"Error: device "<<idev<<" is out of range. Found "<<nDevices<<" devices. Exiting."<<std::endl;
#ifdef _MPI
                MPI_Finalize();
#endif
                exit(1);
            }
        }
    }
    std::string s = "Running on devices : ";
    for (auto &idev:opt.device_list) s+=std::to_string(idev)+", ";
    Log()<<s<<std::endl;

    Log()<<"Starting time and compute samplers"<<std::endl;
    auto t1 = NewTimer();
    auto s1 = NewComputeSampler(0.01);
    size_t shsize = 0;
    pu_gpuStream_t stream = 0;
    std::vector<float> xvec(opt.Nentries), yvec(opt.Nentries);
    std::random_device rd; // obtain a random number from hardware
    std::mt19937 gen(rd()); // seed the generator
    std::normal_distribution<> distr(0,1);
    float sum = 0;
    Log()<<"Starting random number generation for "<<opt.Nentries<<std::endl;
    // CUDA too old for vector auto iteration parallelisation
#if defined(_OPENMP) && !defined(_CUDA) 
#pragma omp parallel for reduction(+:sum) schedule(static)
#endif
    for (auto &x:xvec) {
        x = distr(gen);
        sum += x;
    }

    Log()<<opt.Nentries<<" entries with initial sum "<<sum<<std::endl;
    Log()<<"Starting GPU computation"<<std::endl;
    auto tgpumem = NewTimerHostOnly();
    for (auto &idev:opt.device_list) 
    {
        pu_gpuErrorCheck(pu_gpuSetDevice(idev));
        auto time_local = NewTimer();
        auto sgpu = NewComputeSampler(0.01);
        auto nbytes = opt.Nentries * sizeof(float);
        pu_gpuErrorCheck(pu_gpuMalloc(&xdev[idev], nbytes));
        pu_gpuErrorCheck(pu_gpuMalloc(&ydev[idev], nbytes));
        pu_gpuErrorCheck(pu_gpuMemcpy(xdev[idev], xvec.data(), nbytes, pu_gpuMemcpyHostToDevice));
        for (auto iter=0;iter<opt.Niter;iter++) {
            if (opt.kernel_type == 1) {
                pu_gpuLaunchKernel(vec_add,
                    dim3(opt.blocksize), dim3(opt.threadsperblock), 
                    shsize, stream,
                    opt.Nentries, xdev[idev], ydev[idev]);
            } 
            else {
                pu_gpuLaunchKernel(vec_add_silly,
                    dim3(opt.blocksize), dim3(opt.threadsperblock), 
                    shsize, stream,
                    opt.Nentries, opt.Niter, xdev[idev], ydev[idev]);
            }
        }
        pu_gpuErrorCheck(pu_gpuMemcpy(yvec.data(), ydev[idev], nbytes, pu_gpuMemcpyDeviceToHost));
        pu_gpuErrorCheck(pu_gpuDeviceSynchronize());
        LogTimeTakenOnDevice(time_local);
        LogGPUStatistics(sgpu);
        sum = 0;
        for (auto &x:yvec) sum += x;
        Log()<<opt.Nentries<<" entries with final sum "<<sum<<" on gpu "<<idev<<std::endl;
    }
    LogTimeTaken(tgpumem);

#ifdef _MPI
    MPI_Finalize();
#endif 

}